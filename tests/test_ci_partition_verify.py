"""Tests for scripts/ci_partition_verify.py, the second partition-drift CI step.

The script decides which active manifests get a full ``freeze-partition
--verify`` and refuses the held-out CGGA cohort unless --allow-held-out is
given. Everything here is synthetic and lives under tmp_path; the real CGGA
path is never touched (a path merely *named* like it is created under tmp_path).
The ``gliostrata`` call itself is replaced by a recording fake.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "ci_partition_verify.py"
CI_YAML = ROOT / ".github" / "workflows" / "ci.yml"

_spec = importlib.util.spec_from_file_location("ci_partition_verify", SCRIPT)
verify = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(verify)


def _manifest(tmp_path, name, cohort):
    path = tmp_path / f"{name}.json"
    path.write_text(json.dumps({"cohort_source": str(cohort)}), encoding="utf-8")
    return path.as_posix()


class Recorder:
    def __init__(self, code=0):
        self.calls = []
        self.code = code

    def __call__(self, cmd):
        self.calls.append(cmd)
        return self.code


def _go(manifests, *, allow=False, code=0):
    out, rec = [], Recorder(code)
    status = verify.verify_active(manifests, allow_held_out=allow, say=out.append, run=rec)
    return status, "\n".join(out), rec


# --- the path test -------------------------------------------------------------
@pytest.mark.parametrize(
    "path",
    [
        "data/processed/cgga/eligibility.parquet",
        "data/raw/cgga/x.csv",
        "data\\processed\\cgga\\x.parquet",
        "DATA/Raw/CGGA/x",
        "data/raw/tcga_gbm/../cgga/x",
        "/tmp/anything/data/processed/cgga/x.parquet",
    ],
)
def test_cgga_paths_are_held_out(path):
    assert verify.is_held_out(path)


@pytest.mark.parametrize(
    "path",
    [
        "data/raw/tcga_gbm/clinical_cohort.parquet",
        "data/processed/tcga_gbm/cgga_like.parquet",
        "data/raw/cgga_other/x",
        "cgga/data/raw/x",
    ],
)
def test_other_paths_are_not_held_out(path):
    assert not verify.is_held_out(path)


# --- behaviour kept from the shell loop ----------------------------------------
def test_absent_cohort_is_skipped(tmp_path):
    m = _manifest(tmp_path, "v2", tmp_path / "data" / "raw" / "tcga_gbm" / "c.parquet")
    status, out, rec = _go([m])
    assert status == 0 and rec.calls == []
    assert "Skipping: cohort" in out and "is not present in this checkout" in out


def test_present_cohort_runs_the_same_verify_command(tmp_path):
    cohort = tmp_path / "data" / "raw" / "tcga_gbm" / "c.parquet"
    cohort.parent.mkdir(parents=True)
    cohort.write_bytes(b"")
    m = _manifest(tmp_path, "v2", cohort)
    status, out, rec = _go([m])
    assert status == 0
    assert rec.calls == [
        ["gliostrata", "freeze-partition", "--input", str(cohort), "--manifest", m, "--verify"]
    ]
    assert f"Verifying {m} against {cohort}" in out


def test_a_failing_verification_stops_with_its_status(tmp_path):
    cohort = tmp_path / "c.parquet"
    cohort.write_bytes(b"")
    m1, m2 = _manifest(tmp_path, "a", cohort), _manifest(tmp_path, "b", cohort)
    status, _, rec = _go([m1, m2], code=3)
    assert status == 3 and len(rec.calls) == 1


def test_an_unreadable_manifest_is_an_error(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("not json", encoding="utf-8")
    status, out, rec = _go([bad.as_posix()])
    assert status == 2 and rec.calls == [] and "Cannot read cohort_source" in out
    status, _, _ = _go([(tmp_path / "missing.json").as_posix()])
    assert status == 2


# --- the held-out guard ---------------------------------------------------------
def _cgga_named_cohort(tmp_path, kind):
    cohort = tmp_path / "data" / kind / "cgga" / "eligibility.parquet"
    cohort.parent.mkdir(parents=True)
    cohort.write_bytes(b"")
    return cohort


@pytest.mark.parametrize("kind", ["processed", "raw"])
def test_a_present_cgga_cohort_is_refused_without_the_flag(tmp_path, kind):
    m = _manifest(tmp_path, "ext", _cgga_named_cohort(tmp_path, kind))
    status, out, rec = _go([m])
    assert status == 1
    assert rec.calls == [], "the held-out cohort must not be verified"
    assert "REFUSED" in out and "--allow-held-out" in out


@pytest.mark.parametrize("kind", ["processed", "raw"])
def test_the_flag_allows_it(tmp_path, kind):
    m = _manifest(tmp_path, "ext", _cgga_named_cohort(tmp_path, kind))
    status, _, rec = _go([m], allow=True)
    assert status == 0 and len(rec.calls) == 1


def test_an_absent_cgga_cohort_is_skipped_as_in_ci(tmp_path):
    cohort = tmp_path / "data" / "processed" / "cgga" / "never_created.parquet"
    m = _manifest(tmp_path, "ext", cohort)
    status, out, rec = _go([m])
    assert status == 0 and rec.calls == [] and "Skipping" in out


def test_a_refusal_does_not_hide_the_other_manifests(tmp_path):
    tcga = tmp_path / "t.parquet"
    tcga.write_bytes(b"")
    m_ext = _manifest(tmp_path, "ext", _cgga_named_cohort(tmp_path, "processed"))
    m_tcga = _manifest(tmp_path, "v2", tcga)
    status, _, rec = _go([m_ext, m_tcga])
    assert status == 1 and len(rec.calls) == 1 and rec.calls[0][3] == str(tcga)


def test_main_passes_the_environment_and_the_flag(tmp_path, monkeypatch):
    seen = []

    def fake(manifests, allow_held_out=False):
        seen.append((manifests, allow_held_out))
        return 0

    monkeypatch.setattr(verify, "verify_active", fake)
    monkeypatch.setenv("ACTIVE_MANIFESTS", "a/b.json  c/d.json")
    assert verify.main([]) == 0
    assert verify.main(["--allow-held-out"]) == 0
    assert seen == [(["a/b.json", "c/d.json"], False), (["a/b.json", "c/d.json"], True)]


# --- the workflow -----------------------------------------------------------------
def test_the_ci_step_calls_the_script_and_no_longer_a_copy_of_the_loop():
    doc = yaml.safe_load(CI_YAML.read_text(encoding="utf-8"))
    runs = [
        str(step.get("run", ""))
        for job in doc["jobs"].values()
        for step in job.get("steps", [])
        if isinstance(step, dict)
    ]
    assert any("python scripts/ci_partition_verify.py" in r for r in runs)
    assert not any("for manifest in" in r for r in runs), "loop left inline in ci.yml"
    assert not any("--allow-held-out" in r for r in runs), "CI must never pass the flag"
