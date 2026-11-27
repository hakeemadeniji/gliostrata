"""Executes the partition-drift CI job script (scripts/ci_partition_drift.py).

What is at stake
----------------
That job is the guardrail against a frozen split drifting. Until 2026-11-27 its
logic was written inline in ci.yml, and the only test of it read the YAML as
text. These tests run the script itself, as the job does, against:

  1. a manifest listed in ACTIVE_MANIFESTS that does not exist  -> the job stops;
  2. a superseded manifest: sealed, valid, not listed            -> labelled, not failed;
  3. an external all-test manifest: every patient held out       -> verified;
  4. an unsealed draft, an edited id list and an unreadable file -> failed;
  5. the repository's own two manifests, with the CI setting      -> verified.

Every manifest here is synthetic and lives under tmp_path. The repository's
own manifests are read, never written.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

from gliostrata.evaluation import freeze_partition

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "ci_partition_drift.py"
CI_YAML = ROOT / ".github" / "workflows" / "ci.yml"
CI_ACTIVE = ["data/partitions/tcga_gbm_clinical_v2.json",
             "data/partitions/cgga_external_v1.json",
             "data/partitions/cgga_t1_subtype_v1.json"]


def _run(cwd, active):
    """Run the job script as CI does: by path, with ACTIVE_MANIFESTS set."""
    env = {**os.environ, "ACTIVE_MANIFESTS": " ".join(active)}
    return subprocess.run(
        [sys.executable, str(SCRIPT)], cwd=cwd, env=env, capture_output=True, text=True
    )


def _ids(prefix, n):
    return [f"{prefix}-{i:03d}" for i in range(n)]


def _freeze(root, name, ids, test_fraction, *, seal=True):
    """Seal (by default) a partition over synthetic ids inside a temporary workspace."""
    root = Path(root)
    source = root / f"{name}_cohort.parquet"
    cohort = pd.DataFrame({"patient_id": list(ids)})
    cohort.to_parquet(source)
    freeze_partition(
        cohort,
        cohort_source=source,
        manifest_path=root / "data" / "partitions" / f"{name}.json",
        partition_name=name,
        seed=20251126,
        test_fraction=test_fraction,
        seal=seal,
    )
    return root / "data" / "partitions" / f"{name}.json"


# --- 1. a listed manifest that is missing -------------------------------------
def test_a_manifest_listed_as_active_but_missing_stops_the_job(tmp_path):
    _freeze(tmp_path, "v1", _ids("A", 12), 0.5)
    p = _run(tmp_path, ["data/partitions/v2.json"])
    assert p.returncode != 0
    assert "ACTIVE_MANIFESTS lists data/partitions/v2.json, which does not exist." in (
        p.stdout + p.stderr
    )


def test_a_directory_with_no_manifest_is_not_a_pass(tmp_path):
    (tmp_path / "data" / "partitions").mkdir(parents=True)
    p = _run(tmp_path, [])
    assert p.returncode != 0
    assert "No partition manifests found" in p.stdout + p.stderr


# --- 2. a superseded manifest -------------------------------------------------
def test_a_superseded_manifest_is_labelled_and_the_active_one_is_checked(tmp_path):
    _freeze(tmp_path, "v1", _ids("A", 12), 0.5)
    _freeze(tmp_path, "v2", _ids("A", 12) + _ids("B", 6), 0.5)
    p = _run(tmp_path, ["data/partitions/v2.json"])
    assert p.returncode == 0, p.stdout + p.stderr
    assert "v1.json  [superseded (historical record)]" in p.stdout
    assert "v2.json  [ACTIVE]" in p.stdout
    assert "cohort       : not checked -- superseded manifests" in p.stdout
    assert "All partition manifests verified." in p.stdout


def test_both_manifests_listed_are_both_treated_as_active(tmp_path):
    _freeze(tmp_path, "v1", _ids("A", 12), 0.5)
    _freeze(tmp_path, "v2", _ids("A", 12) + _ids("B", 6), 0.5)
    p = _run(tmp_path, ["data/partitions/v1.json", "data/partitions/v2.json"])
    assert p.returncode == 0, p.stdout + p.stderr
    assert p.stdout.count("[ACTIVE]") == 2
    assert "superseded" not in p.stdout


# --- 3. an external all-test manifest -----------------------------------------
def test_an_external_all_test_manifest_verifies(tmp_path):
    ids = _ids("EXT", 30)
    _freeze(tmp_path, "ext_v1", ids, 1.0)
    p = _run(tmp_path, ["data/partitions/ext_v1.json"])
    assert p.returncode == 0, p.stdout + p.stderr
    assert "[ACTIVE]" in p.stdout
    assert "n_train/test : 0 / 30" in p.stdout
    assert f"re-derived   : all {len(ids)} ids reproduce" in p.stdout
    assert "full --verify below" in p.stdout


# --- 4. what must fail --------------------------------------------------------
def test_an_unsealed_draft_fails_the_job(tmp_path):
    _freeze(tmp_path, "draft", _ids("A", 12), 0.5, seal=False)
    p = _run(tmp_path, [])
    assert p.returncode == 1
    assert "UNSEALED (draft)" in p.stdout and "PARTITION DRIFT:" in p.stdout


def test_a_sealed_manifest_with_an_edited_id_list_fails(tmp_path):
    path = _freeze(tmp_path, "v2", _ids("B", 12), 0.5)
    assert _run(tmp_path, ["data/partitions/v2.json"]).returncode == 0
    manifest = json.loads(path.read_text(encoding="utf-8"))
    ids = sorted(manifest["train_ids"])
    # Same count, one training patient swapped for one that hashes to the test
    # split: the list no longer matches its digest, and that id no longer
    # re-derives the split it is listed under.
    manifest["train_ids"] = ids[1:] + [manifest["test_ids"][0]]
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    p = _run(tmp_path, ["data/partitions/v2.json"])
    assert p.returncode == 1
    assert "train_ids digest mismatch" in p.stdout
    assert "no longer hash to their frozen split" in p.stdout


def test_an_unreadable_manifest_fails_instead_of_passing(tmp_path):
    (tmp_path / "data" / "partitions").mkdir(parents=True)
    (tmp_path / "data" / "partitions" / "broken.json").write_text("not json", encoding="utf-8")
    p = _run(tmp_path, [])
    assert p.returncode == 1
    assert "could not verify the manifests" in p.stderr


# --- 5. the repository itself, and the workflow that runs this file ------------
def test_the_repository_manifests_pass_with_the_ci_setting():
    p = _run(ROOT, CI_ACTIVE)
    assert p.returncode == 0, p.stdout + p.stderr
    assert "tcga_gbm_clinical_v1.json  [superseded (historical record)]" in p.stdout
    assert "tcga_gbm_clinical_v2.json  [ACTIVE]" in p.stdout
    assert "re-derived   : all 594 ids reproduce" in p.stdout
    assert "cgga_external_v1.json  [ACTIVE]" in p.stdout
    assert "re-derived   : all 218 ids reproduce" in p.stdout
    assert "cgga_t1_subtype_v1.json  [ACTIVE]" in p.stdout
    assert "re-derived   : all 108 ids reproduce" in p.stdout
    # The cohort table is git-ignored: present on the development machine, absent
    # in CI and in any clean checkout. Each case has its own line; both must pass.
    manifest = json.loads((ROOT / CI_ACTIVE[0]).read_text(encoding="utf-8"))
    if (ROOT / manifest["cohort_source"]).exists():
        assert "cohort       : present at" in p.stdout
    else:
        assert "cohort       : NOT in this checkout" in p.stdout
        assert "full --verify not run in CI" in p.stdout
    assert "All partition manifests verified." in p.stdout


def test_the_ci_job_runs_this_file_and_no_longer_a_copy_of_it():
    """The workflow and these tests must execute the same code, not two copies."""
    text = CI_YAML.read_text(encoding="utf-8")
    doc = yaml.safe_load(text)
    runs = [
        str(step.get("run", ""))
        for job in doc["jobs"].values()
        for step in job.get("steps", [])
        if isinstance(step, dict)
    ]
    assert any("python scripts/ci_partition_drift.py" in run for run in runs)
    assert not any("sha256_over_ids" in run for run in runs), "logic left inline in ci.yml"
    for job in doc["jobs"].values():
        for step in job.get("steps", []):
            if isinstance(step, dict) and "ci_partition_drift" in str(step.get("run", "")):
                assert not step.get("continue-on-error", False), "the guardrail must block"
