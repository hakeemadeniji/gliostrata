"""Tests for the held-out cohort read guard (gliostrata.evaluation.read_guard).

No test here reads CGGA. The repository ledger (data/access_ledger.jsonl) is
read as text; every access it describes is re-created on temporary files.

What is at stake
----------------
Four things would each produce a plausible-looking run that the
pre-registration does not allow, and each has a test:

  1. A value read before the registration is deposited.
  2. A value read of a metric that was never registered.
  3. A value read with neither a sealed partition nor a registered scoring rule.
  4. A second value read of an analysis that has already been read.

A further test keeps the ledger and the pre-registration from disagreeing
about how many times the external cohort has been read.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from gliostrata.cli import cli
from gliostrata.evaluation import freeze_partition, registered_metrics
from gliostrata.evaluation.read_guard import (
    KINDS,
    LEDGER_PATH,
    REGISTRATION_PATH,
    VALUES,
    Access,
    ReadRefusedError,
    check_access,
    format_ledger,
    hash_inputs,
    is_deposited,
    load_ledger,
    record_access,
    registration_state,
    summarize,
)

ROOT = Path(__file__).resolve().parents[1]
METRIC = "T4_optimizer_regret_reduction_pct"
MANIFEST = "data/partitions/ext_v1.json"

# Synthetic: it exists only in the temporary files a test writes.
FIXTURE_DOI = "10.0000/osf.io/test-fixture"
NOT_DEPOSITED = dict.fromkeys(("osf_registration_doi", "deposited_at", "registered_commit"))


def _values(**overrides):
    """A value read that satisfies every gate once the fixtures are in place."""
    fields = dict(
        cohort="cgga",
        kind=VALUES,
        analysis=METRIC,
        purpose="fixture",
        inputs={"data/raw/cgga/clinical_cohort.parquet": "0" * 64},
        partition=MANIFEST,
    )
    fields.update(overrides)
    return Access(**fields)


def _deposited():
    return {**NOT_DEPOSITED, "osf_registration_doi": FIXTURE_DOI}


def _seal(tmp_path, *, seal=True):
    """Write a (by default sealed) all-test manifest over 40 synthetic ids."""
    cohort = pd.DataFrame({"patient_id": [f"P{i:03d}" for i in range(40)]})
    source = tmp_path / "cohort.parquet"
    cohort.to_parquet(source)
    path = tmp_path / MANIFEST
    freeze_partition(
        cohort,
        cohort_source=source,
        manifest_path=path,
        partition_name="ext_v1",
        seed=20251126,
        test_fraction=1.0,
        seal=seal,
    )
    return path


# ---------------------------------------------------------------------------
# The gates
# ---------------------------------------------------------------------------
def test_only_a_value_access_is_gated(tmp_path):
    ledger = tmp_path / "ledger.jsonl"
    for kind in KINDS:
        if kind != VALUES:
            record_access(
                Access("cgga", kind, "look", "counts only"),
                ledger_path=ledger,
                registration_path=tmp_path / "missing.yaml",
                root=tmp_path,
            )
    summary = summarize(load_ledger(ledger))
    assert summary["entries"] == len(KINDS) - 1
    assert summary["value_reads"] == 0
    assert summary["before_deposit"] == len(KINDS) - 1


def test_a_value_read_is_refused_before_the_deposit(tmp_path):
    _seal(tmp_path)
    reasons = check_access(_values(), ledger=[], registration=NOT_DEPOSITED, root=tmp_path)
    assert len(reasons) == 1
    assert "REGISTRATION GATE" in reasons[0]


def test_placeholder_text_does_not_count_as_a_deposit():
    assert not is_deposited({**NOT_DEPOSITED, "osf_registration_doi": "TODO"})
    assert not is_deposited({**NOT_DEPOSITED, "osf_registration_doi": "pending"})
    assert not is_deposited({**NOT_DEPOSITED, "osf_registration_doi": 20261126})
    assert is_deposited(_deposited())


def test_an_unregistered_metric_is_refused(tmp_path):
    _seal(tmp_path)
    reasons = check_access(
        _values(analysis="T9_something_new"), ledger=[], registration=_deposited(),
        root=tmp_path,
    )
    assert len(reasons) == 1
    assert "not one of the registered metrics" in reasons[0]
    assert "T9_something_new" in reasons[0]


@pytest.mark.parametrize("metric", registered_metrics())
def test_every_registered_metric_passes_the_registration_gate(tmp_path, metric):
    _seal(tmp_path)
    assert (
        check_access(
            _values(analysis=metric), ledger=[], registration=_deposited(), root=tmp_path
        )
        == []
    )


def test_an_exploratory_read_needs_its_own_preregistration(tmp_path):
    _seal(tmp_path)
    access = _values(analysis="exploratory:calibration-slope")
    reasons = check_access(access, ledger=[], registration=_deposited(), root=tmp_path)
    assert len(reasons) == 1
    assert "exploratory" in reasons[0]
    cited = _values(
        analysis="exploratory:calibration-slope",
        exploratory_preregistration="osf.io/example (a citation that exists only in this test)",
    )
    assert check_access(cited, ledger=[], registration=_deposited(), root=tmp_path) == []


def test_partition_gate(tmp_path):
    missing = check_access(_values(), ledger=[], registration=_deposited(), root=tmp_path)
    assert [r for r in missing if "PARTITION GATE" in r and "does not exist" in r]

    _seal(tmp_path, seal=False)
    draft = check_access(_values(), ledger=[], registration=_deposited(), root=tmp_path)
    assert [r for r in draft if "PARTITION GATE" in r and "unsealed draft" in r]

    neither = check_access(
        _values(partition=None), ledger=[], registration=_deposited(), root=tmp_path
    )
    assert [r for r in neither if "PARTITION GATE" in r]

    # T1 and T3 are fixed by a registered rule instead of a sealed partition.
    by_rule = check_access(
        _values(partition=None, scoring_rule="docs/preregistration.md section 2.2 item 9"),
        ledger=[],
        registration=_deposited(),
        root=tmp_path,
    )
    assert by_rule == []


def test_second_value_read_needs_a_repeat_reason(tmp_path):
    _seal(tmp_path)
    first = {"cohort": "cgga", "analysis": METRIC, "kind": VALUES, "date": "2026-11-26"}
    reasons = check_access(_values(), ledger=[first], registration=_deposited(), root=tmp_path)
    assert len(reasons) == 1
    assert "ONCE PER ANALYSIS" in reasons[0] and "2026-11-26" in reasons[0]

    again = _values(repeat_reason="disclosed re-run after a documented amendment")
    assert check_access(again, ledger=[first], registration=_deposited(), root=tmp_path) == []

    # Another analysis is not spent by this read, and neither is a dry run.
    dry_run = {"cohort": "cgga", "analysis": METRIC, "kind": "availability"}
    other = _values(analysis="T5_endtoend_rmst_success_fraction")
    assert check_access(other, ledger=[dry_run], registration=_deposited(), root=tmp_path) == []


def test_value_read_without_input_checksums_is_refused(tmp_path):
    _seal(tmp_path)
    reasons = check_access(
        _values(inputs={}), ledger=[], registration=_deposited(), root=tmp_path
    )
    assert len(reasons) == 1
    assert "SHA-256" in reasons[0]


def test_all_reasons_are_reported_at_once(tmp_path):
    reasons = check_access(
        _values(partition=None), ledger=[], registration=NOT_DEPOSITED, root=tmp_path
    )
    joined = " | ".join(reasons)
    assert "REGISTRATION GATE" in joined
    assert "PARTITION GATE" in joined


def test_empty_purpose_unknown_cohort_or_unknown_kind_is_refused():
    assert check_access(Access("cgga", "schema", "x", " "))
    assert check_access(Access("cgga", "schema", " ", "p"))
    assert check_access(Access("tcga", "schema", "x", "development data"))
    assert check_access(Access("cgga", "guess", "x", "p"))
    assert check_access(Access("cgga", "schema", "x", "counts only")) == []


# ---------------------------------------------------------------------------
# Recording
# ---------------------------------------------------------------------------
def test_record_access_appends_and_a_refused_read_is_not_written(tmp_path):
    _seal(tmp_path)
    config = tmp_path / "configs" / "registration.yaml"
    config.parent.mkdir()
    config.write_text(f"osf_registration_doi: '{FIXTURE_DOI}'\n", encoding="utf-8")
    ledger = tmp_path / "ledger.jsonl"
    kwargs = dict(ledger_path=ledger, registration_path=config, root=tmp_path)

    entry = record_access(_values(), **kwargs)
    assert entry["seq"] == 0
    assert entry["backfilled"] is False
    assert entry["before_deposit"] is False

    second = record_access(
        _values(analysis="T5_endtoend_rmst_success_fraction"), **kwargs
    )
    assert second["seq"] == 1

    with pytest.raises(ReadRefusedError, match="ONCE PER ANALYSIS"):
        record_access(_values(), **kwargs)
    assert len(load_ledger(ledger)) == 2


def test_a_value_read_is_refused_while_not_deposited(tmp_path):
    _seal(tmp_path)
    ledger = tmp_path / "ledger.jsonl"
    with pytest.raises(ReadRefusedError, match="REGISTRATION GATE"):
        record_access(
            _values(),
            ledger_path=ledger,
            registration_path=tmp_path / "absent.yaml",
            root=tmp_path,
        )
    assert not ledger.exists()


def test_hash_inputs_hashes_bytes_only(tmp_path):
    path = tmp_path / "table.parquet"
    pd.DataFrame({"patient_id": ["a", "b"], "value": [1.5, 2.5]}).to_parquet(path)
    expected = hashlib.sha256(path.read_bytes()).hexdigest()
    assert hash_inputs([path]) == {path.as_posix(): expected}


def test_a_ledger_line_that_was_edited_is_a_finding(tmp_path):
    ledger = tmp_path / "ledger.jsonl"
    ledger.write_text('{"seq": 0}\nnot json\n', encoding="utf-8")
    with pytest.raises(ValueError, match="append-only"):
        load_ledger(ledger)
    assert load_ledger(tmp_path / "absent.jsonl") == []


def test_format_ledger_numbers_value_reads_in_order():
    kinds = ("schema", VALUES, "availability", VALUES)
    entries = [
        {"seq": i, "date": "2026-11-26", "kind": k, "analysis": "a", "purpose": "p"}
        for i, k in enumerate(kinds)
    ]
    text = format_ledger(entries)
    assert "value read 1" in text
    assert "value read 2" in text
    assert "Value reads of the external cohort: 2 (a, a)." in text
    assert "empty" in format_ledger([])


# ---------------------------------------------------------------------------
# The repository's own ledger and registration state
# ---------------------------------------------------------------------------
def test_the_repository_is_not_registered_yet():
    state = registration_state(ROOT / REGISTRATION_PATH)
    assert state["config_found"] is True
    assert not is_deposited(state)


def test_the_repository_ledger_agrees_with_the_registration():
    entries = load_ledger(ROOT / LEDGER_PATH)
    assert [entry["seq"] for entry in entries] == list(range(len(entries)))
    assert all(entry["cohort"] == "cgga" and entry["kind"] in KINDS for entry in entries)
    assert all(entry["before_deposit"] for entry in entries)

    reads = summarize(entries)
    assert reads["value_reads"] == 2
    assert reads["value_read_analyses"] == ["T2_survival_cindex_cgga", "T2_brier_18mo"]

    cohort_file = "data/raw/cgga/clinical_cohort.parquet"
    digests = {entry["inputs"][cohort_file] for entry in entries if entry["kind"] == VALUES}
    assert len(digests) == 1
    assert len(next(iter(digests))) == 64

    preregistration = (ROOT / "docs" / "preregistration.md").read_text(encoding="utf-8")
    assert "CGGA has been read twice, both for registered rows of T2" in preregistration


def test_registration_config_and_preregistration_header_agree():
    """Deposit day: this fails until the header and the config are both updated."""
    header = (ROOT / "docs" / "preregistration.md").read_text(encoding="utf-8")[:4000]
    doi = registration_state(ROOT / REGISTRATION_PATH)["osf_registration_doi"]
    if doi:
        assert doi in header, "the DOI is in the config but not in the header"
    else:
        assert "TO BE COMPLETED BY AUTHOR: OSF DOI" in header
        assert "not yet deposited" in header


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def test_cli_prints_the_ledger_and_the_registration_state(tmp_path):
    ledger = tmp_path / "ledger.jsonl"
    record_access(
        Access("cgga", "schema", "look", "counts only"),
        ledger_path=ledger,
        registration_path=tmp_path / "absent.yaml",
        root=tmp_path,
    )
    result = CliRunner().invoke(cli, ["access-ledger", "--ledger", str(ledger)])
    assert result.exit_code == 0, result.output
    assert "1 entries" in result.output
    assert "#0" in result.output
    assert "counts only" in result.output
    assert "Registration:" in result.output


def test_cli_prints_an_empty_ledger(tmp_path):
    ledger = tmp_path / "ledger.jsonl"
    ledger.write_text("", encoding="utf-8")
    result = CliRunner().invoke(cli, ["access-ledger", "--ledger", str(ledger)])
    assert result.exit_code == 0, result.output
    assert "never been accessed" in result.output
