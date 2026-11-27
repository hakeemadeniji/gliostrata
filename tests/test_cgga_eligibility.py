"""Tests for scripts/cgga_external_eligibility.py (the T4/T5 eligibility rule).

Everything runs on synthetic frames. No file under data/ is opened.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "cgga_external_eligibility.py"


def _load_module():
    """Load the script by path (scripts/ is not a package)."""
    assert SCRIPT_PATH.is_file(), f"script not found: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location("gliostrata_cgga_eligibility_script", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


elig = _load_module()

#: (patient_id, cgga_dataset, overall_survival_months, vital_status)
ROWS = [
    ("P1", "mRNAseq_693", 14.0, 1),     # eligible
    ("P2", "mRNAseq_325", 3.5, 0),      # eligible
    ("P3", "mRNAseq_693", None, 1),     # RNA-seq, no survival time
    ("P4", "mRNAseq_325", 9.0, None),   # RNA-seq, no vital status
    ("P5", "mRNAseq_693", 20.0, 0),     # RNA-seq row with NO expression column
    ("P6", "mRNA_array_301", 11.0, 1),  # labelled microarray sample: not RNA-seq
    ("P7", "mRNA_array_301", 12.0, 0),  # microarray, no expression
    (None, "mRNAseq_693", 6.0, 1),      # no identifier
]
#: Columns of the synthetic counts matrix. Every RNA-seq sample EXCEPT P5 has one.
EXPR_IDS = {"gene_name", "P1", "P2", "P3", "P4"}


@pytest.fixture
def cohort() -> pd.DataFrame:
    frame = pd.DataFrame(
        ROWS,
        columns=["patient_id", "cgga_dataset", "overall_survival_months", "vital_status"],
    )
    # A covariate that must never reach the sealed table.
    frame["age"] = [55, 60, 61, 62, 63, 64, 65, 66]
    return frame


def test_only_patients_with_rna_seq_expression_and_survival_are_eligible(cohort):
    out, report = elig.eligible_patients(cohort, EXPR_IDS)
    assert list(out["patient_id"]) == ["P1", "P2"]
    assert report["eligible_patients"] == 2
    assert report["eligible_rows"] == 2


def test_microarray_rows_are_excluded_even_when_everything_else_is_present(cohort):
    # P6 is a labelled microarray sample with complete survival. It is still
    # not an RNA-seq patient, so it is not eligible for the T4/T5 partition.
    out, _ = elig.eligible_patients(cohort, EXPR_IDS | {"P6", "P7"})
    assert "P6" not in set(out["patient_id"]) and "P7" not in set(out["patient_id"])


def test_one_row_per_patient_and_duplicates_are_counted(cohort):
    doubled = pd.concat([cohort, cohort.iloc[[0, 1]]], ignore_index=True)
    out, report = elig.eligible_patients(doubled, EXPR_IDS)
    assert list(out["patient_id"]) == ["P1", "P2"]
    assert report["eligible_rows"] == 4
    assert report["duplicate_rows_collapsed"] == 2


def test_output_carries_identifiers_only(cohort):
    out, _ = elig.eligible_patients(cohort, EXPR_IDS)
    assert list(out.columns) == ["patient_id", "cgga_dataset"]
    for forbidden in ("overall_survival_months", "vital_status", "age"):
        assert forbidden not in out.columns
    assert out["patient_id"].is_monotonic_increasing


def test_report_agrees_with_the_ingest_overlap_counts(cohort):
    """The 'with expression and survival' count IS the ingest's own number."""
    from types import SimpleNamespace

    ingest = elig.load_ingest_module()
    mine = elig.eligible_patients(cohort, EXPR_IDS)[1]
    theirs = ingest.expression_overlap(cohort, SimpleNamespace(columns=sorted(EXPR_IDS)))
    assert mine["with_expression"] == theirs["with_expression"]
    assert mine["with_survival"] == theirs["with_survival"]
    assert mine["with_expression_and_survival"] == theirs["with_expression_and_survival"]
    assert mine["with_expression"] == 4      # P1..P4
    assert mine["with_survival"] == 5        # P1, P2, P5, P6, P7; a row without an
    #                                          identifier cannot be counted as a patient
    assert mine["with_expression_and_survival"] == 2
    assert mine["rna_seq_rows_without_survival"] == 2   # P3, P4
    assert mine["rna_seq_rows_without_expression_column"] == 2  # P5 and the row with no id
    assert mine["by_batch"] == {"mRNAseq_325": 1, "mRNAseq_693": 1}


def test_missing_required_column_raises(cohort):
    with pytest.raises(ValueError, match="no 'vital_status' column"):
        elig.eligible_patients(cohort.drop(columns=["vital_status"]), EXPR_IDS)


def test_batch_column_is_optional(cohort):
    out, report = elig.eligible_patients(cohort.drop(columns=["cgga_dataset"]), EXPR_IDS)
    assert list(out.columns) == ["patient_id"]
    assert list(out["patient_id"]) == ["P1", "P2"]
    assert report["by_batch"] == {}


def test_dry_run_writes_nothing_and_the_real_run_writes_identifiers(tmp_path, cohort, capsys):
    cohort_path = tmp_path / "cohort.parquet"
    counts_path = tmp_path / "counts.parquet"
    out_path = tmp_path / "nested" / "eligible.parquet"
    cohort.to_parquet(cohort_path, index=False)
    columns = {"gene_name": ["A", "B"]}
    for pid in sorted(p for p in EXPR_IDS if p != "gene_name"):
        columns[pid] = [1, 2]
    pd.DataFrame(columns).to_parquet(counts_path, index=False)
    common = ["--cohort", str(cohort_path), "--counts", str(counts_path),
              "--out", str(out_path), "--ledger", str(tmp_path / "ledger.jsonl")]

    # The schema-only reader saw the column names; it needs no expression value.
    assert elig.expression_sample_ids(counts_path) == EXPR_IDS

    elig.main([*common, "--dry-run"])
    assert not out_path.exists()
    assert not out_path.parent.exists()
    dry = capsys.readouterr().out
    assert "DRY RUN" in dry
    assert "Eligible rows / eligible patients  : 2 / 2" in dry
    # No clinical value may be printed: none of the survival times appear.
    for value in ROWS:
        if value[2] is not None and value[2] not in (1.0, 2.0):
            assert f"{value[2]:.1f}" not in dry

    elig.main(common)
    real = capsys.readouterr().out
    assert out_path.is_file()
    assert "sha256" in real
    written = pd.read_parquet(out_path)
    assert list(written.columns) == ["patient_id", "cgga_dataset"]
    assert list(written["patient_id"]) == ["P1", "P2"]


# --- read discipline and sealing ------------------------------------------------------------
def _write_inputs(tmp_path, cohort):
    cohort_path = tmp_path / "cohort.parquet"
    counts_path = tmp_path / "counts.parquet"
    cohort.to_parquet(cohort_path, index=False)
    columns = {"gene_name": ["A", "B"]}
    for pid in sorted(p for p in EXPR_IDS if p != "gene_name"):
        columns[pid] = [1, 2]
    pd.DataFrame(columns).to_parquet(counts_path, index=False)
    return cohort_path, counts_path


def test_every_run_is_recorded_in_the_ledger_before_reading(tmp_path, cohort):
    import json

    from gliostrata.evaluation.read_guard import hash_inputs

    cohort_path, counts_path = _write_inputs(tmp_path, cohort)
    ledger = tmp_path / "ledger.jsonl"
    elig.main(["--cohort", str(cohort_path), "--counts", str(counts_path),
               "--out", str(tmp_path / "e.parquet"), "--ledger", str(ledger), "--dry-run"])
    entries = [json.loads(line) for line in ledger.read_text(encoding="utf-8").splitlines()]
    assert len(entries) == 1
    e = entries[0]
    assert (e["cohort"], e["kind"], e["analysis"]) == (
        "cgga", "availability", "external-partition-eligibility")
    assert e["inputs"] == hash_inputs([cohort_path, counts_path])
    assert "Dry run" in e["purpose"] and "never an expression value" in e["purpose"]
    assert "T" in e["recorded_at"]  # a time of day, not a date only


def test_only_the_columns_the_rule_needs_are_read(tmp_path, cohort):
    cohort_path, _ = _write_inputs(tmp_path, cohort)
    frame = elig.read_cohort_columns(cohort_path)
    assert set(frame.columns) == set(elig.COHORT_COLUMNS)
    assert "age" not in frame.columns


def test_seal_freezes_an_all_test_partition_and_verifies_it(tmp_path, cohort, capsys):
    import json

    from gliostrata.evaluation.partitions import read_manifest

    cohort_path, counts_path = _write_inputs(tmp_path, cohort)
    ledger = tmp_path / "ledger.jsonl"
    manifest = tmp_path / "partitions" / "ext.json"
    out = tmp_path / "e.parquet"
    args = ["--cohort", str(cohort_path), "--counts", str(counts_path), "--out", str(out),
            "--ledger", str(ledger), "--seal", str(manifest), "--name", "ext_v1"]
    elig.main(args)
    m = read_manifest(manifest)
    assert m.sealed and m.partition_name == "ext_v1"
    assert sorted(m.test_ids) == ["P1", "P2"] and list(m.train_ids) == []
    assert m.seed == elig.DEFAULT_SEED and m.test_fraction == 1.0
    assert "ledger entry 0" in m.notes
    printed = capsys.readouterr().out
    assert "Sealed ext_v1: 2 test, 0 train" in printed and "verified" in printed
    entry = json.loads(ledger.read_text(encoding="utf-8").splitlines()[0])
    assert "all-test partition" in entry["purpose"]

    # A manifest is never overwritten, and the refusal comes before any access.
    with pytest.raises(SystemExit, match="never overwritten"):
        elig.main(args)
    assert len(ledger.read_text(encoding="utf-8").splitlines()) == 1


def test_seal_cannot_be_combined_with_a_dry_run(tmp_path, cohort):
    cohort_path, counts_path = _write_inputs(tmp_path, cohort)
    with pytest.raises(SystemExit):
        elig.main(["--cohort", str(cohort_path), "--counts", str(counts_path),
                   "--ledger", str(tmp_path / "l.jsonl"), "--dry-run",
                   "--seal", str(tmp_path / "m.json")])
    assert not (tmp_path / "l.jsonl").exists()
