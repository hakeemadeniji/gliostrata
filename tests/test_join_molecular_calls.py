"""Unit tests for scripts/join_molecular_calls.py (no network access).

What is at stake here
---------------------
The join is the step where a missing molecular call could silently become a
wildtype call. That would shrink the IDH-mutant fraction, invent patients
who were never tested, and make the pre-specified IDH sensitivity analysis
exclude the wrong patients. So the tests are organised around:

  1. Is "no call" kept distinct from "wildtype", both for a patient who is
     absent from a call table and for the far more common case on real data,
     a patient who IS in the table with a null call?
  2. Is nothing imputed (flags False, call values NA, mgmt_score NA)?
  3. Are the output paths derived from --out, and is a re-run idempotent?
  4. Do the manifest counts match the coverage recorded in the project
     docs, when the real call tables are present?
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "join_molecular_calls.py"

REAL_COHORT = ROOT / "data/raw/tcga_gbm/clinical_cohort.parquet"
REAL_IDH = ROOT / "data/raw/tcga_gbm/mutations/idh_calls.parquet"
REAL_MGMT = ROOT / "data/raw/tcga_gbm/methylation/mgmt_calls.parquet"
REAL_EXPRESSION = ROOT / "data/raw/tcga_gbm/expression/tpm.parquet"

#: Coverage of the call tables, as recorded in docs/mgmt_idh_calling.md
#: section 6.3. These are the ground truth the
#: join must reproduce. If one trips, the call tables or the join changed:
#: confirm the change was intended, regenerate the call tables, and update
#: the number here -- never edit the number to match a broken join.
EXPECTED_MANIFEST = {
    "n_clinical_cohort": 594,
    "n_with_idh_call": 368,
    "n_with_mgmt_call": 415,
    "n_with_both": 314,
    "n_with_neither": 125,
    "n_expression": 284,
    "n_expression_with_idh": 248,
    "n_expression_with_mgmt": 230,
    "n_expression_with_both": 206,
    "n_expression_with_neither": 12,
    "n_idh_mutant": 23,
}


def _load_module():
    """Load the script by path (scripts/ is not a package)."""
    assert SCRIPT_PATH.is_file(), f"script not found: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location(
        "gliostrata_join_molecular_calls_script", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def J():
    return _load_module()


# ---------------------------------------------------------------------------
# Fixtures: the real layout, in miniature
# ---------------------------------------------------------------------------
def _cohort(ids):
    return pd.DataFrame({"patient_id": list(ids), "age": [60.0] * len(ids)})


def _idh_calls(rows):
    """rows: (patient_id, idh_mutated) with None for no mutation data."""
    return pd.DataFrame(
        {
            "patient_id": [r[0] for r in rows],
            "idh_mutated": pd.array([r[1] for r in rows], dtype="Int64"),
            "call_basis": ["x"] * len(rows),
        }
    )


def _mgmt_calls(rows):
    """rows: (patient_id, mgmt_methylated, mgmt_score), None where no data."""
    return pd.DataFrame(
        {
            "patient_id": [r[0] for r in rows],
            "mgmt_methylated": pd.array([r[1] for r in rows], dtype="Int64"),
            "mgmt_score": pd.array([r[2] for r in rows], dtype="Float64"),
            "call_basis": ["x"] * len(rows),
        }
    )


def _mini():
    """Four cohort patients; every one of them has a row in each call table.

    P1  IDH wildtype, MGMT methylated      (both called)
    P2  IDH mutant,    MGMT unmethylated   (both called)
    P3  no IDH data,   MGMT unmethylated   (one called)
    P4  no IDH data,   no MGMT data        (neither called)
    PX  outside the cohort: IDH mutant with expression, MGMT no data
    """
    cohort = _cohort(["P1", "P2", "P3", "P4"])
    idh = _idh_calls(
        [("P1", 0), ("P2", 1), ("P3", None), ("P4", None), ("PX", 1)]
    )
    mgmt = _mgmt_calls(
        [
            ("P1", 1, 0.61),
            ("P2", 0, 0.04),
            ("P3", 0, 0.02),
            ("P4", None, None),
            ("PX", None, None),
        ]
    )
    return cohort, idh, mgmt


# ---------------------------------------------------------------------------
# 1. Missing versus wildtype
# ---------------------------------------------------------------------------
def test_patient_with_no_idh_data_is_na_not_wildtype(J):
    cohort, idh, mgmt = _mini()
    joined = J.build_joined(cohort, idh, mgmt)
    row = joined.set_index("patient_id").loc["P3"]
    assert pd.isna(row["idh_mutated"]), "no mutation data must not become 0"
    assert bool(row["has_idh_call"]) is False


def test_patient_absent_from_the_call_table_is_na(J):
    cohort, idh, mgmt = _mini()
    cohort = _cohort(["P1", "P2", "P3", "P4", "P9"])  # P9 is in neither table
    joined = J.build_joined(cohort, idh, mgmt)
    row = joined.set_index("patient_id").loc["P9"]
    assert pd.isna(row["idh_mutated"])
    assert bool(row["has_idh_call"]) is False
    assert pd.isna(row["mgmt_methylated"])
    assert bool(row["has_mgmt_call"]) is False


def test_wildtype_is_zero_with_a_call(J):
    cohort, idh, mgmt = _mini()
    joined = J.build_joined(cohort, idh, mgmt)
    row = joined.set_index("patient_id").loc["P1"]
    assert int(row["idh_mutated"]) == 0
    assert bool(row["has_idh_call"]) is True


def test_a_row_in_the_table_is_not_a_call(J):
    """Regression guard: every cohort patient HAS a row in both call tables.

    Were the flags derived from row presence, all four patients would show
    has_idh_call = True. They must not.
    """
    cohort, idh, mgmt = _mini()
    assert set(cohort["patient_id"]) <= set(idh["patient_id"])
    assert set(cohort["patient_id"]) <= set(mgmt["patient_id"])
    joined = J.build_joined(cohort, idh, mgmt)
    assert joined["has_idh_call"].tolist() == [True, True, False, False]
    assert joined["has_mgmt_call"].tolist() == [True, True, True, False]


def test_no_call_means_no_score(J):
    cohort, idh, mgmt = _mini()
    joined = J.build_joined(cohort, idh, mgmt)
    row = joined.set_index("patient_id").loc["P4"]
    assert pd.isna(row["mgmt_methylated"]), "no methylation data must not become 0"
    assert pd.isna(row["mgmt_score"]), "a score without a call is not a value"
    assert bool(row["has_mgmt_call"]) is False


# ---------------------------------------------------------------------------
# 2. Shape, schema, nothing imputed
# ---------------------------------------------------------------------------
def test_one_row_per_cohort_patient_in_cohort_order(J):
    cohort, idh, mgmt = _mini()
    joined = J.build_joined(cohort, idh, mgmt)
    assert joined["patient_id"].tolist() == cohort["patient_id"].tolist()
    assert "PX" not in set(joined["patient_id"]), "non-cohort patient leaked in"


def test_output_columns_and_dtypes(J):
    cohort, idh, mgmt = _mini()
    joined = J.build_joined(cohort, idh, mgmt)
    assert tuple(joined.columns) == J.JOINED_COLUMNS
    assert str(joined["idh_mutated"].dtype) == "Int8"
    assert str(joined["mgmt_methylated"].dtype) == "Int8"
    assert str(joined["mgmt_score"].dtype) == "float64"
    assert str(joined["has_idh_call"].dtype) == "bool"
    assert str(joined["has_mgmt_call"].dtype) == "bool"


def test_values_are_copied_not_recomputed(J):
    cohort, idh, mgmt = _mini()
    joined = J.build_joined(cohort, idh, mgmt)
    by_id = joined.set_index("patient_id")
    assert by_id.loc["P1", "mgmt_score"] == pytest.approx(0.61)
    assert by_id.loc["P2", "mgmt_score"] == pytest.approx(0.04)
    assert by_id.loc["P1", "idh_mutated"] == 0
    assert by_id.loc["P2", "idh_mutated"] == 1


def test_duplicate_patient_in_a_call_table_is_refused(J):
    cohort, idh, mgmt = _mini()
    bad = _idh_calls([("P1", 0), ("P1", 1)])
    with pytest.raises(ValueError, match="duplicate patient_id"):
        J.build_joined(cohort, bad, mgmt)


def test_duplicate_patient_in_the_cohort_is_refused(J):
    _, idh, mgmt = _mini()
    with pytest.raises(ValueError, match="duplicate patient_id"):
        J.build_joined(_cohort(["P1", "P1"]), idh, mgmt)


def test_call_value_outside_zero_one_is_refused(J):
    cohort, _, mgmt = _mini()
    bad = _idh_calls([("P1", 2)])
    with pytest.raises(ValueError, match="must be 0, 1 or"):
        J.build_joined(cohort, bad, mgmt)


def test_missing_required_column_is_refused(J):
    cohort, idh, mgmt = _mini()
    with pytest.raises(ValueError, match="required column"):
        J.build_joined(cohort, idh.drop(columns=["idh_mutated"]), mgmt)


def test_score_without_a_call_is_refused(J):
    cohort, _, _ = _mini()
    idh = _idh_calls([("P1", 0), ("P2", 1), ("P3", None), ("P4", None)])
    bad = _mgmt_calls(
        [("P1", 1, 0.6), ("P2", 0, 0.04), ("P3", 0, 0.02), ("P4", None, 0.9)]
    )
    with pytest.raises(ValueError, match="mgmt_score is present without"):
        J.build_joined(cohort, idh, bad)


# ---------------------------------------------------------------------------
# 3. Manifest arithmetic, paths, idempotency
# ---------------------------------------------------------------------------
def test_manifest_counts_on_the_mini_fixture(J):
    cohort, idh, mgmt = _mini()
    joined = J.build_joined(cohort, idh, mgmt)
    m = J.build_manifest(joined, ["P1", "P2", "P3", "P4", "PX"], idh, mgmt)
    assert tuple(m) == J.MANIFEST_KEYS
    assert m["n_clinical_cohort"] == 4
    assert m["n_with_idh_call"] == 2
    assert m["n_with_mgmt_call"] == 3
    assert m["n_with_both"] == 2
    assert m["n_with_neither"] == 1
    # PX is outside the cohort: counted in the expression rows only.
    assert m["n_expression"] == 5
    assert m["n_expression_with_idh"] == 3
    assert m["n_expression_with_mgmt"] == 3
    assert m["n_expression_with_both"] == 2
    assert m["n_expression_with_neither"] == 1
    assert m["n_idh_mutant"] == 1, "only cohort patients count here"
    assert m["n_idh_mutant_with_expression"] == 2, "P2 and PX"
    assert m["n_joined_rows"] == 4


def test_manifest_is_json_serialisable_with_int_values(J):
    cohort, idh, mgmt = _mini()
    joined = J.build_joined(cohort, idh, mgmt)
    m = J.build_manifest(joined, ["P1"], idh, mgmt)
    assert all(type(v) is int for v in m.values())
    assert json.loads(json.dumps(m)) == m


def test_output_paths_follow_the_out_argument(J):
    joined, manifest = J.output_paths("elsewhere/nested")
    assert joined == Path("elsewhere/nested") / "molecular_calls_joined.parquet"
    assert manifest == Path("elsewhere/nested") / "molecular_calls_manifest.json"


def _write_mini_inputs(tmp_path, samples):
    cohort, idh, mgmt = _mini()
    paths = {
        "cohort": tmp_path / "cohort.parquet",
        "idh": tmp_path / "idh.parquet",
        "mgmt": tmp_path / "mgmt.parquet",
        "expression": tmp_path / "tpm.parquet",
    }
    cohort.to_parquet(paths["cohort"], index=False)
    idh.to_parquet(paths["idh"], index=False)
    mgmt.to_parquet(paths["mgmt"], index=False)
    # Sample columns carry the TCGA barcode prefix the real matrix has; the
    # script reads patient ids from the columns with that prefix.
    pd.DataFrame({f"TCGA-{s}": [1.0] for s in samples}).to_parquet(
        paths["expression"]
    )
    return paths


def _argv(paths, out):
    return [
        "--cohort", str(paths["cohort"]),
        "--idh-calls", str(paths["idh"]),
        "--mgmt-calls", str(paths["mgmt"]),
        "--expression", str(paths["expression"]),
        "--out", str(out),
    ]


def test_run_writes_both_outputs_under_out(J, tmp_path):
    paths = _write_mini_inputs(tmp_path, ["P1"])
    out = tmp_path / "derived"
    m = J.run(
        paths["cohort"], paths["idh"], paths["mgmt"], paths["expression"], out
    )
    joined_path, manifest_path = J.output_paths(out)
    assert joined_path.is_file() and manifest_path.is_file()
    assert json.loads(manifest_path.read_text(encoding="ascii")) == m
    assert m["n_expression"] == 1

    # Round trip: the NA-versus-wildtype distinction survives the parquet.
    back = pd.read_parquet(joined_path).set_index("patient_id")
    assert pd.isna(back.loc["P3", "idh_mutated"])
    assert back.loc["P1", "idh_mutated"] == 0
    assert not back.loc["P3", "has_idh_call"]


def test_cli_is_idempotent(J, tmp_path):
    paths = _write_mini_inputs(tmp_path, ["P1", "PX"])
    argv = _argv(paths, tmp_path / "out")
    J.main(argv)
    joined_path, manifest_path = J.output_paths(tmp_path / "out")
    first_manifest = manifest_path.read_bytes()
    first_joined = pd.read_parquet(joined_path)

    J.main(argv)

    assert manifest_path.read_bytes() == first_manifest, "manifest not idempotent"
    pd.testing.assert_frame_equal(pd.read_parquet(joined_path), first_joined)


# ---------------------------------------------------------------------------
# 4. The real call tables (skipped cleanly when data/raw is absent)
# ---------------------------------------------------------------------------
real_data = pytest.mark.skipif(
    not all(
        p.exists()
        for p in (REAL_COHORT, REAL_IDH, REAL_MGMT, REAL_EXPRESSION)
    ),
    reason="TCGA-GBM call tables not present",
)


@real_data
def test_real_call_tables_reproduce_the_documented_coverage(J):
    cohort = pd.read_parquet(REAL_COHORT, columns=["patient_id"])
    idh = pd.read_parquet(REAL_IDH, columns=["patient_id", "idh_mutated"])
    mgmt = pd.read_parquet(
        REAL_MGMT, columns=["patient_id", "mgmt_methylated", "mgmt_score"]
    )
    joined = J.build_joined(cohort, idh, mgmt)
    expr = J.expression_patient_ids(REAL_EXPRESSION)

    assert len(joined) == EXPECTED_MANIFEST["n_clinical_cohort"]
    assert joined["patient_id"].tolist() == cohort["patient_id"].astype(str).tolist()

    m = J.build_manifest(joined, expr, idh, mgmt)
    for key, expected in EXPECTED_MANIFEST.items():
        assert m[key] == expected, f"{key}: expected {expected}, got {m[key]}"

    assert m["n_idh_mutant"] <= m["n_with_idh_call"]
    assert 0 < m["n_idh_mutant_with_expression"] <= m["n_idh_mutant"]
    assert m["n_joined_rows"] == 594


@real_data
def test_real_joined_table_keeps_the_missing_distinct_from_wildtype(J):
    """On the real tables every cohort patient has a row in both call files.

    If the flags were ever derived from row presence, has_idh_call would be
    True for all 594. It must be True for exactly the called patients.
    """
    cohort = pd.read_parquet(REAL_COHORT, columns=["patient_id"])
    idh = pd.read_parquet(REAL_IDH, columns=["patient_id", "idh_mutated"])
    mgmt = pd.read_parquet(
        REAL_MGMT, columns=["patient_id", "mgmt_methylated", "mgmt_score"]
    )
    assert set(cohort["patient_id"].astype(str)) <= set(idh["patient_id"])
    assert set(cohort["patient_id"].astype(str)) <= set(mgmt["patient_id"])

    joined = J.build_joined(cohort, idh, mgmt)

    no_idh = joined[~joined["has_idh_call"]]
    assert len(no_idh) == 594 - EXPECTED_MANIFEST["n_with_idh_call"]
    assert no_idh["idh_mutated"].isna().all(), "a missing call became a call"

    no_mgmt = joined[~joined["has_mgmt_call"]]
    assert no_mgmt["mgmt_methylated"].isna().all()
    assert no_mgmt["mgmt_score"].isna().all()

    called = joined[joined["has_idh_call"]]
    assert called["idh_mutated"].isin([0, 1]).all()
    assert int((called["idh_mutated"] == 0).sum()) == (
        EXPECTED_MANIFEST["n_with_idh_call"] - EXPECTED_MANIFEST["n_idh_mutant"]
    )
