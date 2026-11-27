"""Tests for scripts/build_tcga_feature_store.py (synthetic inputs only).

What is at stake here
---------------------
This script produces the first feature store built from real data. Five things
can go wrong, and each would produce a plausible parquet file:

  1. A cohort column shadowing a computed value. The harmonised cohort
     carries placeholder columns under the names of the L2 outputs. Tests
     check that all-null placeholders are removed and that a column with
     values stops the script.
  2. A patient dropped without being counted. A test checks every exclusion
     reason against the summary, and that they add up to the cohort.
  3. A missing value turned into a number. Tests check that the defaults
     exclude, and that ``zero`` is an explicit choice that is still reported.
  4. A field absent for everyone read as "every patient is in the reference
     group". A test checks that it is reported, not hidden.
  5. A read of the held-out cohort. A test checks that a CGGA path is refused.

TCGA is not read either: every table below is synthetic. The expression matrix
is a small gene x sample matrix holding the genes of the three signatures,
indexed by HGNC symbol (annotation=None).
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gliostrata.features.store import FEATURE_COLUMNS
from gliostrata.profiling.deconvolve import NEFTEL_SIGNATURES, TME_SIGNATURES
from gliostrata.profiling.verhaak import SUBTYPES, VERHAAK_SIGNATURES

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "build_tcga_feature_store.py"
STATES = ("npc", "opc", "ac", "mes")
TMES = ("immune", "stromal", "endothelial")


def _load_module():
    assert SCRIPT_PATH.is_file(), f"script not found: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location(
        "gliostrata_build_tcga_feature_store_script", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def script():
    return _load_module()


def _inputs(seed=0, n=11):
    """Cohort, expression, annotation and joined calls for ``n`` patients.

    Patients 0 to n-2 have an expression sample; the last patient has none.
    The expression matrix holds one sample outside the cohort ("OUT").
    The joined table has rows for patients 0 to n-3 only: patient n-3 has no
    MGMT call and patient n-2 has no row at all. Resection extent is absent
    for everyone, as it is in the TCGA-GBM clinical table.
    """
    rng = np.random.default_rng(seed)
    ids = [f"TCGA-SYN-{i:02d}" for i in range(n)]
    genes = set()
    for signatures in (VERHAAK_SIGNATURES, NEFTEL_SIGNATURES, TME_SIGNATURES):
        for panel in signatures.values():
            genes.update(panel)
    genes = sorted(genes)
    n_expr = n - 1
    cohort = pd.DataFrame({
        "patient_id": ids,
        "age": np.linspace(30.0, 75.0, n),
        "sex": (["M", "F"] * n)[:n],
        "kps": np.tile([80.0, 70.0, 90.0, 100.0], n)[:n],
        "extent_of_resection": [None] * n,
        "overall_survival_months": rng.uniform(1.0, 60.0, n),
        "vital_status": rng.integers(0, 2, n).astype(float),
    })
    samples = ids[:n_expr] + ["OUT"]
    values = rng.lognormal(mean=2.0, sigma=1.5, size=(len(genes), len(samples)))
    expression = pd.DataFrame(values, index=genes, columns=samples)
    n_joined = n_expr - 1
    mgmt = ([1, 0] * n_joined)[:n_joined - 1] + [pd.NA]
    idh = [int(i % 3 == 1) for i in range(n_joined)]
    joined = pd.DataFrame({
        "patient_id": ids[:n_joined],
        "mgmt_methylated": pd.array(mgmt, dtype="Int8"),
        "idh_mutated": pd.array(idh, dtype="Int8"),
    })
    return cohort, expression, None, joined


def _ids(store):
    return store["patient_id"].tolist()


def _write_inputs(tmp_path):
    cohort, expression, _, joined = _inputs()
    # An Ensembl-indexed matrix with its annotation, as the real inputs are.
    symbols = expression.index
    expression = expression.copy()
    expression.index = [f"ENSG-SYN-{i:06d}" for i in range(len(expression))]
    frames = {
        "cohort": cohort,
        "expression": expression,
        "joined": joined,
        "annotation": pd.DataFrame({"gene_id": expression.index, "gene_name": symbols}),
    }
    paths = {}
    for name, frame in frames.items():
        paths[name] = tmp_path / f"{name}.parquet"
        frame.to_parquet(paths[name])
    return paths


# --- 2. every exclusion is counted -------------------------------------------
def test_every_exclusion_is_counted_and_they_add_up(script):
    store, r = script.build(*_inputs())
    assert _ids(store) == [f"TCGA-SYN-{i:02d}" for i in range(8)]
    assert r["n_cohort"] == 11
    assert r["n_expression_samples"] == 11
    assert r["n_expression_samples_outside_cohort"] == 1
    assert r["n_scoring_set"] == 10
    assert r["n_excluded_without_expression"] == 1
    assert r["n_excluded_missing_demographic"] == 0
    assert r["n_excluded_missing_molecular"] == 2
    assert r["n_missing"] == {
        "without_joined_row": 1, "mgmt_call": 2, "idh_call": 1,
        "molecular_any": 2, "age": 0, "kps": 0, "sex": 0, "extent_of_resection": 10,
    }
    assert r["n_patients"] == len(store) == 8
    assert (r["n_patients"] + r["n_excluded_without_expression"]
            + r["n_excluded_missing_demographic"]
            + r["n_excluded_missing_molecular"]) == r["n_cohort"]


def test_store_holds_the_feature_columns_and_the_positive_calls(script):
    store, r = script.build(*_inputs())
    for col in ("patient_id", *FEATURE_COLUMNS, "subtype",
                "overall_survival_months", "vital_status"):
        assert col in store.columns
    assert store["mgmt_methylated"].isin([0, 1]).all()
    assert store["idh_mutated"].isin([0, 1]).all()
    assert int(store["mgmt_methylated"].sum()) == r["n_positive"]["mgmt_methylated"] == 4
    assert int(store["idh_mutated"].sum()) == r["n_positive"]["idh_mutated"] == 3
    assert sum(r["subtype_counts"].values()) == len(store)


# --- 3. nothing is imputed silently ------------------------------------------
def test_zero_policy_keeps_the_patient_and_is_reported(script):
    store, r = script.build(*_inputs(), missing_molecular="zero")
    assert r["n_patients"] == 10
    assert r["n_excluded_missing_molecular"] == 0
    row = store.set_index("patient_id").loc["TCGA-SYN-08"]
    assert (row["mgmt_methylated"], row["idh_mutated"]) == (0, 0)
    # The absent calls are still counted, whichever policy is used.
    assert r["n_missing"]["mgmt_call"] == 2
    assert r["policies"]["missing_molecular"] == "zero"


def test_missing_kps_stops_by_default_and_can_be_dropped(script):
    cohort, expression, annotation, joined = _inputs()
    cohort.loc[0, "kps"] = np.nan
    with pytest.raises(ValueError, match="kps"):
        script.build(cohort, expression, annotation, joined)
    store, r = script.build(cohort, expression, annotation, joined, missing_demographic="drop")
    assert "TCGA-SYN-00" not in _ids(store)
    assert r["n_excluded_missing_demographic"] == 1
    assert r["n_missing"]["kps"] == 1
    assert r["n_patients"] == 7


@pytest.mark.parametrize("kwargs", [{"missing_molecular": "impute"},
                                    {"missing_demographic": "ignore"}])
def test_unknown_policy_is_refused(script, kwargs):
    with pytest.raises(ValueError, match="must be one of"):
        script.build(*_inputs(), **kwargs)


# --- values are what the L2 and L3 code say they are --------------------------
def test_values_are_finite_and_normalised(script):
    store, r = script.build(*_inputs())
    assert np.isfinite(store[list(FEATURE_COLUMNS)].to_numpy(dtype=float)).all()
    assert store[list(FEATURE_COLUMNS)].notna().to_numpy().all()
    np.testing.assert_allclose(store[[f"state_{s}" for s in STATES]].sum(axis=1), 1.0)
    np.testing.assert_allclose(store[[f"tme_{s}" for s in TMES]].sum(axis=1), 1.0)
    np.testing.assert_allclose(store[[f"subtype_prob_{s}" for s in SUBTYPES]].sum(axis=1), 1.0)
    assert r["row_sum_range"]["state"] == pytest.approx([1.0, 1.0])
    assert r["row_sum_range"]["tme"] == pytest.approx([1.0, 1.0])
    assert set(store["subtype"]) <= set(SUBTYPES)
    assert r["subtype_counts"] == {
        st: int((store["subtype"] == st).sum()) for st in SUBTYPES
    }
    for col, rng in r["feature_ranges"].items():
        assert rng["min"] <= rng["median"] <= rng["max"]
    assert set(r["feature_ranges"]) == set(FEATURE_COLUMNS)


def test_build_is_deterministic(script):
    first, _ = script.build(*_inputs(seed=3))
    second, _ = script.build(*_inputs(seed=3))
    pd.testing.assert_frame_equal(first, second)


def test_expression_samples_outside_the_cohort_do_not_affect_the_store(script):
    cohort, expression, annotation, joined = _inputs()
    store, _ = script.build(cohort, expression, annotation, joined)
    assert "OUT" not in _ids(store)
    without, _ = script.build(cohort, expression.drop(columns=["OUT"]), annotation, joined)
    pd.testing.assert_frame_equal(without, store)


def test_subtype_probabilities_depend_on_the_scoring_set(script):
    cohort, expression, annotation, joined = _inputs()
    store, _ = script.build(cohort, expression, annotation, joined)
    smaller, r = script.build(cohort.drop(index=[5, 6]), expression, annotation, joined)
    assert r["n_scoring_set"] == 8
    common = [i for i in _ids(smaller) if i in set(_ids(store))]
    assert len(common) == 6
    probs = [f"subtype_prob_{s}" for s in SUBTYPES]
    left = store.set_index("patient_id").loc[common, probs].to_numpy()
    right = smaller.set_index("patient_id").loc[common, probs].to_numpy()
    # A within-cohort ranking: the same patient gets different probabilities
    # when other patients leave the scoring set.
    assert not np.allclose(left, right)


def test_ensembl_indexed_matrix_with_annotation_gives_the_same_store(script):
    cohort, expression, _, joined = _inputs()
    plain, _ = script.build(cohort, expression, None, joined)
    annotated = expression.copy()
    annotated.index = [f"ENSG-SYN-{i:06d}" for i in range(len(annotated))]
    annotation = pd.DataFrame({"gene_id": annotated.index, "gene_name": expression.index})
    store, _ = script.build(cohort, annotated, annotation, joined)
    pd.testing.assert_frame_equal(store, plain)


# --- 1. a cohort column must not shadow a computed value -----------------------
def test_all_null_placeholders_are_removed_and_change_nothing(script):
    cohort, expression, annotation, joined = _inputs()
    plain, _ = script.build(cohort, expression, annotation, joined)
    for col in script.SHADOWED_COLUMNS:
        cohort[col] = None
    assert set(script.SHADOWED_COLUMNS) <= set(cohort.columns)
    store, _ = script.build(cohort, expression, annotation, joined)
    pd.testing.assert_frame_equal(store, plain)
    assert store["subtype"].notna().all()
    derived = [f"state_{s}" for s in STATES] + [f"tme_{s}" for s in TMES]
    assert store[derived].notna().to_numpy().all()


@pytest.mark.parametrize("column, value",
                         [("subtype", "Classical"), ("state_npc", 0.5),
                          ("mgmt_methylated", 1)])
def test_a_cohort_column_with_values_under_a_computed_name_is_refused(script, column, value):
    cohort, expression, annotation, joined = _inputs()
    cohort[column] = value
    with pytest.raises(ValueError, match="shadow"):
        script.build(cohort, expression, annotation, joined)


# --- 4. a field absent for everyone is reported, not hidden --------------------
def test_a_field_absent_for_everyone_is_reported_not_dropped(script):
    store, r = script.build(*_inputs())
    assert r["not_measured_in_source"] == ["extent_of_resection"]
    assert r["constant_features"] == ["extent_of_resection_STR",
                                      "extent_of_resection_GTR"]
    assert (store["extent_of_resection_STR"] == 0).all()
    assert (store["extent_of_resection_GTR"] == 0).all()
    assert r["extent_of_resection_counts"] == {"None": 8}
    assert r["n_excluded_missing_demographic"] == 0


def test_resection_is_encoded_when_it_is_present(script):
    cohort, expression, annotation, joined = _inputs()
    cohort["extent_of_resection"] = (["GTR", "STR", "Biopsy"] * 4)[:len(cohort)]
    store, r = script.build(cohort, expression, annotation, joined)
    rows = store.set_index("patient_id")
    gtr, strr = "extent_of_resection_GTR", "extent_of_resection_STR"
    assert (rows.loc["TCGA-SYN-00", gtr], rows.loc["TCGA-SYN-00", strr]) == (1, 0)
    assert (rows.loc["TCGA-SYN-01", gtr], rows.loc["TCGA-SYN-01", strr]) == (0, 1)
    assert (rows.loc["TCGA-SYN-02", gtr], rows.loc["TCGA-SYN-02", strr]) == (0, 0)
    assert r["not_measured_in_source"] == []
    assert r["extent_of_resection_counts"] == {"GTR": 3, "STR": 3, "Biopsy": 2}
    assert strr not in r["constant_features"]


# --- malformed inputs are refused, never repaired ------------------------------
def test_duplicate_cohort_patient_is_refused(script):
    cohort, expression, annotation, joined = _inputs()
    with pytest.raises(ValueError, match="duplicate"):
        script.build(pd.concat([cohort, cohort.iloc[[0]]]), expression, annotation, joined)


def test_missing_required_column_is_refused(script):
    cohort, expression, annotation, joined = _inputs()
    with pytest.raises(ValueError, match="required column"):
        script.build(cohort.drop(columns=["kps"]), expression, annotation, joined)
    with pytest.raises(ValueError, match="required column"):
        script.build(cohort, expression, annotation, joined.drop(columns=["mgmt_methylated"]))


def test_duplicate_expression_sample_is_refused(script):
    cohort, expression, annotation, joined = _inputs()
    doubled = pd.concat([expression, expression[["OUT"]]], axis=1)
    with pytest.raises(ValueError, match="duplicate sample"):
        script.build(cohort, doubled, annotation, joined)


def test_scorer_returning_other_patients_is_refused(script):
    def wrong(expression, annotation):
        out = script.classify_verhaak(expression, annotation)
        out["patient_id"] = out["patient_id"] + "x"
        return out

    with pytest.raises(ValueError, match="do not match"):
        script.build(*_inputs(), classify=wrong)


def test_a_cohort_without_any_expression_sample_is_refused(script):
    cohort, expression, annotation, joined = _inputs()
    with pytest.raises(ValueError, match="no cohort patient has an expression sample"):
        script.build(cohort, expression[["OUT"]], annotation, joined)


# --- 5. the held-out cohort is never read --------------------------------------
def test_the_held_out_cohort_is_refused(script, tmp_path):
    with pytest.raises(ValueError, match="TCGA only"):
        script.refuse_external("data/raw/CGGA/counts.parquet")
    absent = tmp_path / "data" / "raw" / "cgga" / "cohort.parquet"
    # Refused before anything is opened: the file does not exist.
    with pytest.raises(ValueError, match="TCGA only"):
        script.run(absent, "expression", "annotation", "joined",
              tmp_path / "store.parquet", tmp_path / "report.json")


def test_sha256_hashes_the_whole_file_across_chunks(script, tmp_path):
    small = tmp_path / "abc.txt"
    small.write_bytes(b"abc")
    assert script.sha256_file(small) == (
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    )
    payload = bytes(range(256)) * ((1 << 20) // 256 + 3)
    big = tmp_path / "big.bin"
    big.write_bytes(payload)
    assert len(payload) > (1 << 20)
    assert script.sha256_file(big) == hashlib.sha256(payload).hexdigest()


# --- the files on disk ---------------------------------------------------------
def test_run_writes_the_store_and_a_report_that_names_its_inputs(script, tmp_path):
    p = _write_inputs(tmp_path)
    out = tmp_path / "out" / "store.parquet"
    report_path = tmp_path / "out" / "report.json"
    r = script.run(p["cohort"], p["expression"], p["annotation"], p["joined"],
              out, report_path)
    store = pd.read_parquet(out)
    assert len(store) == r["n_patients"] == 8
    assert r["feature_columns"] == list(FEATURE_COLUMNS)
    for name in ("cohort", "expression", "annotation", "joined"):
        digest = hashlib.sha256(p[name].read_bytes()).hexdigest()
        assert r["inputs"][name]["sha256"] == digest
        assert r["inputs"][name]["bytes"] == p[name].stat().st_size
    assert r["output"]["sha256"] == hashlib.sha256(out.read_bytes()).hexdigest()
    text = report_path.read_text(encoding="ascii")
    assert json.loads(text) == r
    assert len(r["not_validated"]) == 2
    assert set(r["software"]) == {"gliostrata", "python", "pandas"}


def test_main_prints_the_counts_in_ascii(script, tmp_path, capsys):
    p = _write_inputs(tmp_path)
    script.main(["--cohort", str(p["cohort"]), "--expression", str(p["expression"]),
            "--annotation", str(p["annotation"]), "--joined", str(p["joined"]),
            "--out", str(tmp_path / "store.parquet"),
            "--report", str(tmp_path / "report.json")])
    out = capsys.readouterr().out
    assert out.isascii()
    assert "patients in the feature store       8" in out
    assert "extent_of_resection" in out
    assert "Nothing here is validated" in out
    assert (tmp_path / "store.parquet").is_file()
