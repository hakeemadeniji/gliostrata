"""Unit tests for scripts/idh_sensitivity_subtype.py (no network access).

What is at stake here
---------------------
This script re-reports the Verhaak subtype distribution with the IDH-mutant
expression samples removed. Three things can go wrong:

  1. The wrong samples are removed. Only idh_mutated == 1 is excluded; a
     sample with no IDH call stays.
  2. The aggregate hides the sample-level story. The committed classifier
     labels a sample relative to the other samples in the matrix, so a
     retained sample can change label when the cohort shrinks. The script
     must report how many did, or a small percentage-point change reads as
     more stable than it is.
  3. A sample disappears silently. The classifier's output is checked
     against the sample list it was given, every time.

The committed classifier is slow on 60k genes, so the logic tests inject a
stub classifier; the real classifier is exercised on the real matrix when
it is present.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "idh_sensitivity_subtype.py"

REAL_EXPRESSION = ROOT / "data/raw/tcga_gbm/expression/tpm.parquet"
REAL_ANNOTATION = ROOT / "data/raw/tcga_gbm/expression/gene_annotation.parquet"
REAL_JOINED = ROOT / "data/processed/tcga_gbm/molecular_calls_joined.parquet"
REAL_IDH = ROOT / "data/raw/tcga_gbm/mutations/idh_calls.parquet"

#: Called-subtype distribution of the committed classifier on the 284
#: expression samples (Classical / Mesenchymal / Proneural / Neural). A
#: canary: if this trips, the expression matrix, the annotation or the
#: classifier changed -- confirm the change was intended, never widen the
#: tolerance to keep it green.
COMMITTED_PERCENT = {
    "Classical": 31.0,
    "Mesenchymal": 28.9,
    "Proneural": 23.6,
    "Neural": 16.5,
}


def _load_module():
    """Load the script by path (scripts/ is not a package)."""
    assert SCRIPT_PATH.is_file(), f"script not found: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location(
        "gliostrata_idh_sensitivity_subtype_script", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def S():
    return _load_module()


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
def _expression(samples):
    return pd.DataFrame(
        [[1.0] * len(samples), [2.0] * len(samples)], columns=list(samples),
        index=pd.Index(["ENSG_A", "ENSG_B"], name="gene_id"),
    )


def _joined(ids, mutant_ids, no_call_ids=()):
    status = pd.Series(0, index=pd.Index(ids), dtype="Int64")
    status.loc[list(mutant_ids)] = 1
    status.loc[list(no_call_ids)] = pd.NA
    return pd.DataFrame({"patient_id": status.index, "idh_mutated": status.to_numpy()})


def _stub(labels, flip=None):
    """Classifier stub: a fixed label per sample, optionally cohort-dependent.

    `flip` maps a sample to the label it gets whenever every `trigger`
    sample is still present in the matrix. That reproduces, in miniature,
    the within-cohort dependence of the real classifier.
    """
    seen = []

    def classify(expression, annotation):
        cols = [str(c) for c in expression.columns]
        seen.append(list(cols))
        out = []
        for s in cols:
            lab = labels[s]
            if flip and s in flip and set(flip[s]["triggers"]) <= set(cols):
                lab = flip[s]["to"]
            out.append(lab)
        return pd.DataFrame({"patient_id": cols, "subtype": out})

    classify.seen = seen
    return classify


# ---------------------------------------------------------------------------
# 1. The right samples are excluded
# ---------------------------------------------------------------------------
_LABELS = {
    "A1": "Classical", "A2": "Classical", "A3": "Classical",
    "B1": "Mesenchymal", "B2": "Mesenchymal",
    "C1": "Proneural",
    "M1": "Neural", "M2": "Proneural",       # IDH mutants
    "N1": "Neural", "N2": "Neural",          # no IDH call
}


def _base_joined():
    return _joined(_LABELS, mutant_ids=["M1", "M2"], no_call_ids=["N1", "N2"])


def test_mutants_are_removed_before_the_second_classification(S):
    stub = _stub(_LABELS)
    r = S.analyze(_expression(_LABELS), None, _base_joined(), classify=stub)

    assert stub.seen[0] == list(_LABELS), "the full run sees every sample"
    assert stub.seen[1] == [s for s in _LABELS if s not in ("M1", "M2")]
    # No-call samples are NOT removed: an absent call proves nothing.
    assert "N1" in stub.seen[1] and "N2" in stub.seen[1]
    assert r["n_expression_full"] == 10
    assert r["n_idh_mutant_excluded"] == 2
    assert r["n_expression_without_idh"] == 8
    assert r["idh_mutant_samples_excluded"] == ["M1", "M2"]


def test_a_missing_call_is_not_treated_as_mutant(S):
    stub = _stub(_LABELS)
    joined = _joined(_LABELS, mutant_ids=[], no_call_ids=["M1", "M2"])
    r = S.analyze(_expression(_LABELS), None, joined, classify=stub)
    assert r["n_idh_mutant_excluded"] == 0
    assert r["distribution_without_idh"] == r["distribution_full"]


def test_distribution_counts_and_percentages(S):
    stub = _stub(_LABELS)
    r = S.analyze(_expression(_LABELS), None, _base_joined(), classify=stub)
    d = r["distribution_full"]
    assert list(d) == list(S.SUBTYPES)
    assert d["Classical"] == {"n": 3, "percent": 30.0}
    assert d["Mesenchymal"] == {"n": 2, "percent": 20.0}
    assert d["Proneural"] == {"n": 2, "percent": 20.0}
    assert d["Neural"] == {"n": 3, "percent": 30.0}
    assert sum(v["n"] for v in d.values()) == 10
    assert r["distribution_without_idh"]["Neural"] == {"n": 2, "percent": 25.0}
    assert r["distribution_without_idh"]["Proneural"] == {"n": 1, "percent": 12.5}


def test_pp_change_is_signed_and_its_maximum_is_absolute(S):
    stub = _stub(_LABELS)
    r = S.analyze(_expression(_LABELS), None, _base_joined(), classify=stub)
    change = r["pp_change_per_subtype"]
    assert change["Classical"] == pytest.approx(7.5), "3/8 vs 3/10"
    assert change["Mesenchymal"] == pytest.approx(5.0)
    assert change["Proneural"] == pytest.approx(-7.5)
    assert change["Neural"] == pytest.approx(-5.0)
    # Percentage-point changes must net to zero across the subtypes.
    assert sum(change.values()) == pytest.approx(0.0, abs=0.05)
    assert r["max_pp_change"] == pytest.approx(7.5)


def test_sample_level_relabeling_is_reported(S):
    """A retained sample whose label depends on the cohort must be counted.

    C1 is Classical while M1 is in the matrix and Proneural once M1 is
    gone, although C1 itself was never excluded.
    """
    stub = _stub(_LABELS, flip={"C1": {"triggers": ["M1"], "to": "Classical"}})
    r = S.analyze(_expression(_LABELS), None, _base_joined(), classify=stub)
    assert r["n_relabelled_among_retained"] == 1

    # Full matrix: C1 counts as Classical.
    assert r["distribution_full"]["Classical"] == {"n": 4, "percent": 40.0}
    assert r["distribution_full"]["Proneural"] == {"n": 1, "percent": 10.0}
    # Reduced matrix: M2 is gone, and C1 now counts as Proneural.
    assert r["distribution_without_idh"]["Classical"] == {"n": 3, "percent": 37.5}
    assert r["distribution_without_idh"]["Proneural"] == {"n": 1, "percent": 12.5}
    assert r["pp_change_per_subtype"]["Proneural"] == pytest.approx(2.5)

    # With every retained sample keeping its full-cohort label (C1 stays
    # Classical), the only movement is the removal of M2: Proneural 0/8
    # against 1/10.
    frozen = r["distribution_retained_with_full_cohort_labels"]
    assert frozen["Classical"]["n"] == 4 and frozen["Proneural"]["n"] == 0
    assert r["pp_change_per_subtype_with_full_cohort_labels"]["Proneural"] == (
        pytest.approx(-10.0)
    )
    gone = r["distribution_excluded_with_full_cohort_labels"]
    assert gone["Neural"]["n"] == 1 and gone["Proneural"]["n"] == 1
    for subtype in S.SUBTYPES:
        assert frozen[subtype]["n"] + gone[subtype]["n"] == (
            r["distribution_full"][subtype]["n"]
        )


def test_no_relabeling_when_labels_are_cohort_independent(S):
    stub = _stub(_LABELS)
    r = S.analyze(_expression(_LABELS), None, _base_joined(), classify=stub)
    assert r["n_relabelled_among_retained"] == 0
    assert r["pp_change_per_subtype_with_full_cohort_labels"] == (
        r["pp_change_per_subtype"]
    )


# ---------------------------------------------------------------------------
# 2. No sample disappears, and nothing is invented
# ---------------------------------------------------------------------------
def test_a_sample_outside_the_joined_table_is_kept_and_listed(S):
    labels = dict(_LABELS, X1="Neural")
    stub = _stub(labels)
    r = S.analyze(_expression(labels), None, _base_joined(), classify=stub)
    assert "X1" in stub.seen[1], "no IDH status is not a reason to drop a sample"
    assert r["n_expression_full"] == 11
    assert r["n_retained_without_joined_row"] == 1
    assert r["retained_without_joined_row"] == ["X1"]
    assert r["distribution_full"]["Neural"]["n"] == 4


def test_classifier_output_is_checked_against_the_samples(S):
    labels = dict(_LABELS, X1="Neural")

    def drops_one(expression, annotation):
        return _stub(labels)(expression, annotation).iloc[1:]

    with pytest.raises(ValueError, match="does not match the samples"):
        S.analyze(_expression(labels), None, _base_joined(), classify=drops_one)

    def invents_a_label(expression, annotation):
        out = _stub(labels)(expression, annotation)
        out["subtype"] = out["subtype"].replace({"Neural": "NotASubtype"})
        return out

    with pytest.raises(ValueError, match="unknown subtype labels"):
        S.analyze(_expression(labels), None, _base_joined(), classify=invents_a_label)


def test_bad_inputs_are_refused(S):
    stub = _stub(_LABELS)
    with pytest.raises(ValueError, match="required column"):
        S.analyze(
            _expression(_LABELS), None,
            _base_joined().drop(columns=["idh_mutated"]), classify=stub,
        )
    with pytest.raises(ValueError, match="leaves nothing"):
        S.analyze(
            _expression(_LABELS), None,
            _joined(_LABELS, mutant_ids=list(_LABELS)), classify=stub,
        )
    with pytest.raises(ValueError, match="duplicate sample columns"):
        S.analyze(_expression(["A1", "A1"]), None, _base_joined(), classify=stub)


def test_result_contract(S):
    stub = _stub(_LABELS)
    r = S.analyze(_expression(_LABELS), None, _base_joined(), classify=stub)
    assert tuple(list(r)[: len(S.RESULT_KEYS)]) == S.RESULT_KEYS
    assert type(r["max_pp_change"]) is float
    biggest = max(abs(v) for v in r["pp_change_per_subtype"].values())
    assert r["max_pp_change"] == biggest
    assert json.loads(json.dumps(r)) == r


def test_run_writes_the_result_and_is_idempotent(S, tmp_path):
    _base_joined().to_parquet(tmp_path / "joined.parquet", index=False)
    _expression(_LABELS).to_parquet(tmp_path / "tpm.parquet")
    pd.DataFrame({"gene_id": ["ENSG_A"], "gene_name": ["A"]}).to_parquet(
        tmp_path / "ann.parquet"
    )
    paths = (tmp_path / "tpm.parquet", tmp_path / "ann.parquet",
             tmp_path / "joined.parquet")
    out = tmp_path / "nested" / "result.json"
    r = S.run(*paths, out, classify=_stub(_LABELS))
    text = out.read_text(encoding="ascii")
    assert json.loads(text) == r
    assert r["n_expression_without_idh"] == 8

    S.run(*paths, out, classify=_stub(_LABELS))
    assert out.read_text(encoding="ascii") == text


# ---------------------------------------------------------------------------
# 3. The committed classifier on the real matrix (skipped when absent)
# ---------------------------------------------------------------------------
real_data = pytest.mark.skipif(
    not all(p.exists() for p in (REAL_EXPRESSION, REAL_ANNOTATION, REAL_JOINED)),
    reason="TCGA-GBM expression or joined molecular calls not present",
)


@real_data
def test_real_classifier_reproduces_the_committed_distribution(S):
    from gliostrata.profiling.verhaak import classify_verhaak

    expression = pd.read_parquet(REAL_EXPRESSION)
    annotation = pd.read_parquet(REAL_ANNOTATION)
    joined = pd.read_parquet(REAL_JOINED)
    r = S.analyze(expression, annotation, joined, classify=classify_verhaak)

    assert r["n_expression_full"] == 284
    assert r["n_idh_mutant_excluded"] == 16
    assert r["n_expression_without_idh"] == 268

    for subtype, expected in COMMITTED_PERCENT.items():
        got = r["distribution_full"][subtype]["percent"]
        assert got == pytest.approx(expected, abs=0.5), (
            f"{subtype}: committed split is {expected}%, got {got}%"
        )

    for key, n in (
        ("distribution_full", r["n_expression_full"]),
        ("distribution_without_idh", r["n_expression_without_idh"]),
    ):
        assert sum(v["n"] for v in r[key].values()) == n
        assert sum(v["percent"] for v in r[key].values()) == pytest.approx(
            100.0, abs=0.1
        )

    # The frozen-label counts are a partition of the full run: for every
    # subtype, retained + excluded must rebuild the full count exactly.
    frozen = r["distribution_retained_with_full_cohort_labels"]
    gone = r["distribution_excluded_with_full_cohort_labels"]
    for subtype in S.SUBTYPES:
        assert frozen[subtype]["n"] + gone[subtype]["n"] == (
            r["distribution_full"][subtype]["n"]
        ), subtype
    assert sum(v["n"] for v in gone.values()) == r["n_idh_mutant_excluded"]
    assert 0 <= r["n_relabelled_among_retained"] <= r["n_expression_without_idh"]
    assert r["max_pp_change"] == max(
        abs(v) for v in r["pp_change_per_subtype"].values()
    )


@real_data
@pytest.mark.skipif(not REAL_IDH.exists(), reason="idh_calls.parquet not present")
def test_real_sample_without_a_joined_row_is_not_idh_mutant(S):
    """TCGA-28-2510 has expression but no clinical-cohort row.

    The script keeps it because the joined table holds no IDH status for
    it. That is only correct if it is not in fact a mutant, which the raw
    call table settles.
    """
    import pyarrow.parquet as pq

    joined_ids = set(
        pd.read_parquet(REAL_JOINED, columns=["patient_id"])["patient_id"].astype(str)
    )
    expression_ids = [
        c for c in pq.ParquetFile(REAL_EXPRESSION).schema.names
        if str(c).startswith("TCGA-")
    ]
    orphans = [s for s in expression_ids if s not in joined_ids]
    assert orphans == ["TCGA-28-2510"], (
        "the set of expression patients outside the clinical cohort changed; "
        "re-check the keep-and-list rule in scripts/idh_sensitivity_subtype.py"
    )

    calls = pd.read_parquet(REAL_IDH, columns=["patient_id", "idh_mutated"])
    status = calls.set_index("patient_id")["idh_mutated"].loc[orphans]
    assert status.notna().all(), "IDH status unknown: the keep decision is unverified"
    assert not (status == 1).any(), (
        "an expression sample outside the cohort is IDH-mutant; the script "
        "keeps it because the joined table has no row for it"
    )
