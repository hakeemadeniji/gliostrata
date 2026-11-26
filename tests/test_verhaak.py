"""Tests for `gliostrata.profiling.verhaak` -- signature-based subtype scoring.

These tests can prove three things and deliberately do not pretend to prove a
fourth:

1. The *mechanics* are right -- gene-id resolution, version stripping,
   duplicate-symbol collapsing, missing-gene reporting, the compositional
   contract on the probability vector, determinism.
2. The *scoring* is right, by construction: synthetic matrices are built where
   exactly one subtype's signature genes are planted high in a known set of
   samples, and the classifier must recover that assignment.
3. The real 284-sample TCGA matrix is processed end to end and its canonical
   markers behave as expected.

What they cannot test is subtype *accuracy*: there are no subtype labels in
this project (clinical_cohort.parquet's `subtype` column is 100% NULL), so no
test here measures agreement with a ground-truth label, and none should be
read as doing so.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gliostrata.profiling import verhaak as V

DATA = Path(__file__).resolve().parents[1] / "data/raw/tcga_gbm/expression"
TPM_PATH = DATA / "tpm.parquet"
ANN_PATH = DATA / "gene_annotation.parquet"

PROB_COLUMNS = [f"subtype_prob_{s}" for s in V.SUBTYPES]


# ---------------------------------------------------------------------------
# Synthetic matrix builder
# ---------------------------------------------------------------------------
def _background_genes(n: int = 400) -> list[str]:
    """Filler gene symbols that are in no signature."""
    used = {g for genes in V.VERHAAK_SIGNATURES.values() for g in genes}
    out = []
    i = 0
    while len(out) < n:
        g = f"BGGENE{i}"
        if g not in used:
            out.append(g)
        i += 1
    return out


def synthetic_matrix(
    assignments: dict[str, str],
    *,
    high: float = 500.0,
    low: float = 10.0,
    seed: int = 0,
    signatures: dict[str, tuple[str, ...]] | None = None,
) -> pd.DataFrame:
    """Symbol-indexed TPM-like matrix with planted subtype signal.

    ``assignments`` maps sample id -> the subtype whose signature genes are
    planted high in that sample. Every other signature, and all background
    genes, sit at ``low``. Noise is small relative to the high/low gap so the
    planted truth is unambiguous.
    """
    sigs = V.VERHAAK_SIGNATURES if signatures is None else signatures
    all_sig_genes = sorted({g for genes in sigs.values() for g in genes})
    genes = all_sig_genes + _background_genes()
    rng = np.random.default_rng(seed)
    samples = list(assignments)
    mat = pd.DataFrame(
        rng.uniform(low * 0.8, low * 1.2, size=(len(genes), len(samples))),
        index=pd.Index(genes, name="gene_name"),
        columns=samples,
    )
    for sample, subtype in assignments.items():
        for g in sigs[subtype]:
            mat.loc[g, sample] = high * rng.uniform(0.9, 1.1)
    return mat


def _annotation_for(symbols, *, version: int = 7) -> pd.DataFrame:
    """Fake gene_annotation table giving each symbol a versioned Ensembl id."""
    return pd.DataFrame(
        {
            "gene_id": [f"ENSG{i:011d}.{version}" for i in range(len(symbols))],
            "gene_name": list(symbols),
            "gene_type": ["protein_coding"] * len(symbols),
        }
    )


# ---------------------------------------------------------------------------
# Gene identifier handling
# ---------------------------------------------------------------------------
def test_strip_ensembl_version():
    assert V.strip_ensembl_version(["ENSG00000146648.18"]) == ["ENSG00000146648"]
    # Already unversioned, and PAR_Y-style suffixes, are left alone.
    assert V.strip_ensembl_version(["ENSG00000146648"]) == ["ENSG00000146648"]
    assert V.strip_ensembl_version(["ENSG00000182378.14_PAR_Y"]) == [
        "ENSG00000182378.14_PAR_Y"
    ]
    assert V.strip_ensembl_version(["SYMBOL", None, np.nan]) == ["SYMBOL", "", ""]


def test_resolve_symbols_exact_and_version_stripped():
    ann = pd.DataFrame(
        {
            "gene_id": ["ENSG00000146648.18", "ENSG00000134121.15"],
            "gene_name": ["EGFR", "CHL1"],
            "gene_type": ["protein_coding", "protein_coding"],
        }
    )
    # Exact match.
    assert V.resolve_symbols(["ENSG00000146648.18"], ann) == ["EGFR"]
    # Different patch version still resolves via the stripped id.
    assert V.resolve_symbols(["ENSG00000146648.22"], ann) == ["EGFR"]
    # Unversioned input resolves too.
    assert V.resolve_symbols(["ENSG00000146648"], ann) == ["EGFR"]
    # Unknown ids come back empty -- never guessed.
    assert V.resolve_symbols(["ENSG99999999999.1"], ann) == [""]


def test_resolve_symbols_requires_schema():
    with pytest.raises(ValueError, match="missing required column"):
        V.resolve_symbols(["ENSG1"], pd.DataFrame({"gene_id": ["ENSG1"]}))


def test_to_symbol_matrix_maps_and_drops_unmapped():
    ann = pd.DataFrame(
        {
            "gene_id": ["ENSG00000146648.18"],
            "gene_name": ["EGFR"],
            "gene_type": ["protein_coding"],
        }
    )
    expr = pd.DataFrame(
        {"S1": [100.0, 5.0]},
        index=pd.Index(["ENSG00000146648.18", "ENSG00000000000.1"], name="gene_id"),
    )
    sym = V.to_symbol_matrix(expr, ann)
    assert list(sym.index) == ["EGFR"]
    assert sym.loc["EGFR", "S1"] == 100.0


def test_to_symbol_matrix_collapses_duplicates_by_max_mean():
    ann = pd.DataFrame(
        {
            "gene_id": ["ENSG1.1", "ENSG2.1"],
            "gene_name": ["DUP", "DUP"],
            "gene_type": ["protein_coding", "lncRNA"],
        }
    )
    expr = pd.DataFrame(
        {"S1": [1.0, 50.0], "S2": [3.0, 70.0]},
        index=pd.Index(["ENSG1.1", "ENSG2.1"], name="gene_id"),
    )
    sym = V.to_symbol_matrix(expr, ann)
    assert list(sym.index) == ["DUP"]
    # The higher-mean row (ENSG2) wins.
    assert sym.loc["DUP", "S1"] == 50.0


def test_to_symbol_matrix_duplicate_collapse_is_deterministic():
    ann = pd.DataFrame(
        {
            "gene_id": [f"ENSG{i}.1" for i in range(5)],
            "gene_name": ["DUP"] * 5,
            "gene_type": ["protein_coding"] * 5,
        }
    )
    expr = pd.DataFrame(
        {"S1": [1.0, 1.0, 9.0, 1.0, 1.0]},
        index=pd.Index([f"ENSG{i}.1" for i in range(5)], name="gene_id"),
    )
    first = V.to_symbol_matrix(expr, ann)
    for _ in range(5):
        pd.testing.assert_frame_equal(first, V.to_symbol_matrix(expr, ann))


# ---------------------------------------------------------------------------
# Signature mapping report
# ---------------------------------------------------------------------------
def test_map_signatures_reports_missing_genes_explicitly():
    mat = pd.DataFrame({"S1": [1.0, 2.0]}, index=["EGFR", "NES"])
    mappings = V.map_signatures(mat, {"Classical": V.VERHAAK_SIGNATURES["Classical"]})
    m = mappings["Classical"]
    assert m.found == ("EGFR", "NES")
    assert "AKT2" in m.missing
    assert m.n_requested == len(V.VERHAAK_SIGNATURES["Classical"])
    assert m.mapping_rate == pytest.approx(2 / m.n_requested)
    report = V.mapping_report(mappings)
    assert report.loc[0, "n_found"] == 2
    assert "AKT2" in report.loc[0, "missing_genes"]


def test_score_signatures_refuses_starved_signature():
    mat = pd.DataFrame({"S1": [1.0], "S2": [2.0]}, index=["EGFR"])
    with pytest.raises(ValueError, match="too few genes"):
        V.score_signatures(mat, {"Classical": V.VERHAAK_SIGNATURES["Classical"]})


def test_all_four_signatures_are_disjoint():
    seen: dict[str, str] = {}
    for name, genes in V.VERHAAK_SIGNATURES.items():
        assert len(set(genes)) == len(genes), f"{name} has duplicate genes"
        for g in genes:
            assert g not in seen, f"{g} is in both {seen.get(g)} and {name}"
            seen[g] = name


# ---------------------------------------------------------------------------
# Scoring correctness on planted signal
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("method", ["ssgsea", "zmean"])
def test_planted_subtype_wins(method):
    assignments = {}
    for i, st in enumerate(V.SUBTYPES):
        for rep in range(4):
            assignments[f"{st[:2].upper()}{rep}"] = st
    mat = synthetic_matrix(assignments)
    calls = V.classify_verhaak(mat, method=method)
    truth = pd.Series(assignments)
    got = calls.set_index("patient_id")["subtype"]
    assert (got.loc[truth.index] == truth).all(), got.to_dict()


def test_score_matrix_peaks_on_planted_signature():
    assignments = {"A": "Classical", "B": "Mesenchymal", "C": "Proneural", "D": "Neural"}
    mat = synthetic_matrix(assignments)
    sc = V.verhaak_score_matrix(mat)
    assert sc.scores.shape == (4, 4)
    assert sc.overall_mapping_rate == 1.0
    for sample, subtype in assignments.items():
        assert sc.scores.loc[sample].idxmax() == subtype


def test_probabilities_sum_to_one_and_match_argmax():
    assignments = {f"S{i}": V.SUBTYPES[i % 4] for i in range(12)}
    calls = V.classify_verhaak(synthetic_matrix(assignments))
    probs = calls[PROB_COLUMNS].to_numpy()
    assert np.allclose(probs.sum(axis=1), 1.0)
    assert probs.min() >= 0.0 and probs.max() <= 1.0
    argmax_label = [V.SUBTYPES[i] for i in probs.argmax(axis=1)]
    assert list(calls["subtype"]) == argmax_label
    assert np.allclose(calls["subtype_call_confidence"], probs.max(axis=1))


def test_output_columns_match_stub_contract():
    """The stub's downstream consumers expect these exact column names."""
    from gliostrata.profiling.subtype import SUBTYPES as STUB_SUBTYPES

    assert list(V.SUBTYPES) == list(STUB_SUBTYPES)
    assignments = {f"S{i}": V.SUBTYPES[i % 4] for i in range(8)}
    calls = V.classify_verhaak(synthetic_matrix(assignments))
    for col in ["patient_id", "subtype", "subtype_classifier_version", *PROB_COLUMNS]:
        assert col in calls.columns
    # And it is honest about not being the supervised model.
    assert bool(calls["subtype_is_supervised"].iloc[0]) is False


def test_temperature_flattens_but_does_not_change_argmax():
    assignments = {f"S{i}": V.SUBTYPES[i % 4] for i in range(12)}
    mat = synthetic_matrix(assignments)
    sharp = V.classify_verhaak(mat, temperature=0.5)
    flat = V.classify_verhaak(mat, temperature=4.0)
    assert list(sharp["subtype"]) == list(flat["subtype"])
    assert flat["subtype_call_confidence"].mean() < sharp["subtype_call_confidence"].mean()


def test_bad_temperature_and_method_rejected():
    mat = synthetic_matrix({f"S{i}": V.SUBTYPES[i % 4] for i in range(8)})
    with pytest.raises(ValueError, match="temperature"):
        V.classify_verhaak(mat, temperature=0.0)
    with pytest.raises(ValueError, match="unknown method"):
        V.classify_verhaak(mat, method="magic")


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------
def test_determinism_across_repeated_calls_and_seeds():
    """Fully deterministic: `seed` is accepted for parity and must not matter."""
    assignments = {f"S{i}": V.SUBTYPES[i % 4] for i in range(12)}
    mat = synthetic_matrix(assignments)
    a = V.classify_verhaak(mat, seed=42)
    b = V.classify_verhaak(mat, seed=42)
    c = V.classify_verhaak(mat, seed=12345)
    pd.testing.assert_frame_equal(a, b)
    pd.testing.assert_frame_equal(a, c)


def test_ssgsea_is_cohort_independent():
    """ssGSEA is per-sample: a sample's raw score must not depend on its peers."""
    assignments = {f"S{i}": V.SUBTYPES[i % 4] for i in range(12)}
    mat = synthetic_matrix(assignments)
    full = V.verhaak_score_matrix(mat).scores
    subset = V.verhaak_score_matrix(mat[["S0", "S1", "S2"]]).scores
    np.testing.assert_allclose(
        full.loc[["S0", "S1", "S2"]].to_numpy(), subset.to_numpy(), rtol=1e-12
    )


def test_row_shuffle_does_not_change_scores():
    assignments = {f"S{i}": V.SUBTYPES[i % 4] for i in range(8)}
    mat = synthetic_matrix(assignments)
    shuffled = mat.sample(frac=1.0, random_state=7)
    pd.testing.assert_frame_equal(
        V.verhaak_score_matrix(mat).scores,
        V.verhaak_score_matrix(shuffled).scores,
    )


# ---------------------------------------------------------------------------
# Degenerate inputs
# ---------------------------------------------------------------------------
def test_all_zero_sample_gets_uniform_probabilities():
    """An all-zero column is pure ties, so every ssGSEA score is exactly 0."""
    assignments = {f"S{i}": V.SUBTYPES[i % 4] for i in range(8)}
    mat = synthetic_matrix(assignments)
    mat["ZERO"] = 0.0
    sc = V.verhaak_score_matrix(mat)
    np.testing.assert_allclose(sc.scores.loc["ZERO"].to_numpy(), 0.0, atol=1e-12)
    calls = V.classify_verhaak(mat).set_index("patient_id")
    row = calls.loc["ZERO", PROB_COLUMNS].to_numpy(dtype=float)
    assert row.sum() == pytest.approx(1.0)
    # The zero sample is the least confident call in the cohort.
    assert calls.loc["ZERO", "subtype_call_confidence"] == pytest.approx(
        calls["subtype_call_confidence"].min()
    )


def test_constant_sample_also_scores_zero():
    assignments = {f"S{i}": V.SUBTYPES[i % 4] for i in range(8)}
    mat = synthetic_matrix(assignments)
    mat["FLAT"] = 42.0
    sc = V.verhaak_score_matrix(mat)
    np.testing.assert_allclose(sc.scores.loc["FLAT"].to_numpy(), 0.0, atol=1e-12)


def test_empty_sample_set_returns_empty_frame_with_contract_columns():
    mat = synthetic_matrix({f"S{i}": V.SUBTYPES[i % 4] for i in range(4)})
    empty = mat.iloc[:, :0]
    calls = V.classify_verhaak(empty)
    assert len(calls) == 0
    for col in ["patient_id", "subtype", *PROB_COLUMNS]:
        assert col in calls.columns
    dist = V.subtype_distribution(calls)
    assert list(dist["n"]) == [0, 0, 0, 0]


def test_single_sample_works_with_ssgsea_and_is_rejected_by_zmean():
    mat = synthetic_matrix({"ONLY": "Mesenchymal"})
    calls = V.classify_verhaak(mat, method="ssgsea")
    assert len(calls) == 1
    assert calls.loc[0, "subtype"] == "Mesenchymal"
    assert calls[PROB_COLUMNS].to_numpy().sum() == pytest.approx(1.0)
    with pytest.raises(ValueError, match="at least 2 samples"):
        V.classify_verhaak(mat, method="zmean")


def test_nan_input_is_rejected():
    mat = synthetic_matrix({f"S{i}": V.SUBTYPES[i % 4] for i in range(4)})
    mat.iloc[0, 0] = np.nan
    with pytest.raises(ValueError, match="NaN"):
        V.verhaak_score_matrix(mat)


# ---------------------------------------------------------------------------
# Label-free validation helpers
# ---------------------------------------------------------------------------
def test_marker_concordance_detects_planted_markers():
    assignments = {}
    for st in V.SUBTYPES:
        for rep in range(4):
            assignments[f"{st[:2].upper()}{rep}"] = st
    mat = synthetic_matrix(assignments)
    calls = V.classify_verhaak(mat)
    conc = V.marker_concordance(mat, calls)
    assert conc["concordant"].all(), conc.to_dict("records")
    assert (conc["fold_change"] > 1.0).all()


def test_marker_concordance_reports_absent_marker_rather_than_crashing():
    assignments = {f"S{i}": V.SUBTYPES[i % 4] for i in range(8)}
    mat = synthetic_matrix(assignments)
    calls = V.classify_verhaak(mat)
    conc = V.marker_concordance(mat, calls, markers={"NOT_A_GENE": "Classical"})
    assert len(conc) == 1
    assert conc.loc[0, "concordant"] is False or conc.loc[0, "concordant"] == False  # noqa: E712
    assert np.isnan(conc.loc[0, "mean_in_expected"])


def test_subtype_distribution_is_not_degenerate_on_planted_data():
    assignments = {}
    for st in V.SUBTYPES:
        for rep in range(5):
            assignments[f"{st[:2].upper()}{rep}"] = st
    calls = V.classify_verhaak(synthetic_matrix(assignments))
    dist = V.subtype_distribution(calls)
    assert (dist["n"] > 0).all()
    assert dist["fraction"].sum() == pytest.approx(1.0)


def test_silhouette_positive_on_well_separated_planted_data():
    assignments = {}
    for st in V.SUBTYPES:
        for rep in range(6):
            assignments[f"{st[:2].upper()}{rep}"] = st
    mat = synthetic_matrix(assignments)
    sc = V.verhaak_score_matrix(mat)
    calls = V.classify_verhaak(mat)
    sil = V.score_silhouette(sc.scores, calls["subtype"])
    assert sil > 0.3, sil


def test_silhouette_nan_with_one_class():
    scores = pd.DataFrame(np.ones((5, 4)), columns=list(V.SUBTYPES))
    assert np.isnan(V.score_silhouette(scores, ["Classical"] * 5))


def test_neural_is_flagged_as_contamination_axis():
    """Wang 2017: Neural is non-tumour admixture. The module must say so."""
    assert V.NEURAL_IS_CONTAMINATION_AXIS is True
    assert "contamination" in V.__doc__.lower()


# ---------------------------------------------------------------------------
# Real-data smoke test (skipped when the parquet is absent)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not (TPM_PATH.exists() and ANN_PATH.exists()),
    reason="TCGA-GBM expression parquet not present",
)
def test_real_tcga_gbm_end_to_end():
    tpm = pd.read_parquet(TPM_PATH)
    ann = pd.read_parquet(ANN_PATH)

    sc = V.verhaak_score_matrix(tpm, ann)
    # GENCODE v36 contains every gene in every signature.
    assert sc.overall_mapping_rate == 1.0, sc.report().to_dict("records")

    calls = V.classify_verhaak(tpm, ann)
    assert len(calls) == tpm.shape[1]
    assert np.allclose(calls[PROB_COLUMNS].to_numpy().sum(axis=1), 1.0)

    # Not degenerate: every subtype is called, and no subtype runs away.
    dist = V.subtype_distribution(calls)
    assert (dist["n"] >= 10).all(), dist.to_dict("records")
    assert dist["fraction"].max() < 0.55

    # Internal consistency: canonical markers are highest in their own group.
    conc = V.marker_concordance(tpm, calls, ann)
    assert conc["concordant"].all(), conc.to_dict("records")
    assert (conc["fold_change"] > 2.0).all(), conc.to_dict("records")

    # The partition is geometrically coherent in score space.
    assert V.score_silhouette(sc.scores, calls["subtype"]) > 0.05


@pytest.mark.skipif(
    not (TPM_PATH.exists() and ANN_PATH.exists()),
    reason="TCGA-GBM expression parquet not present",
)
def test_real_data_has_no_subtype_labels_to_validate_against():
    """Guards the central honesty claim: there is nothing to score accuracy on.

    If this test ever fails because labels appeared, that is the signal to
    build the supervised classifier the protocol actually specifies -- and to
    stop relying on this module's unsupervised scores.
    """
    cohort_path = TPM_PATH.parents[1] / "clinical_cohort.parquet"
    if not cohort_path.exists():
        pytest.skip("clinical cohort parquet not present")
    cohort = pd.read_parquet(cohort_path)
    assert "subtype" in cohort.columns
    assert cohort["subtype"].notna().sum() == 0, (
        "subtype labels now exist -- train the supervised model and measure T1"
    )
