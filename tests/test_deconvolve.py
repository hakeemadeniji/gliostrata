"""Tests for `gliostrata.profiling.deconvolve` -- Neftel state and TME scoring.

What these tests establish, and what they deliberately do not
------------------------------------------------------------
They establish that the *mechanics* are correct: gene-id resolution and
Ensembl version-suffix stripping, explicit reporting of missing signature
genes, the compositional contract (every row in [0, 1] and summing to 1),
determinism, correct behaviour on degenerate samples, and -- by planting known
marker signal in synthetic matrices -- that the intended state or TME
compartment actually wins.

They do NOT and cannot establish that the output approximates real cell-state
or cell-type proportions. The module normalises marker-enrichment scores to
sum to 1; that is not deconvolution, there is no reference profile and no
mixing model, and no ground truth exists in this project to compare against.
Protocol test T1 (state-deconvolution Pearson r >= 0.75 against ground truth)
is therefore NOT addressed by this module or by this test file, and nothing
here should be cited as evidence for it. Validating T1 needs either a
CIBERSORTx run on these samples or pseudo-bulk mixtures built from the Neftel
atlas with known mixing weights.

One test (`test_rescaling_a_sample_does_not_change_its_fractions`) and one
real-data test (`test_real_data_fractions_shift_with_misc_rna_batch`) exist
specifically because this cohort carries a severe library-prep batch effect:
a large subset of samples puts the majority of their TPM mass in `misc_RNA`
transcripts. Marker-enrichment scoring reads *relative* abundance, so it is
exposed to exactly that kind of distortion. Those two tests pin down which
distortions the module is immune to and which it is not.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gliostrata.profiling import deconvolve as D
from gliostrata.profiling.verhaak import map_signatures

DATA = Path(__file__).resolve().parents[1] / "data/raw/tcga_gbm/expression"
TPM_PATH = DATA / "tpm.parquet"
ANN_PATH = DATA / "gene_annotation.parquet"

HAVE_REAL = TPM_PATH.exists() and ANN_PATH.exists()
needs_real = pytest.mark.skipif(not HAVE_REAL, reason="TCGA-GBM expression parquet not present")

STATES = list(D.STATES)
TME = list(D.TME_FRACTIONS)
STATE_NAMES = list(D.NEFTEL_SIGNATURES)
TME_NAMES = list(D.TME_SIGNATURES)


# ---------------------------------------------------------------------------
# Synthetic matrix builders
# ---------------------------------------------------------------------------
def _background_genes(signatures: dict[str, tuple[str, ...]], n: int = 400) -> list[str]:
    """Filler symbols that are in no panel used by this module."""
    used = {g for genes in signatures.values() for g in genes}
    used |= {g for genes in D.NEFTEL_SIGNATURES.values() for g in genes}
    used |= {g for genes in D.TME_SIGNATURES.values() for g in genes}
    used |= set(D.BROAD_IMMUNE_PANEL)
    out: list[str] = []
    i = 0
    while len(out) < n:
        g = f"BGGENE{i}"
        if g not in used:
            out.append(g)
        i += 1
    return out


def synthetic_matrix(
    signatures: dict[str, tuple[str, ...]],
    assignments: dict[str, str],
    *,
    high: float = 500.0,
    low: float = 10.0,
    seed: int = 0,
    extra_genes: tuple[str, ...] = (),
) -> pd.DataFrame:
    """Symbol-indexed TPM-like matrix with one panel planted high per sample.

    ``assignments`` maps sample id -> the panel name whose genes are planted
    high in that sample. All other panels and all background genes sit near
    ``low``. The high/low gap is large relative to the noise, so the planted
    truth is unambiguous and the expected argmax is not a coin flip.
    """
    all_sig_genes = sorted({g for genes in signatures.values() for g in genes})
    genes = all_sig_genes + list(extra_genes) + _background_genes(signatures)
    genes = list(dict.fromkeys(genes))
    rng = np.random.default_rng(seed)
    samples = list(assignments)
    mat = pd.DataFrame(
        rng.uniform(low * 0.8, low * 1.2, size=(len(genes), len(samples))),
        index=pd.Index(genes, name="gene_name"),
        columns=samples,
    )
    for sample, panel in assignments.items():
        for g in signatures[panel]:
            mat.loc[g, sample] = high * rng.uniform(0.9, 1.1)
    return mat


def state_matrix(reps: int = 4, **kw) -> tuple[pd.DataFrame, dict[str, str]]:
    assignments = {
        f"{name.split('-')[0]}{rep}": name
        for name in STATE_NAMES
        for rep in range(reps)
    }
    return synthetic_matrix(D.NEFTEL_SIGNATURES, assignments, **kw), assignments


def tme_matrix(reps: int = 4, **kw) -> tuple[pd.DataFrame, dict[str, str]]:
    assignments = {
        f"{name[:3].upper()}{rep}": name for name in TME_NAMES for rep in range(reps)
    }
    return synthetic_matrix(D.TME_SIGNATURES, assignments, **kw), assignments


def _annotation_for(symbols, *, version: int = 7) -> pd.DataFrame:
    """Fake gene_annotation table giving each symbol a versioned Ensembl id."""
    return pd.DataFrame(
        {
            "gene_id": [f"ENSG{i:011d}.{version}" for i in range(len(symbols))],
            "gene_name": list(symbols),
            "gene_type": ["protein_coding"] * len(symbols),
        }
    )


def _as_ensembl(mat: pd.DataFrame, *, version: int = 7) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Re-index a symbol matrix by versioned Ensembl ids and return the annotation."""
    ann = _annotation_for(list(mat.index), version=version)
    out = mat.copy()
    out.index = pd.Index(ann["gene_id"].tolist(), name="gene_id")
    return out, ann


# ---------------------------------------------------------------------------
# Panel hygiene -- a gene must not vote twice
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "signatures,label",
    [(D.NEFTEL_SIGNATURES, "Neftel"), (D.TME_SIGNATURES, "TME")],
)
def test_panels_are_internally_disjoint_and_duplicate_free(signatures, label):
    seen: dict[str, str] = {}
    for name, genes in signatures.items():
        assert len(set(genes)) == len(genes), f"{label}/{name} has duplicate genes"
        for g in genes:
            assert g not in seen, f"{label}: {g} is in both {seen[g]} and {name}"
            seen[g] = name


def test_panel_names_map_onto_the_downstream_column_contract():
    from gliostrata.profiling.states import STATES as STUB_STATES
    from gliostrata.profiling.tme import TME_FRACTIONS as STUB_TME

    assert STATES == list(STUB_STATES)
    assert TME == list(STUB_TME)
    assert len(D.NEFTEL_SIGNATURES) == len(STATES)
    assert len(D.TME_SIGNATURES) == len(TME)


def test_broad_immune_panel_is_nearly_independent_of_the_scoring_panel():
    """A cross-check panel that shared most of its genes would be tautological."""
    overlap = set(D.BROAD_IMMUNE_PANEL) & set(D.TME_SIGNATURES["immune"])
    assert overlap == {"PTPRC"}, overlap
    assert len(D.BROAD_IMMUNE_PANEL) >= 15


# ---------------------------------------------------------------------------
# Gene identifier handling
# ---------------------------------------------------------------------------
def test_ensembl_version_suffix_is_stripped_when_resolving():
    """ENSG00000146648.18 -> EGFR, and a different patch version still resolves."""
    ann = pd.DataFrame(
        {
            "gene_id": ["ENSG00000146648.18", "ENSG00000073756.12"],
            "gene_name": ["EGFR", "PTGS2"],
            "gene_type": ["protein_coding", "protein_coding"],
        }
    )
    expr = pd.DataFrame(
        {"S1": [100.0, 7.0], "S2": [120.0, 9.0]},
        # Deliberately a DIFFERENT patch version than the annotation.
        index=pd.Index(["ENSG00000146648.22", "ENSG00000073756.12"], name="gene_id"),
    )
    sym = D.to_symbol_matrix(expr, ann)
    assert list(sym.index) == ["EGFR", "PTGS2"]
    assert sym.loc["EGFR", "S1"] == 100.0


def test_scoring_through_ensembl_ids_matches_scoring_through_symbols():
    mat, assignments = state_matrix()
    by_symbol = D.deconvolve_states_expression(mat)
    ens, ann = _as_ensembl(mat)
    by_ensembl = D.deconvolve_states_expression(ens, ann)
    pd.testing.assert_frame_equal(by_symbol, by_ensembl)


def test_missing_panel_genes_are_reported_not_silently_ignored():
    mat, _ = state_matrix()
    # Remove three AC-like genes from the matrix entirely.
    dropped = ["GFAP", "AQP4", "SLC1A2"]
    trimmed = mat.drop(index=dropped)
    res = D.state_score_matrix(trimmed)
    ac = res.mappings["AC-like"]
    assert set(dropped) <= set(ac.missing)
    assert ac.n_found == ac.n_requested - 3
    assert ac.mapping_rate == pytest.approx(ac.n_found / ac.n_requested)
    report = res.report().set_index("signature")
    for g in dropped:
        assert g in report.loc["AC-like", "missing_genes"]
    # Other panels are untouched.
    assert res.mappings["MES-like"].mapping_rate == 1.0
    assert res.overall_mapping_rate < 1.0


def test_panel_starved_of_genes_raises_rather_than_degrading_silently():
    mat, _ = state_matrix()
    # Leave AC-like with a single gene present.
    keep = ["GFAP"]
    drop = [g for g in D.NEFTEL_SIGNATURES["AC-like"] if g not in keep]
    with pytest.raises(ValueError, match="too few genes"):
        D.deconvolve_states_expression(mat.drop(index=drop))


def test_unmapped_gene_ids_are_dropped_and_show_up_as_missing():
    """Unmapped ids never reach scoring; the mapping report names every one.

    ``map_signatures`` is the reporting path and allows a panel with zero genes
    found, so the full report is visible. ``score_signatures`` then REFUSES to
    score such a panel rather than returning a degraded number -- both
    behaviours are asserted, because either alone is the wrong contract.
    """
    mat, _ = state_matrix()
    ens, ann = _as_ensembl(mat)
    # Blank out the annotation rows for the whole MES-like panel.
    mes = set(D.NEFTEL_SIGNATURES["MES-like"])
    ann = ann[~ann["gene_name"].isin(mes)].reset_index(drop=True)

    sym = D.to_symbol_matrix(ens, ann)
    assert not (set(sym.index) & mes), "unmapped ids must not reach the symbol matrix"

    # map_signatures lives in verhaak.py (shared scoring machinery).
    mappings = map_signatures(sym, D.NEFTEL_SIGNATURES)
    assert mappings["MES-like"].n_found == 0
    assert set(mappings["MES-like"].missing) == mes
    assert mappings["MES-like"].mapping_rate == 0.0
    # Other panels are untouched.
    assert mappings["AC-like"].mapping_rate == 1.0

    report = D.mapping_report(mappings).set_index("signature")
    reported = report.loc["MES-like", "missing_genes"].split(",")
    # Every missing gene is named, in the panel's own order -- no guessing.
    assert reported == [g for g in D.NEFTEL_SIGNATURES["MES-like"]]
    assert report.loc["MES-like", "n_found"] == 0

    # And scoring is refused rather than silently degraded.
    with pytest.raises(ValueError, match="too few genes"):
        D.state_score_matrix(ens, ann)


# ---------------------------------------------------------------------------
# Planted signal: the right compartment must win
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("score_method", ["ssgsea", "zmean"])
def test_planted_state_wins(score_method):
    mat, assignments = state_matrix()
    frac = D.deconvolve_states_expression(mat, score_method=score_method).set_index("patient_id")
    for sample, panel in assignments.items():
        expected = D._STATE_TO_COLUMN[panel]
        got = frac.loc[sample, STATES].astype(float).idxmax()
        assert got == expected, f"{sample}: planted {panel}, got {got}"


@pytest.mark.parametrize("score_method", ["ssgsea", "zmean"])
def test_planted_tme_compartment_wins(score_method):
    mat, assignments = tme_matrix()
    frac = D.deconvolve_tme_expression(mat, score_method=score_method).set_index("patient_id")
    for sample, panel in assignments.items():
        expected = D._TME_TO_COLUMN[panel]
        got = frac.loc[sample, TME].astype(float).idxmax()
        assert got == expected, f"{sample}: planted {panel}, got {got}"


def test_gfap_high_sample_gets_the_highest_ac_fraction():
    """The brief's explicit check, on synthetic data where GFAP is not saturated."""
    mat, assignments = state_matrix()
    frac = D.deconvolve_states_expression(mat).set_index("patient_id")
    ac_samples = [s for s, p in assignments.items() if p == "AC-like"]
    others = [s for s in frac.index if s not in ac_samples]
    assert mat.loc["GFAP", ac_samples].min() > mat.loc["GFAP", others].max()
    assert frac.loc[ac_samples, "state_ac"].min() > frac.loc[others, "state_ac"].max()


def test_ptprc_high_sample_gets_the_highest_immune_fraction():
    """The brief's explicit check for the TME side."""
    mat, assignments = tme_matrix()
    frac = D.deconvolve_tme_expression(mat).set_index("patient_id")
    imm = [s for s, p in assignments.items() if p == "immune"]
    others = [s for s in frac.index if s not in imm]
    assert mat.loc["PTPRC", imm].min() > mat.loc["PTPRC", others].max()
    assert frac.loc[imm, "tme_immune"].min() > frac.loc[others, "tme_immune"].max()


def test_raw_score_matrix_peaks_on_the_planted_panel():
    mat, assignments = state_matrix(reps=2)
    scores = D.state_score_matrix(mat).scores
    for sample, panel in assignments.items():
        assert scores.loc[sample].idxmax() == panel


# ---------------------------------------------------------------------------
# Compositional contract
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("fraction_method", ["softmax_z", "proportional"])
@pytest.mark.parametrize(
    "fn,cols,builder",
    [
        (D.deconvolve_states_expression, STATES, state_matrix),
        (D.deconvolve_tme_expression, TME, tme_matrix),
    ],
)
def test_fractions_are_in_range_and_sum_to_one(fn, cols, builder, fraction_method):
    mat, _ = builder()
    out = fn(mat, fraction_method=fraction_method)
    summary = D.check_fractions(out, tuple(cols))
    assert summary["n"] == mat.shape[1]
    assert summary["min"] >= 0.0
    assert summary["max"] <= 1.0
    assert summary["max_abs_sum_error"] < 1e-9
    assert np.allclose(out[cols].to_numpy().sum(axis=1), 1.0)


def test_check_fractions_rejects_out_of_range_values():
    bad = pd.DataFrame({"patient_id": ["A"], **{c: [0.0] for c in STATES}})
    bad.loc[0, "state_npc"] = 1.5
    bad.loc[0, "state_opc"] = -0.5
    with pytest.raises(ValueError, match=r"out of \[0, 1\]"):
        D.check_fractions(bad, D.STATES)


def test_check_fractions_rejects_rows_that_do_not_sum_to_one():
    bad = pd.DataFrame({"patient_id": ["A"], **{c: [0.1] for c in STATES}})
    with pytest.raises(ValueError, match="do not sum to 1"):
        D.check_fractions(bad, D.STATES)


def test_check_fractions_rejects_nan_and_missing_columns():
    nan_frame = pd.DataFrame({"patient_id": ["A"], **{c: [0.25] for c in STATES}})
    nan_frame.loc[0, "state_mes"] = np.nan
    with pytest.raises(ValueError, match="NaN"):
        D.check_fractions(nan_frame, D.STATES)
    with pytest.raises(ValueError, match="missing fraction columns"):
        D.check_fractions(pd.DataFrame({"patient_id": ["A"]}), D.STATES)


def test_output_columns_match_the_stub_contract_exactly():
    mat, _ = state_matrix(reps=2)
    states = D.deconvolve_states_expression(mat)
    for col in ["patient_id", *STATES, "state_deconvolution_method", "state_deconvolution_version"]:
        assert col in states.columns

    mat2, _ = tme_matrix(reps=2)
    tme = D.deconvolve_tme_expression(mat2)
    for col in ["patient_id", *TME, "tme_deconvolution_method", "tme_deconvolution_version"]:
        assert col in tme.columns


def test_module_declares_it_is_not_true_deconvolution():
    """The honesty requirement, enforced mechanically.

    If someone later swaps in a real deconvolution method they must flip these
    flags deliberately, and this test is where they will be reminded to.
    """
    mat, _ = state_matrix(reps=2)
    out = D.deconvolve_states_expression(mat)
    assert bool(out["is_true_deconvolution"].iloc[0]) is False
    assert out["state_deconvolution_method"].iloc[0] == D.METHOD_NAME
    assert "not-deconvolution" in D.METHOD_VERSION
    doc = D.deconvolve_states_expression.__doc__ or ""
    assert "not cell-state proportions" in doc
    assert "T1" in doc
    assert "not cell-type proportions" in (D.deconvolve_tme_expression.__doc__ or "")
    module_doc = (D.__doc__ or "").lower()
    assert "not cell-type proportions" in module_doc
    assert "cibersortx" in module_doc


def test_scores_to_fractions_rejects_bad_arguments():
    scores = pd.DataFrame(np.ones((3, 4)), columns=STATE_NAMES)
    with pytest.raises(ValueError, match="unknown method"):
        D.scores_to_fractions(scores, method="nnls")
    with pytest.raises(ValueError, match="temperature"):
        D.scores_to_fractions(scores, temperature=0.0)


def test_proportional_falls_back_to_uniform_when_all_scores_nonpositive():
    scores = pd.DataFrame(-np.ones((2, 4)), columns=STATE_NAMES, index=["A", "B"])
    frac = D.scores_to_fractions(scores, method="proportional")
    assert np.allclose(frac.to_numpy(), 0.25)


def test_temperature_flattens_fractions_without_changing_the_argmax():
    mat, _ = state_matrix()
    sharp = D.deconvolve_states_expression(mat, temperature=0.5).set_index("patient_id")
    flat = D.deconvolve_states_expression(mat, temperature=4.0).set_index("patient_id")
    assert list(sharp[STATES].idxmax(axis=1)) == list(flat[STATES].idxmax(axis=1))
    assert flat[STATES].max(axis=1).mean() < sharp[STATES].max(axis=1).mean()


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "fn,builder", [(D.deconvolve_states_expression, state_matrix), (D.deconvolve_tme_expression, tme_matrix)]
)
def test_determinism_and_seed_is_inert(fn, builder):
    """Deterministic by construction: `seed` exists for parity and must not matter."""
    mat, _ = builder()
    a = fn(mat, seed=42)
    b = fn(mat, seed=42)
    c = fn(mat, seed=999_999)
    pd.testing.assert_frame_equal(a, b)
    pd.testing.assert_frame_equal(a, c)


def test_gene_row_order_does_not_change_the_result():
    mat, _ = state_matrix()
    shuffled = mat.sample(frac=1.0, random_state=11)
    pd.testing.assert_frame_equal(
        D.deconvolve_states_expression(mat),
        D.deconvolve_states_expression(shuffled),
    )


def test_duplicate_symbols_collapse_deterministically():
    """GENCODE v36 has 1233 duplicated symbols; the max-mean rule must be stable."""
    mat, _ = state_matrix(reps=2)
    ens, ann = _as_ensembl(mat)
    # Add five extra rows all annotated as GFAP, one of them the highest-mean row.
    extra_ids = [f"ENSGDUP{i:05d}.3" for i in range(5)]
    extra = pd.DataFrame(
        np.tile(np.array([1.0, 1.0, 9999.0, 1.0, 1.0])[:, None], (1, ens.shape[1])),
        index=pd.Index(extra_ids, name="gene_id"),
        columns=ens.columns,
    )
    ens2 = pd.concat([ens, extra])
    ann2 = pd.concat(
        [ann, pd.DataFrame({"gene_id": extra_ids, "gene_name": ["GFAP"] * 5,
                            "gene_type": ["protein_coding"] * 5})],
        ignore_index=True,
    )
    first = D.deconvolve_states_expression(ens2, ann2)
    for _ in range(4):
        pd.testing.assert_frame_equal(first, D.deconvolve_states_expression(ens2, ann2))
    # The max-mean row (9999) won, so GFAP is now flat and AC-like loses its edge.
    sym = D.to_symbol_matrix(ens2, ann2)
    assert float(sym.loc["GFAP"].min()) == pytest.approx(9999.0)


# ---------------------------------------------------------------------------
# Library-size / batch-effect behaviour
#
# This cohort has a severe library-prep batch effect: a large subset of samples
# carries the majority of its TPM mass in misc_RNA transcripts. Marker-
# enrichment scoring reads relative abundance, so these two tests pin down
# which distortions it survives and which it does not.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("factor", [0.01, 3.0, 1000.0])
@pytest.mark.parametrize("fraction_method", ["softmax_z", "proportional"])
def test_rescaling_a_sample_does_not_change_its_fractions(factor, fraction_method):
    """SCALE INVARIANCE HOLDS, exactly, for the default ssGSEA path.

    ssGSEA scores depend only on the *within-sample ranking* of genes, and
    multiplying a whole column by a positive constant is a strictly monotone
    transform that leaves every rank untouched. So a pure library-size or
    global-normalisation difference between samples cannot move the output at
    all -- not approximately, bit-for-bit. Both fraction rules inherit this,
    because both consume the score matrix rather than expression.

    What this does NOT buy us is immunity to the misc_RNA batch effect, which
    is compositional rather than a pure rescaling -- see
    `test_compositional_distortion_does_move_fractions` below.
    """
    mat, _ = state_matrix()
    base = D.deconvolve_states_expression(mat, fraction_method=fraction_method)
    scaled = mat.copy()
    scaled.iloc[:, 0] = scaled.iloc[:, 0] * factor
    after = D.deconvolve_states_expression(scaled, fraction_method=fraction_method)
    np.testing.assert_allclose(
        base[STATES].to_numpy(), after[STATES].to_numpy(), rtol=0, atol=0
    )


def test_compositional_distortion_does_move_fractions():
    """SCALE INVARIANCE IS NOT COMPOSITIONAL INVARIANCE -- documented, not fixed.

    The misc_RNA batch effect is not a rescaling. It injects a large block of
    highly-expressed non-marker transcripts, which in TPM space necessarily
    compresses everything else *and* pushes those transcripts above the marker
    genes in the within-sample ranking. ssGSEA ranks against the whole gene
    universe, so marker ranks fall and enrichment scores move.

    This test injects exactly that distortion into half of a synthetic cohort
    and asserts the fractions DO change measurably. That is a real limitation
    of marker-enrichment scoring on this cohort, and the honest handling is to
    pin it here rather than to claim robustness the method does not have.

    Mitigations a later, reviewed step should consider -- none are implemented
    here -- are restricting the ranking universe to protein-coding genes,
    re-deriving TPM over protein-coding genes only, or regressing the
    misc_RNA mass fraction out of the scores. The batch variable is partly
    confounded with tissue source site, so none of them is free.
    """
    mat, assignments = state_matrix(reps=6)
    base = D.deconvolve_states_expression(mat).set_index("patient_id")

    # Inject a misc_RNA-like block into half the samples: 2000 non-marker
    # transcripts expressed well above the marker panels, then renormalise
    # every column to 1e6 as real TPM would be.
    distorted = mat.copy()
    junk = pd.DataFrame(
        0.0, index=[f"MISCRNA{i}" for i in range(2000)], columns=distorted.columns
    )
    distorted = pd.concat([distorted, junk])
    hit = list(distorted.columns)[::2]
    distorted.loc[junk.index, hit] = 5000.0
    distorted = distorted / distorted.sum(axis=0) * 1e6

    after = D.deconvolve_states_expression(distorted).set_index("patient_id")
    delta = (after[STATES] - base.loc[after.index, STATES]).abs().to_numpy().max()
    assert delta > 0.01, (
        "expected the compositional distortion to move the fractions; if this "
        "now fails, the scoring was made compositionally robust and this test "
        "should be rewritten to assert that instead"
    )

    # The distortion degrades but does not destroy the planted signal: the
    # argmax is still right for most samples. This bounds the damage.
    correct = sum(
        after.loc[s, STATES].astype(float).idxmax() == D._STATE_TO_COLUMN[p]
        for s, p in assignments.items()
    )
    assert correct >= 0.5 * len(assignments), f"{correct}/{len(assignments)} argmax survived"


# ---------------------------------------------------------------------------
# Degenerate inputs
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "fn,cols,builder,uniform",
    [
        (D.deconvolve_states_expression, STATES, state_matrix, 0.25),
        (D.deconvolve_tme_expression, TME, tme_matrix, 1 / 3),
    ],
)
def test_all_zero_sample_gets_uniform_fractions(fn, cols, builder, uniform):
    """An all-zero column is entirely ties, so every ssGSEA score is exactly 0.

    The right answer for a sample carrying no information is "no preference",
    not an arbitrary winner. The compositional contract still holds.
    """
    mat, _ = builder()
    mat = mat.copy()
    mat["ZERO"] = 0.0
    out = fn(mat).set_index("patient_id")
    D.check_fractions(out.reset_index(), tuple(cols))
    row = out.loc["ZERO", cols].astype(float).to_numpy()
    np.testing.assert_allclose(row, uniform, atol=1e-9)
    assert row.sum() == pytest.approx(1.0)


def test_constant_nonzero_sample_also_gets_uniform_fractions():
    mat, _ = state_matrix()
    mat = mat.copy()
    mat["FLAT"] = 77.0
    scores = D.state_score_matrix(mat).scores
    np.testing.assert_allclose(scores.loc["FLAT"].to_numpy(), 0.0, atol=1e-12)
    out = D.deconvolve_states_expression(mat).set_index("patient_id")
    np.testing.assert_allclose(out.loc["FLAT", STATES].astype(float).to_numpy(), 0.25, atol=1e-9)


@pytest.mark.parametrize(
    "fn,cols,builder",
    [
        (D.deconvolve_states_expression, STATES, state_matrix),
        (D.deconvolve_tme_expression, TME, tme_matrix),
    ],
)
def test_empty_sample_set_returns_empty_frame_with_contract_columns(fn, cols, builder):
    mat, _ = builder()
    out = fn(mat.iloc[:, :0])
    assert len(out) == 0
    for col in ["patient_id", *cols]:
        assert col in out.columns
    assert D.check_fractions(out, tuple(cols))["n"] == 0
    assert len(D.fraction_distribution(out, tuple(cols))) == 0


def test_single_sample_works_with_ssgsea_and_zmean_is_rejected():
    mat = synthetic_matrix(D.NEFTEL_SIGNATURES, {"ONLY": "MES-like"})
    out = D.deconvolve_states_expression(mat)
    assert len(out) == 1
    assert out[STATES].iloc[0].astype(float).idxmax() == "state_mes"
    assert out[STATES].to_numpy().sum() == pytest.approx(1.0)
    with pytest.raises(ValueError, match="at least 2 samples"):
        D.deconvolve_states_expression(mat, score_method="zmean")


def test_nan_expression_is_rejected():
    mat, _ = state_matrix(reps=2)
    mat = mat.copy()
    mat.iloc[0, 0] = np.nan
    with pytest.raises(ValueError, match="NaN"):
        D.deconvolve_states_expression(mat)


def test_unknown_methods_are_rejected():
    mat, _ = state_matrix(reps=2)
    with pytest.raises(ValueError, match="unknown method"):
        D.deconvolve_states_expression(mat, score_method="cibersortx")
    with pytest.raises(ValueError, match="unknown method"):
        D.deconvolve_states_expression(mat, fraction_method="nnls")


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------
def test_fraction_distribution_summarises_every_column():
    mat, assignments = state_matrix()
    out = D.deconvolve_states_expression(mat)
    dist = D.fraction_distribution(out, D.STATES)
    assert list(dist["column"]) == STATES
    assert dist["n_argmax"].sum() == len(out)
    # Balanced planted design -> each state wins its own four samples.
    assert (dist["n_argmax"] == 4).all()
    assert ((dist["mean"] > 0) & (dist["mean"] < 1)).all()


def test_marker_fraction_concordance_detects_planted_markers():
    mat, _ = state_matrix(reps=6)
    out = D.deconvolve_states_expression(mat)
    conc = D.marker_fraction_concordance(mat, out, D.STATE_CONCORDANCE_MARKERS)
    assert len(conc) == len(D.STATE_CONCORDANCE_MARKERS)
    assert conc["concordant"].all(), conc.to_dict("records")
    assert (conc["mean_in_marker_high"] > conc["mean_in_rest"]).all()
    assert (conc["spearman_r"] > 0).all()


def test_marker_fraction_concordance_handles_tme_and_absent_markers():
    mat, _ = tme_matrix(reps=6)
    out = D.deconvolve_tme_expression(mat)
    conc = D.marker_fraction_concordance(mat, out, D.TME_CONCORDANCE_MARKERS)
    assert conc["concordant"].all(), conc.to_dict("records")

    absent = D.marker_fraction_concordance(mat, out, {"NOT_A_GENE": "tme_immune"})
    assert len(absent) == 1
    assert not bool(absent.loc[0, "concordant"])
    assert np.isnan(absent.loc[0, "spearman_r"])


def test_marker_fraction_concordance_on_empty_input():
    mat, _ = state_matrix(reps=2)
    empty = D.deconvolve_states_expression(mat.iloc[:, :0])
    conc = D.marker_fraction_concordance(mat, empty, D.STATE_CONCORDANCE_MARKERS)
    assert len(conc) == 0


def test_immune_signature_agreement_is_positive_on_planted_data():
    mat, _ = tme_matrix(reps=6, extra_genes=D.BROAD_IMMUNE_PANEL)
    # Plant the broad panel high in exactly the immune-high samples, so the
    # two independent panels should agree.
    imm = [c for c in mat.columns if c.startswith("IMM")]
    mat.loc[list(D.BROAD_IMMUNE_PANEL), imm] = 500.0
    out = D.deconvolve_tme_expression(mat)
    agree = D.immune_signature_agreement(mat, out)
    assert agree["n"] == mat.shape[1]
    assert agree["mapping_rate"] == 1.0
    assert agree["spearman_r"] > 0.5, agree
    assert agree["pearson_r"] > 0.4, agree


def test_immune_signature_agreement_on_empty_input():
    mat, _ = tme_matrix(reps=2, extra_genes=D.BROAD_IMMUNE_PANEL)
    empty = D.deconvolve_tme_expression(mat.iloc[:, :0])
    agree = D.immune_signature_agreement(mat.iloc[:, :0], empty)
    assert agree["n"] == 0


def test_gfap_saturation_caveat_is_documented_in_the_marker_dict():
    """GFAP is a poor discriminator on real GBM; the module must not hide that."""
    assert D.STATE_CONCORDANCE_MARKERS["GFAP"] == "state_ac"
    # Non-saturated astrocytic transporters are offered alongside it.
    assert D.STATE_CONCORDANCE_MARKERS["AQP4"] == "state_ac"
    assert D.STATE_CONCORDANCE_MARKERS["SLC1A2"] == "state_ac"


# ---------------------------------------------------------------------------
# Real-data smoke tests (skipped when the parquet is absent)
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def real_data():
    if not HAVE_REAL:
        pytest.skip("TCGA-GBM expression parquet not present")
    return pd.read_parquet(TPM_PATH), pd.read_parquet(ANN_PATH)


@needs_real
def test_real_tcga_gbm_states_end_to_end(real_data):
    tpm, ann = real_data
    res = D.state_score_matrix(tpm, ann)
    assert res.overall_mapping_rate == 1.0, res.report().to_dict("records")

    out = D.deconvolve_states_expression(tpm, ann)
    assert len(out) == tpm.shape[1]
    D.check_fractions(out, D.STATES)

    dist = D.fraction_distribution(out, D.STATES)
    # Not degenerate: every state wins somewhere, none runs away with the cohort.
    assert (dist["n_argmax"] >= 10).all(), dist.to_dict("records")
    assert dist["n_argmax"].max() < 0.6 * len(out)
    assert (dist["sd"] > 0.01).all()


@needs_real
def test_real_tcga_gbm_tme_end_to_end(real_data):
    tpm, ann = real_data
    res = D.tme_score_matrix(tpm, ann)
    assert res.overall_mapping_rate == 1.0, res.report().to_dict("records")

    out = D.deconvolve_tme_expression(tpm, ann)
    D.check_fractions(out, D.TME_FRACTIONS)
    dist = D.fraction_distribution(out, D.TME_FRACTIONS)
    assert (dist["n_argmax"] >= 10).all(), dist.to_dict("records")


@needs_real
def test_real_data_immune_axis_agrees_with_an_independent_immune_panel(real_data):
    """Convergent evidence: the immune axis tracks a panel it barely shares genes with."""
    tpm, ann = real_data
    out = D.deconvolve_tme_expression(tpm, ann)
    agree = D.immune_signature_agreement(tpm, out, ann)
    assert agree["n"] == tpm.shape[1]
    assert agree["mapping_rate"] == 1.0
    assert agree["spearman_r"] > 0.4, agree
    assert agree["pearson_r"] > 0.3, agree


@needs_real
def test_real_data_marker_concordance_holds_except_for_saturated_gfap(real_data):
    """Measured on n=284: every marker concords EXCEPT GFAP, which is saturated.

    GFAP sits among the top-ranked genes in essentially every GBM sample
    (median ~3290 TPM), so its within-cohort variance carries little
    AC-specific information and the top-quartile check fails for it. The
    non-saturated astrocytic transporters SLC1A2 and AQP4 behave as expected.
    This test pins that exact pattern so a regression in either direction is
    visible.
    """
    tpm, ann = real_data
    out = D.deconvolve_states_expression(tpm, ann)
    conc = D.marker_fraction_concordance(
        tpm, out, D.STATE_CONCORDANCE_MARKERS, ann
    ).set_index("marker")

    for marker in ["AQP4", "SLC1A2", "PDGFRA", "DCX", "CHI3L1"]:
        assert bool(conc.loc[marker, "concordant"]), conc.loc[marker].to_dict()
        assert conc.loc[marker, "spearman_r"] > 0.2, conc.loc[marker].to_dict()

    # GFAP: documented failure, not a silent one.
    assert not bool(conc.loc["GFAP", "concordant"]), (
        "GFAP now concords -- the saturation caveat in STATE_CONCORDANCE_MARKERS "
        "and in the real-data report needs updating"
    )
    assert conc.loc["GFAP", "spearman_r"] < 0.2


@needs_real
def test_real_data_tme_marker_concordance(real_data):
    tpm, ann = real_data
    out = D.deconvolve_tme_expression(tpm, ann)
    conc = D.marker_fraction_concordance(tpm, out, D.TME_CONCORDANCE_MARKERS, ann)
    assert conc["concordant"].all(), conc.to_dict("records")
    assert (conc["spearman_r"] > 0.2).all(), conc.to_dict("records")


@needs_real
def test_real_data_is_scale_invariant(real_data):
    """Library-size differences between real samples cannot move the output."""
    tpm, ann = real_data
    base = D.deconvolve_states_expression(tpm, ann)
    scaled = tpm.copy()
    scaled.iloc[:, 0] = scaled.iloc[:, 0] * 17.0
    scaled.iloc[:, 1] = scaled.iloc[:, 1] * 0.013
    after = D.deconvolve_states_expression(scaled, ann)
    np.testing.assert_allclose(
        base[STATES].to_numpy(), after[STATES].to_numpy(), rtol=0, atol=0
    )


@needs_real
def test_real_data_fractions_shift_with_misc_rna_batch(real_data):
    """The misc_RNA library-prep batch effect DOES shift the real fractions.

    A large subset of this cohort carries the majority of its TPM mass in
    misc_RNA transcripts (median ~50% versus ~0.1% for the rest). Because
    marker-enrichment scoring reads relative abundance and ranks against the
    whole gene universe, the two batches get systematically different state
    fractions. This test measures and pins that, so the effect is on the
    record rather than discovered downstream.

    It is a confounder, not necessarily an artefact of this module alone: the
    batch variable is partly confounded with tissue source site, so a shift
    here is not by itself proof that the biology differs or that it does not.
    Any downstream use of these fractions must adjust for, or stratify on,
    the misc_RNA mass fraction.
    """
    tpm, ann = real_data
    gene_type = ann.set_index("gene_id")["gene_type"].reindex(tpm.index)
    misc_mass = tpm[gene_type.eq("misc_RNA").to_numpy()].sum(axis=0) / tpm.sum(axis=0)
    hi = (misc_mass > 0.3).to_numpy()

    # The batch effect is real and large in this cohort.
    assert 50 < hi.sum() < tpm.shape[1] - 50, f"{hi.sum()} high-misc_RNA samples"
    assert float(misc_mass[hi].median()) > 0.3
    assert float(misc_mass[~hi].median()) < 0.1

    out = D.deconvolve_states_expression(tpm, ann)
    shifts = {c: abs(float(out[c][hi].mean()) - float(out[c][~hi].mean())) for c in STATES}
    assert max(shifts.values()) > 0.02, (
        f"expected a measurable batch shift; got {shifts}. If this now fails, "
        "the scoring was made robust to the misc_RNA batch effect and this "
        "test should be rewritten to assert that instead."
    )
    # Bound the damage: no state is wholly explained by batch.
    assert max(shifts.values()) < 0.25, shifts
