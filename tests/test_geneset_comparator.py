"""T2 comparator re-implementation (option B). Synthetic data only; no cohort is read."""

import numpy as np
import pandas as pd
import pytest

from gliostrata.prediction.geneset_comparator import (
    GeneSet,
    choose_gene_set,
    collapse_to_symbols,
    design_matrix,
    fit_gene_set,
    log2_cpm,
    paired_bootstrap_difference,
    score_gene_set,
    score_published,
    stratified_concordance_index,
)


def _counts(n_genes=30, n_samples=120, seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(n_genes)]
    samples = [f"S{j}" for j in range(n_samples)]
    lam = rng.gamma(2.0, 50.0, size=(n_genes, 1)) * rng.lognormal(0, 0.3, size=(1, n_samples))
    return pd.DataFrame(rng.poisson(lam), index=genes, columns=samples)


def _survival(signal, seed=1):
    rng = np.random.default_rng(seed)
    t = rng.exponential(np.exp(-signal)) * 20 + 0.1
    c = rng.uniform(5, 60, size=len(signal))
    return np.minimum(t, c), (t <= c).astype(int)


def test_gene_set_rejects_empty_duplicate_and_self_pairs():
    with pytest.raises(ValueError):
        GeneSet(pmid="1", first_published="2024-01-01")
    with pytest.raises(ValueError):
        GeneSet(pmid="1", first_published="2024-01-01", genes=("A", "A"))
    with pytest.raises(ValueError):
        GeneSet(pmid="1", first_published="2024-01-01", pairs=(("A", "A"),))
    gs = GeneSet(pmid="1", first_published="2024-01-01", genes=("A",), pairs=(("B", "C"),))
    assert gs.symbols == {"A", "B", "C"}
    assert gs.term_names == ["A", "B>C"]
    assert gs.missing_from(["A", "B"]) == ["C"]


def test_collapse_sums_ensembl_ids_sharing_a_symbol():
    counts = pd.DataFrame({"s1": [1, 2, 5], "s2": [3, 4, 7]}, index=["E1.1", "E2.1", "E3.1"])
    ann = pd.DataFrame({"gene_id": ["E1.1", "E2.1", "E3.1"], "gene_name": ["A", "A", "B"]})
    out = collapse_to_symbols(counts, ann)
    assert out.loc["A"].tolist() == [3, 7] and out.loc["B"].tolist() == [5, 7]


def test_log2_cpm_is_per_sample_and_uses_only_the_universe():
    counts = _counts()
    universe = list(counts.index[:20])
    full = log2_cpm(counts, universe)
    alone = log2_cpm(counts[["S0"]], universe)
    pd.testing.assert_series_equal(full["S0"], alone["S0"])  # no cross-sample quantity
    lib = counts.loc[universe, "S0"].sum()
    assert full.loc["G0", "S0"] == pytest.approx(np.log2(counts.loc["G0", "S0"] / lib * 1e6 + 1))
    assert list(full.index) == universe  # genes outside the universe are not used
    with pytest.raises(ValueError):
        log2_cpm(counts, universe + ["NOT_A_GENE"])


def test_design_matrix_orders_genes_then_pairs_and_ties_score_zero():
    expr = pd.DataFrame({"s1": [1.0, 2.0, 2.0], "s2": [3.0, 1.0, 1.0]}, index=["A", "B", "C"])
    gs = GeneSet(pmid="9", first_published="2024-05-01", genes=("A",), pairs=(("A", "B"), ("B", "C")))
    x_mat = design_matrix(expr, gs, ["s1", "s2"])
    assert x_mat.tolist() == [[1.0, 0.0, 0.0], [3.0, 1.0, 0.0]]
    with pytest.raises(ValueError, match="not in the expression matrix"):
        design_matrix(expr, GeneSet(pmid="9", first_published="2024-05-01", genes=("Z",)), ["s1"])


def test_choice_prefers_the_informative_gene_set_and_breaks_ties_by_date():
    counts = _counts(n_samples=200)
    expr = log2_cpm(counts, list(counts.index))
    samples = list(expr.columns)
    signal = (expr.loc["G3"] - expr.loc["G3"].mean()) * 2.0
    time, event = _survival(signal.to_numpy())
    good = GeneSet(pmid="200", first_published="2024-06-01", genes=("G3", "G4"))
    noise = GeneSet(pmid="100", first_published="2024-01-01", genes=("G10", "G11"))
    table, chosen = choose_gene_set([noise, good], expr, samples, time, event, n_boot=50)
    assert chosen.pmid == "200"
    assert table.loc[table["chosen"], "pmid"].item() == "200"
    # identical sets: the earlier publication wins
    twin = GeneSet(pmid="300", first_published="2024-02-01", genes=("G3", "G4"))
    _, chosen2 = choose_gene_set([good, twin], expr, samples, time, event, n_boot=50)
    assert chosen2.pmid == "300"


def test_a_constant_term_is_reported_not_repaired():
    counts = _counts(n_samples=80)
    expr = log2_cpm(counts, list(counts.index))
    samples = list(expr.columns)
    time, event = _survival(np.zeros(len(samples)))
    expr.loc["G0"] = 5.0  # constant
    bad = GeneSet(pmid="1", first_published="2024-01-01", genes=("G0",))
    ok = GeneSet(pmid="2", first_published="2024-01-01", genes=("G1",))
    table, chosen = choose_gene_set([bad, ok], expr, samples, time, event, n_boot=20)
    assert chosen.pmid == "2"
    assert "constant" in table.set_index("pmid").loc["1", "error"]


def test_fit_and_score_use_each_sample_alone():
    counts = _counts(n_samples=150)
    expr = log2_cpm(counts, list(counts.index))
    samples = list(expr.columns)
    # make the pair vary across samples
    expr.loc["G6"] = expr.loc["G5"] + np.random.default_rng(9).normal(scale=0.5, size=len(samples))
    time, event = _survival((expr.loc["G2"] - expr.loc["G2"].mean()).to_numpy())
    gs = GeneSet(pmid="5", first_published="2024-01-01", genes=("G2",), pairs=(("G5", "G6"),))
    model = fit_gene_set(gs, expr, samples, time, event)
    r_all = score_gene_set(model, gs, expr, samples)
    r_one = score_gene_set(model, gs, expr, samples[7:8])
    assert r_one[0] == pytest.approx(r_all[7])


def test_paired_bootstrap_difference_is_paired_and_reproducible():
    rng = np.random.default_rng(3)
    risk_a = rng.normal(size=150)
    time, event = _survival(risk_a)
    risk_b = risk_a + rng.normal(scale=2.0, size=150)
    out1 = paired_bootstrap_difference(time, event, risk_a, risk_b, n_boot=200, seed=42)
    out2 = paired_bootstrap_difference(time, event, risk_a, risk_b, n_boot=200, seed=42)
    assert out1 == out2
    assert out1["difference"] == pytest.approx(out1["cindex_ours"] - out1["cindex_comparator"])
    assert out1["ci_low"] <= out1["difference"] <= out1["ci_high"]
    same = paired_bootstrap_difference(time, event, risk_a, risk_a, n_boot=50)
    assert same["difference"] == 0 and same["ci_low"] == 0 and same["ci_high"] == 0


def test_stratified_cindex_counts_only_within_stratum_pairs():
    from gliostrata.prediction.cox import concordance_index

    rng = np.random.default_rng(5)
    risk = rng.normal(size=120)
    time, event = _survival(risk)
    one = np.zeros(120, dtype=int)
    assert stratified_concordance_index(time, event, risk, one) == pytest.approx(
        concordance_index(time, event, risk))
    # a shift between strata changes the pooled C-index but not the stratified one
    strata = np.repeat(["b1", "b2"], 60)
    shifted = risk + np.where(strata == "b2", 10.0, 0.0)
    assert stratified_concordance_index(time, event, shifted, strata) == pytest.approx(
        stratified_concordance_index(time, event, risk, strata))
    assert concordance_index(time, event, shifted) != pytest.approx(
        concordance_index(time, event, risk))
    with pytest.raises(ValueError):
        stratified_concordance_index(time, event, risk, np.arange(120))  # no within pairs


def test_paired_bootstrap_with_strata_reports_both_versions():
    rng = np.random.default_rng(6)
    risk_a = rng.normal(size=100)
    time, event = _survival(risk_a)
    strata = np.repeat(["mRNAseq_693", "mRNAseq_325"], 50)
    risk_b = risk_a + np.where(strata == "mRNAseq_325", 3.0, 0.0) + rng.normal(size=100)
    out = paired_bootstrap_difference(time, event, risk_a, risk_b, strata=strata, n_boot=100)
    assert out["stratified"] is True
    assert out["n_by_stratum"] == {"mRNAseq_325": 50, "mRNAseq_693": 50}
    assert out["cindex_comparator"] == pytest.approx(
        stratified_concordance_index(time, event, risk_b, strata))
    assert "unstratified_cindex_comparator" in out


def test_score_published_is_the_printed_linear_formula():
    expr = pd.DataFrame({"s1": [1.0, 2.0, 0.5], "s2": [0.0, 4.0, 1.0]}, index=["A", "B", "C"])
    coefs = {"A": 0.5, "B": -0.25, "C": 2.0}
    out = score_published(expr, coefs, ["s2", "s1"])
    assert out.tolist() == [0.5 * 0 - 0.25 * 4 + 2 * 1, 0.5 * 1 - 0.25 * 2 + 2 * 0.5]
    with pytest.raises(ValueError, match="not in the matrix"):
        score_published(expr, {"Z": 1.0}, ["s1"])


def test_primary_config_matches_the_printed_formula():
    import json
    from pathlib import Path

    cfg = json.loads(Path("configs/t2_comparator/primary.json").read_text(encoding="utf-8"))
    assert cfg["pmid"] == "38174685"
    assert cfg["coefficients"] == {"CPPED1": 0.223340457, "F3": 0.1517013, "G0S2": 0.046311722,
                                   "MAPK1": -0.219913858, "MME": 0.192197058,
                                   "MMP9": 0.016505054, "MPO": 0.414575274}
    for gene, coef in cfg["coefficients"].items():
        assert f"Exp ({gene})" in cfg["formula_as_printed"]
        assert f"{abs(coef)}" in cfg["formula_as_printed"]
