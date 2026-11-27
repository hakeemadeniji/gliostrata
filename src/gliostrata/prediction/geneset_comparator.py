"""T2's comparator: a published 2024 gene set refitted under the project's own conditions.

Pre-registration section 2.2 item 2, as amended 2026-11-27 (option B). The published
*model* is not re-run, because almost no paper states the expression scale its
coefficients expect. Instead, the published *gene set* is refitted the same way for every
candidate:

* **Expression:** log2(CPM + 1), computed for each sample from its own STAR read counts.
  The library size is summed over a fixed gene universe, the symbols present in both
  the TCGA and the CGGA matrices, so no sample's value depends on any other sample's.
* **Pairs:** a gene pair is 1 when its first gene's value exceeds its second's, and 0
  otherwise (a tie is 0).
* **Fit:** the gene set's terms enter :class:`CoxPHModel` (Breslow, ridge 1e-6), the
  estimator of the T2 variants. No clinical covariate is added, and no term is selected
  or dropped.
* **Choice:** 5-fold stratified cross-validation (:func:`cross_validate_cox`, seed 42) on
  the patients of the development training split. The highest pooled out-of-fold
  C-index wins; ties go to the earlier publication, then to the lower PMID.

Nothing here reads CGGA. The evaluation script supplies CGGA's expression only after the
deposit, through the read guard. :func:`paired_bootstrap_difference` is the registered
paired percentile bootstrap of the C-index difference.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from gliostrata.prediction.cox import CoxPHModel, concordance_index, cross_validate_cox

__all__ = [
    "GeneSet",
    "collapse_to_symbols",
    "log2_cpm",
    "design_matrix",
    "cv_cindex",
    "choose_gene_set",
    "fit_gene_set",
    "score_gene_set",
    "paired_bootstrap_difference",
]

CV_SEED = 42
CV_SPLITS = 5


@dataclass(frozen=True)
class GeneSet:
    """One published 2024 gene set, as named in its paper."""

    pmid: str
    first_published: str  # ISO date, for the tie rule
    genes: tuple[str, ...] = ()
    pairs: tuple[tuple[str, str], ...] = ()
    source: str = ""  # where the paper lists the terms
    label: str = ""
    notes: str = field(default="", compare=False)

    def __post_init__(self) -> None:
        if not self.genes and not self.pairs:
            raise ValueError(f"gene set {self.pmid}: no genes and no pairs")
        terms = list(self.genes) + [f"{a}>{b}" for a, b in self.pairs]
        if len(set(terms)) != len(terms):
            raise ValueError(f"gene set {self.pmid}: a term is listed twice")
        for a, b in self.pairs:
            if a == b:
                raise ValueError(f"gene set {self.pmid}: pair {a}|{b} compares a gene with itself")

    @property
    def symbols(self) -> frozenset[str]:
        return frozenset(self.genes) | frozenset(g for p in self.pairs for g in p)

    @property
    def term_names(self) -> list[str]:
        return list(self.genes) + [f"{a}>{b}" for a, b in self.pairs]

    def missing_from(self, universe: Iterable[str]) -> list[str]:
        u = set(universe)
        return sorted(s for s in self.symbols if s not in u)


def collapse_to_symbols(counts: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    """Ensembl-indexed counts (genes x samples) to symbol-indexed counts.

    The counts of every Ensembl id that shares a GENCODE ``gene_name`` are summed, so
    the result is keyed like a symbol-level matrix such as CGGA's. Ids with no
    annotated name are dropped.
    """
    if not {"gene_id", "gene_name"} <= set(annotation.columns):
        raise ValueError("annotation needs gene_id and gene_name")
    names = annotation.set_index("gene_id")["gene_name"]
    missing = counts.index.difference(names.index)
    if len(missing):
        raise ValueError(f"{len(missing)} count rows have no annotation")
    sym = names.reindex(counts.index)
    keep = sym.notna() & (sym.astype(str).str.len() > 0)
    out = counts.loc[keep].groupby(sym[keep].astype(str).to_numpy()).sum()
    out.index.name = "gene_name"
    return out


def log2_cpm(counts_by_symbol: pd.DataFrame, universe: Sequence[str]) -> pd.DataFrame:
    """log2(CPM + 1) per sample, with the library size summed over ``universe``.

    ``counts_by_symbol`` is genes (symbols) x samples. Only the rows in ``universe``
    are used, and every one of them must be present. Each sample's value depends on
    that sample's counts alone.
    """
    u = pd.Index(list(dict.fromkeys(universe)))
    absent = u.difference(counts_by_symbol.index)
    if len(absent):
        raise ValueError(f"{len(absent)} universe genes are not in the matrix, e.g. {list(absent[:5])}")
    sub = counts_by_symbol.loc[u].astype(float)
    if (sub.to_numpy() < 0).any():
        raise ValueError("negative counts")
    lib = sub.sum(axis=0)
    if (lib <= 0).any():
        raise ValueError("a sample has no counts over the universe")
    return np.log2(sub.div(lib, axis=1) * 1e6 + 1.0)


def design_matrix(expr: pd.DataFrame, gene_set: GeneSet, samples: Sequence[str]) -> np.ndarray:
    """Samples x terms: one column per gene, then one 0/1 column per pair."""
    missing = gene_set.missing_from(expr.index)
    if missing:
        raise ValueError(f"gene set {gene_set.pmid}: not in the expression matrix: {missing}")
    sub = expr.loc[:, list(samples)]
    cols = [sub.loc[g].to_numpy(dtype=float) for g in gene_set.genes]
    cols += [(sub.loc[a].to_numpy(dtype=float) > sub.loc[b].to_numpy(dtype=float)).astype(float)
             for a, b in gene_set.pairs]
    return np.column_stack(cols)


def _check_terms_vary(x_mat: np.ndarray, gene_set: GeneSet) -> None:
    const = [n for n, sd in zip(gene_set.term_names, x_mat.std(axis=0)) if not sd > 0]
    if const:
        raise ValueError(f"gene set {gene_set.pmid}: terms constant in the fitting data: {const}")


def cv_cindex(gene_set: GeneSet, expr: pd.DataFrame, samples: Sequence[str],
              time: np.ndarray, event: np.ndarray, *, seed: int = CV_SEED,
              n_splits: int = CV_SPLITS, n_boot: int = 1000):
    """5-fold stratified CV of the refitted gene set; returns the :class:`CVResult`."""
    x_mat = design_matrix(expr, gene_set, samples)
    _check_terms_vary(x_mat, gene_set)
    return cross_validate_cox(x_mat, time, event, n_splits=n_splits, seed=seed, n_boot=n_boot)


def choose_gene_set(gene_sets: Sequence[GeneSet], expr: pd.DataFrame, samples: Sequence[str],
                    time: np.ndarray, event: np.ndarray, *, seed: int = CV_SEED,
                    n_boot: int = 1000) -> tuple[pd.DataFrame, GeneSet]:
    """Score every gene set by CV on the same patients and pick the registered winner.

    Highest pooled out-of-fold C-index; ties go to the earlier publication, then to the
    lower PMID. A gene set that cannot be fitted (a constant term) is reported and
    skipped, never silently repaired.
    """
    if not gene_sets:
        raise ValueError("no eligible gene set")
    rows = []
    for gs in gene_sets:
        try:
            res = cv_cindex(gs, expr, samples, time, event, seed=seed, n_boot=n_boot)
            rows.append({"pmid": gs.pmid, "first_published": gs.first_published,
                         "n_terms": len(gs.term_names), "cv_cindex": res.c_index,
                         "cv_ci_low": res.c_index_ci[0], "cv_ci_high": res.c_index_ci[1],
                         "n_obs": res.n_obs, "n_events": res.n_events, "error": None})
        except ValueError as exc:
            rows.append({"pmid": gs.pmid, "first_published": gs.first_published,
                         "n_terms": len(gs.term_names), "cv_cindex": np.nan,
                         "cv_ci_low": np.nan, "cv_ci_high": np.nan, "n_obs": len(samples),
                         "n_events": int(np.sum(event)), "error": str(exc)})
    table = pd.DataFrame(rows)
    ok = table[table["error"].isna()]
    if ok.empty:
        raise ValueError("no gene set could be fitted")
    ranked = ok.sort_values(["cv_cindex", "first_published", "pmid"],
                            ascending=[False, True, True], kind="mergesort")
    winner = ranked.iloc[0]["pmid"]
    table["rank"] = table["pmid"].map({p: i + 1 for i, p in enumerate(ranked["pmid"])})
    table["chosen"] = table["pmid"] == winner
    return table.sort_values("rank", na_position="last").reset_index(drop=True), \
        next(g for g in gene_sets if g.pmid == winner)


def fit_gene_set(gene_set: GeneSet, expr: pd.DataFrame, samples: Sequence[str],
                 time: np.ndarray, event: np.ndarray) -> CoxPHModel:
    x_mat = design_matrix(expr, gene_set, samples)
    _check_terms_vary(x_mat, gene_set)
    model = CoxPHModel()
    model.fit(x_mat, np.asarray(time, dtype=float), np.asarray(event).astype(int),
              feature_names=gene_set.term_names)
    return model


def score_gene_set(model: CoxPHModel, gene_set: GeneSet, expr: pd.DataFrame,
                   samples: Sequence[str]) -> np.ndarray:
    """Log-hazard risk (higher is worse) for each sample, from its own values only."""
    return model.predict_log_hazard(design_matrix(expr, gene_set, samples))


def paired_bootstrap_difference(time, event, risk_ours, risk_comparator, *, n_boot: int = 1000,
                                seed: int = 42) -> dict:
    """Paired percentile bootstrap of C(ours) - C(comparator) over the same patients.

    Both risk vectors are fixed; only the patients are resampled, with the same indices
    for both models in each resample.
    """
    time = np.asarray(time, dtype=float)
    event = np.asarray(event).astype(int)
    a = np.asarray(risk_ours, dtype=float)
    b = np.asarray(risk_comparator, dtype=float)
    if not (len(time) == len(event) == len(a) == len(b)):
        raise ValueError("time, event and both risks must be the same length")
    c_a = concordance_index(time, event, a)
    c_b = concordance_index(time, event, b)
    rng = np.random.default_rng(seed)
    n = len(time)
    diffs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        try:
            diffs.append(concordance_index(time[idx], event[idx], a[idx])
                         - concordance_index(time[idx], event[idx], b[idx]))
        except ValueError:
            continue
    d = np.asarray(diffs, dtype=float)
    return {"cindex_ours": c_a, "cindex_comparator": c_b, "difference": c_a - c_b,
            "ci_low": float(np.quantile(d, 0.025)), "ci_high": float(np.quantile(d, 0.975)),
            "n": n, "n_events": int(event.sum()), "n_boot_valid": len(diffs),
            "n_boot": n_boot, "seed": seed}
