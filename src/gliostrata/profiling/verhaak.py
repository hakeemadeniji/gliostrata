"""L2 -- Verhaak GBM transcriptional subtype scoring from bulk expression.

WHAT THIS MODULE IS, AND WHAT IT IS NOT
=======================================
`gliostrata.profiling.subtype` is a PASSTHROUGH STUB: `classify_subtypes()`
reads an existing ``subtype`` column and emits noisy probabilities centred on
it. It derives nothing from expression data. This module is built ALONGSIDE
that stub (exactly as `prediction/cox.py` was built alongside
`prediction/survival.py`); integration is a separate, reviewed step.

The GlioStrata protocol specifies a **supervised** random-forest + 1D-CNN
ensemble trained on labelled TCGA-GBM samples, evaluated against CGGA with a
pre-registered accuracy threshold (``t1_subtype_accuracy_cgga = 0.88``).

**That is not implementable here, and this module does not implement it.**
``data/raw/tcga_gbm/clinical_cohort.parquet`` has a ``subtype`` column that is
100% NULL. There are no subtype labels anywhere in this project. A supervised
classifier cannot be trained, and -- just as importantly -- no accuracy can be
computed, because accuracy requires labels. Any number presented as "subtype
accuracy" from this codebase today would be fabricated.

What this module does instead is **unsupervised, signature-based single-sample
subtype scoring**: score each sample against published per-subtype marker
signatures, then turn those scores into a calibrated-looking probability vector
and an argmax label. This is the same family of method used for single-sample
subtype assignment in the literature (ssGSEA-based scoring, Verhaak 2010's own
"ssGSEA of the 840-gene signature" assignment step, and the single-sample
classifier of Wang 2017), so it is scientifically defensible. It is NOT a
trained classifier and carries no validated accuracy.

PROTOCOL TEST T1 IS NOT ADDRESSED BY THIS MODULE.

Biology and literature basis of the signatures
----------------------------------------------
Verhaak et al., *Cancer Cell* 17:98-110 (2010) defined four TCGA-GBM
expression subtypes -- Proneural, Neural, Classical, Mesenchymal -- from an
840-gene signature, each anchored on specific genomic/pathway events. Wang
et al., *Cancer Cell* 32:42-56 (2017) re-analysed tumour-core-enriched IDH
wild-type GBM and showed that the **Neural** subtype is not a tumour-intrinsic
class: it is driven by normal neuronal/oligodendrocytic tissue admixture at the
tumour margin and largely disappears when non-tumour signal is removed. Wang
2017 therefore recognises three subtypes (CL / MES / PN).

This module scores all four, because downstream GlioStrata code (the stub's
``SUBTYPES`` list, the feature store, the CLI) expects four
``subtype_prob_*`` columns -- but the Neural signature here is deliberately
composed of neuronal genes and should be read as a **normal-brain
contamination axis**, not a tumour subtype. ``NEURAL_IS_CONTAMINATION_AXIS``
records this in machine-readable form, and a high Neural score is better
interpreted as "this sample contains a lot of non-tumour brain" than as a
biological class.

Signature construction policy
-----------------------------
Each list below is 14-24 genes of *canonical, individually citable* markers
for that subtype, not a reconstruction of the full published 840-gene
signature (which is not redistributed here and which I will not invent).
Breadth that cannot be sourced is worse than a short defensible list, so every
gene is a named marker or a named pathway member from the papers cited inline.
The four lists are kept mutually disjoint so that a gene cannot vote for two
subtypes at once.

Scoring
-------
Default ``method="ssgsea"``: single-sample rank enrichment in the sense of
Barbie et al., *Nature* 462:108-112 (2009) -- rank genes within each sample,
weight by rank^0.25, and integrate the weighted/unweighted ECDF difference over
the ranked list. It is computed here in closed form (see ``_ssgsea_scores``)
rather than by an explicit random walk, which is exact and much faster. It is
cohort-independent: a sample's score does not change if other samples are
added or removed.

Alternative ``method="zmean"``: mean of gene-wise z-scores of log2(TPM+1)
across the cohort. Cohort-dependent by construction, and closer to how Verhaak
2010 themselves centred data, but it needs at least two samples.

Probabilities are a softmax over the *cohort z-scores* of the per-subtype
scores. This makes them comparable across subtypes (raw enrichment scores for
different signatures are not on a common scale) but it also means they are
**relative within the supplied cohort**, and they are calibrated to nothing.
They are a confidence ordering, not a posterior probability of a label.

Everything in this module is fully deterministic: no random number generator is
used anywhere. ``seed`` is accepted only for signature compatibility with the
stub and has no effect.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import silhouette_score

__all__ = [
    "SUBTYPES",
    "VERHAAK_SIGNATURES",
    "NEURAL_IS_CONTAMINATION_AXIS",
    "CONCORDANCE_MARKERS",
    "SignatureMapping",
    "SubtypeScores",
    "strip_ensembl_version",
    "resolve_symbols",
    "to_symbol_matrix",
    "map_signatures",
    "mapping_report",
    "score_signatures",
    "verhaak_score_matrix",
    "classify_verhaak",
    "subtype_distribution",
    "marker_concordance",
    "score_silhouette",
    "PUBLISHED_TCGA_GBM_PROPORTIONS",
    "SCORING_VERSION",
]

SCORING_VERSION = "0.1.0-signature-ssgsea"

#: Subtype order, matching ``gliostrata.profiling.subtype.SUBTYPES`` exactly so
#: the emitted ``subtype_prob_*`` column set is drop-in compatible.
SUBTYPES: tuple[str, ...] = ("Classical", "Mesenchymal", "Proneural", "Neural")

#: Wang et al. 2017 (Cancer Cell 32:42) showed the Neural class is non-tumour
#: (neuronal) admixture rather than a tumour-intrinsic subtype. Treat a high
#: Neural score as a contamination / tumour-margin indicator.
NEURAL_IS_CONTAMINATION_AXIS = True

# ---------------------------------------------------------------------------
# Signature gene lists. Provenance is given per gene group.
#
# Every symbol is an HGNC symbol as used in the GENCODE v36 annotation shipped
# with this project (gene_annotation.parquet). Where a symbol was renamed after
# the original papers, the GENCODE v36 spelling is used and the legacy name is
# noted (e.g. SEPT11 -> SEPTIN11, IL8 -> CXCL8).
# ---------------------------------------------------------------------------
VERHAAK_SIGNATURES: dict[str, tuple[str, ...]] = {
    # Classical. Verhaak 2010: defined by chromosome 7 gain / chromosome 10
    # loss with high-level EGFR amplification and EGFR overexpression, and by
    # Notch (NOTCH3, JAG1, LFNG) and Sonic-hedgehog (SMO, GAS1, GLI2) pathway
    # activation; the neural-precursor marker NES is highly expressed; AKT2 and
    # FGFR3 are named CL markers. SEC61G is immediately adjacent to EGFR at
    # 7p11.2 and is co-amplified, so it tracks the defining CL lesion; it also
    # appears in the Wang 2017 CL list. CDH4, MEOX2 and SEPTIN11 (SEPT11) are
    # named CL-enriched genes in Verhaak 2010. PTCH1 is included as the
    # canonical SHH receptor completing that pathway module.
    "Classical": (
        "EGFR", "SEC61G", "NES", "AKT2", "FGFR3",
        "NOTCH3", "JAG1", "LFNG", "HEY1",
        "SMO", "GAS1", "GLI2", "PTCH1",
        "CDH4", "MEOX2", "SEPTIN11",
        "PTPRA", "ELOVL2", "ACSBG1", "SOX9",
    ),
    # Mesenchymal. Verhaak 2010: NF1 hemizygous deletion/mutation, high
    # CHI3L1 (YKL-40) and MET, and expression of mesenchymal and astrocytic
    # markers CD44 and MERTK. The subtype is explicitly characterised by
    # activation of the TNF-super-family / NF-kB pathway -- Verhaak 2010 names
    # TRADD, RELB, TNFRSF1A and the TLR/caspase module (TLR2, TLR4, CASP1,
    # CASP4, CASP8) -- and by an angiogenic/invasive programme (SERPINE1,
    # TIMP1, PLAU, PLAUR, LOX, IL6, CXCL8 [= IL8], TGFBI, VDR, BCL3).
    # Collagens are deliberately EXCLUDED here and reserved for the stromal
    # module in `deconvolve.py`, so that stromal content does not masquerade as
    # tumour-intrinsic MES.
    "Mesenchymal": (
        "CHI3L1", "CD44", "MET", "MERTK", "NF1",
        "TRADD", "RELB", "TNFRSF1A", "BCL3",
        "TLR2", "TLR4", "CASP1", "CASP4", "CASP8",
        "SERPINE1", "TIMP1", "PLAU", "PLAUR", "LOX",
        "IL6", "CXCL8", "TGFBI", "VDR",
    ),
    # Proneural. Verhaak 2010: PDGFRA amplification/overexpression and IDH1
    # point mutation are the defining events, with high expression of the
    # oligodendrocytic-development programme -- OLIG2, DLL3, NKX2-2, SOX2 are
    # named explicitly -- and of proneural transcription factors ASCL1 and the
    # SOX family (SOX4, SOX8, SOX11), HES6, MYT1, and the neuronal/receptor
    # genes ERBB3, EPHB1, GRIA2, NCAM1, CNTN1. BCAN and CSPG5 are
    # oligodendrocyte-precursor proteoglycans in the same programme.
    #
    # NOTE: IDH1 *mutation* is the PN-associated lesion. Mutation status is not
    # recoverable from a TPM matrix (IDH1 mRNA level is not informative about
    # R132H), so IDH1 is intentionally NOT scored here. Calling PN from
    # expression alone cannot separate IDH-mutant from IDH-wild-type PN.
    "Proneural": (
        "PDGFRA", "OLIG1", "OLIG2", "DLL1", "DLL3", "NKX2-2",
        "SOX2", "SOX4", "SOX8", "SOX11", "ASCL1", "HES6", "MYT1",
        "ERBB3", "EPHB1", "GRIA2", "NCAM1", "CNTN1",
        "BCAN", "CSPG5", "TCF12", "DCX",
    ),
    # Neural -- SCORED AS A NON-TUMOUR CONTAMINATION AXIS, NOT A SUBTYPE.
    # Verhaak 2010 characterised this class by neuron markers NEFL, GABRA1,
    # SYT1 and SLC12A5. Wang 2017 showed the class is produced by normal
    # neuronal tissue admixture. The list below is therefore a deliberately
    # pure mature-neuron panel (pan-neuronal cytoskeleton, synaptic vesicle
    # machinery, ionotropic/GABAergic receptors), which is the right thing to
    # score if the quantity of interest is "how much normal brain is in this
    # biopsy".
    "Neural": (
        "NEFL", "NEFM", "GABRA1", "GABRB2", "SYT1", "SLC12A5",
        "SNAP25", "SYN1", "SYP", "CPLX1", "NRGN",
        "GRIN1", "GRIN2A", "GRIA1", "SLC17A7",
        "RBFOX3", "MAP2", "STMN2", "CAMK2A", "KIF5A",
    ),
}

#: Single-gene markers used for the internal-consistency check. Each is the
#: single most canonical marker of its subtype from Verhaak 2010, and each is
#: expected to be markedly higher in samples called to that subtype. This is an
#: internal-consistency check, not a validation against labels -- the marker is
#: itself part of the signature that produced the call.
CONCORDANCE_MARKERS: dict[str, str] = {
    "EGFR": "Classical",
    "CHI3L1": "Mesenchymal",
    "PDGFRA": "Proneural",
    "NEFL": "Neural",
}

#: Approximate published TCGA-GBM subtype proportions, for a sanity check that
#: calls are not degenerate. Verhaak 2010 (n=200 core TCGA samples) reported
#: roughly even quarters across the four classes; modern re-analyses that strip
#: non-tumour signal (Wang 2017) shrink Neural substantially and distribute
#: that mass across CL/MES/PN. The ranges below are intentionally loose.
PUBLISHED_TCGA_GBM_PROPORTIONS: dict[str, tuple[float, float]] = {
    "Classical": (0.22, 0.38),
    "Mesenchymal": (0.22, 0.40),
    "Proneural": (0.20, 0.38),
    "Neural": (0.02, 0.28),
}

_SSGSEA_ALPHA = 0.25
_ENSEMBL_VERSION_RE = re.compile(r"\.\d+$")


# ---------------------------------------------------------------------------
# Gene identifier handling
#
# These helpers are generic (not Verhaak-specific) and are imported by
# `gliostrata.profiling.deconvolve`, which is part of the same workstream.
# ---------------------------------------------------------------------------
def strip_ensembl_version(gene_ids) -> list[str]:
    """Strip the trailing version suffix from Ensembl gene ids.

    ``ENSG00000146648.18 -> ENSG00000146648``. Ids without a version suffix,
    and non-Ensembl identifiers, are returned unchanged. Non-string entries
    (e.g. NaN) are coerced to ``""``.
    """
    out: list[str] = []
    for gid in gene_ids:
        if not isinstance(gid, str):
            out.append("")
            continue
        out.append(_ENSEMBL_VERSION_RE.sub("", gid))
    return out


def resolve_symbols(gene_ids, annotation: pd.DataFrame) -> list[str]:
    """Map expression-matrix gene ids to HGNC symbols via the annotation table.

    Matching is attempted first on the exact id (``ENSG...18``) and then on the
    version-stripped id, so an annotation table built against a different
    GENCODE patch release still resolves. Ids that resolve to nothing are
    returned as ``""`` -- never silently dropped and never guessed.

    Parameters
    ----------
    gene_ids:
        Iterable of identifiers from the expression matrix index.
    annotation:
        DataFrame with at least ``gene_id`` and ``gene_name`` columns.
    """
    required = {"gene_id", "gene_name"}
    missing = required - set(annotation.columns)
    if missing:
        raise ValueError(
            f"annotation is missing required column(s): {sorted(missing)}. "
            "Expected the gene_annotation.parquet schema (gene_id, gene_name, gene_type)."
        )
    ann_ids = annotation["gene_id"].astype(str).tolist()
    ann_syms = annotation["gene_name"].astype(str).tolist()

    exact: dict[str, str] = dict(zip(ann_ids, ann_syms))
    stripped: dict[str, str] = {}
    for sid, sym in zip(strip_ensembl_version(ann_ids), ann_syms):
        # First annotation row wins, so the result does not depend on hash order.
        stripped.setdefault(sid, sym)

    ids = [g if isinstance(g, str) else "" for g in gene_ids]
    out: list[str] = []
    for gid, sid in zip(ids, strip_ensembl_version(ids)):
        sym = exact.get(gid)
        if sym is None:
            sym = stripped.get(sid)
        out.append("" if sym is None or sym == "nan" else sym)
    return out


def to_symbol_matrix(
    expression: pd.DataFrame,
    annotation: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Return a symbol-indexed copy of a gene x sample expression matrix.

    If ``annotation`` is ``None`` the matrix index is assumed to already hold
    HGNC symbols (useful for synthetic test matrices).

    Duplicate symbols -- 1233 of them in GENCODE v36 -- are collapsed by the
    standard *max-mean probe* rule: of the rows sharing a symbol, the one with
    the highest mean expression across samples is kept. Ties are broken by
    original row order, so the result is deterministic. Unmapped ids are
    dropped (and are reported by :func:`map_signatures` as missing genes, which
    is where the user sees the consequence).
    """
    if expression.empty and expression.shape[0] == 0:
        return expression.copy()

    df = expression.copy()
    if annotation is not None:
        df.index = pd.Index(resolve_symbols(df.index, annotation), name="gene_name")
        df = df.loc[df.index != ""]
    else:
        df.index = pd.Index([str(i) for i in df.index], name="gene_name")

    if df.index.has_duplicates:
        if df.shape[1] == 0:
            df = df.loc[~df.index.duplicated(keep="first")]
        else:
            means = df.to_numpy(dtype=float, copy=False).mean(axis=1)
            order = np.argsort(-means, kind="stable")
            df = df.iloc[order]
            df = df.loc[~df.index.duplicated(keep="first")]
    return df.sort_index()


# ---------------------------------------------------------------------------
# Signature mapping report
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SignatureMapping:
    """Which genes of one signature were actually found in the matrix."""

    name: str
    requested: tuple[str, ...]
    found: tuple[str, ...]
    missing: tuple[str, ...]

    @property
    def n_requested(self) -> int:
        return len(self.requested)

    @property
    def n_found(self) -> int:
        return len(self.found)

    @property
    def mapping_rate(self) -> float:
        if not self.requested:
            return 0.0
        return len(self.found) / len(self.requested)

    def as_row(self) -> dict:
        return {
            "signature": self.name,
            "n_requested": self.n_requested,
            "n_found": self.n_found,
            "mapping_rate": self.mapping_rate,
            "missing_genes": ",".join(self.missing),
        }


def map_signatures(
    symbol_matrix: pd.DataFrame,
    signatures: dict[str, tuple[str, ...]],
) -> dict[str, SignatureMapping]:
    """Resolve each signature against the matrix index, reporting what is missing.

    Genes absent from the matrix are recorded explicitly rather than silently
    ignored. A signature with zero found genes is allowed here (so the caller
    can see the full report) but :func:`score_signatures` refuses to score it.
    """
    present = set(symbol_matrix.index)
    out: dict[str, SignatureMapping] = {}
    for name, genes in signatures.items():
        found = tuple(g for g in genes if g in present)
        missing = tuple(g for g in genes if g not in present)
        out[name] = SignatureMapping(name=name, requested=tuple(genes), found=found, missing=missing)
    return out


def mapping_report(mappings: dict[str, SignatureMapping]) -> pd.DataFrame:
    """Tidy DataFrame of signature mapping rates, one row per signature."""
    if not mappings:
        return pd.DataFrame(
            columns=["signature", "n_requested", "n_found", "mapping_rate", "missing_genes"]
        )
    return pd.DataFrame([m.as_row() for m in mappings.values()])


# ---------------------------------------------------------------------------
# Single-sample scoring
# ---------------------------------------------------------------------------
def _ssgsea_scores(
    values: np.ndarray,
    member_positions: dict[str, np.ndarray],
    *,
    alpha: float = _SSGSEA_ALPHA,
) -> dict[str, np.ndarray]:
    """Closed-form ssGSEA enrichment scores.

    Parameters
    ----------
    values:
        ``genes x samples`` expression array (the ranking universe).
    member_positions:
        signature name -> integer row indices of its member genes in ``values``.

    Returns
    -------
    dict mapping signature name to a length-``n_samples`` score array.

    Method
    ------
    For one sample, rank all ``N`` genes in ascending expression order so the
    most-expressed gene has rank ``N`` (ties get the average rank). Walking the
    list from the highest-expressed gene downward, ssGSEA accumulates

    * ``P_w(i)``  -- the rank-weighted ECDF of the signature, with weight
      ``rank^alpha`` (Barbie et al. 2009, alpha = 0.25), and
    * ``P_ns(i)`` -- the uniform ECDF of the complement,

    and reports the integral of ``P_w - P_ns`` over the list. Because both
    ECDFs are step functions with steps only at signature / non-signature
    positions, that integral has the closed form

        ES = sum_j w_j * a_j  -  [ N(N+1)/2 - sum_j a_j ] / (N - k)

    where ``a_j`` is the ascending rank of signature gene ``j``, ``w_j`` its
    normalised weight ``a_j^alpha / sum a^alpha``, and ``k`` the signature size.
    This is exact -- not an approximation of the random walk -- and is O(k) per
    sample after the one-off ranking, which matters at 60660 genes x 284
    samples x many signatures.

    ES is divided by ``N`` so scores are comparable between matrices of
    different gene-universe size. A sample whose expression is entirely
    constant (e.g. an all-zero sample) has every gene tied, so ``a_j`` is the
    same for all genes and ES is exactly 0 for every signature -- which is the
    desired "no information" answer rather than an arbitrary one.
    """
    n_genes, n_samples = values.shape
    out = {name: np.zeros(n_samples, dtype=float) for name in member_positions}
    if n_genes == 0 or n_samples == 0:
        return out

    total_rank_sum = n_genes * (n_genes + 1) / 2.0
    for s in range(n_samples):
        asc_rank = rankdata(values[:, s], method="average")
        for name, idx in member_positions.items():
            k = idx.size
            if k == 0 or k >= n_genes:
                continue
            a = asc_rank[idx]
            w = a ** alpha
            w_sum = w.sum()
            if w_sum <= 0:
                continue
            w = w / w_sum
            hit_term = float(np.dot(w, a))
            miss_term = (total_rank_sum - float(a.sum())) / (n_genes - k)
            out[name][s] = (hit_term - miss_term) / n_genes
    return out


def _zmean_scores(
    values: np.ndarray,
    member_positions: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Mean gene-wise cohort z-score of log2(x+1) over each signature.

    Cohort-relative by construction: it answers "is this signature high in this
    sample *relative to the other samples supplied*". Genes with zero variance
    across the cohort contribute 0 rather than NaN. Requires >= 2 samples.
    """
    n_genes, n_samples = values.shape
    out = {name: np.zeros(n_samples, dtype=float) for name in member_positions}
    if n_genes == 0 or n_samples == 0:
        return out
    if n_samples < 2:
        raise ValueError(
            "method='zmean' is cohort-relative and needs at least 2 samples; "
            "use method='ssgsea' for single-sample scoring."
        )
    logged = np.log2(values + 1.0)
    mu = logged.mean(axis=1, keepdims=True)
    sd = logged.std(axis=1, ddof=0, keepdims=True)
    z = np.divide(logged - mu, sd, out=np.zeros_like(logged), where=sd > 0)
    for name, idx in member_positions.items():
        if idx.size == 0:
            continue
        out[name] = z[idx, :].mean(axis=0)
    return out


def score_signatures(
    expression: pd.DataFrame,
    signatures: dict[str, tuple[str, ...]],
    *,
    annotation: pd.DataFrame | None = None,
    method: str = "ssgsea",
    min_genes: int = 3,
) -> tuple[pd.DataFrame, dict[str, SignatureMapping]]:
    """Score every sample against every signature.

    Parameters
    ----------
    expression:
        ``gene x sample`` matrix. TPM (or any non-negative normalised
        abundance). Must not contain NaN.
    signatures:
        name -> tuple of HGNC symbols.
    annotation:
        gene_annotation table; ``None`` means ``expression`` is already
        symbol-indexed.
    method:
        ``"ssgsea"`` (default, per-sample, cohort-independent) or ``"zmean"``
        (cohort-relative mean z-score; needs >= 2 samples).
    min_genes:
        A signature with fewer than this many genes present in the matrix is an
        error, not a silently degraded score.

    Returns
    -------
    (scores, mappings) where ``scores`` is ``samples x signatures`` with the
    signature order of the input dict, and ``mappings`` reports per-signature
    gene recovery.
    """
    if method not in {"ssgsea", "zmean"}:
        raise ValueError(f"unknown method {method!r}; expected 'ssgsea' or 'zmean'")

    sym = to_symbol_matrix(expression, annotation)
    mappings = map_signatures(sym, signatures)

    starved = {
        name: m.n_found for name, m in mappings.items() if m.n_found < min_genes
    }
    if starved and sym.shape[0] > 0:
        raise ValueError(
            "signature(s) have too few genes present in the expression matrix "
            f"(min_genes={min_genes}): {starved}. Check that the matrix is "
            "gene-level and that `annotation` matches it; mapping rates are "
            "available from map_signatures()."
        )

    if sym.shape[1] == 0 or sym.shape[0] == 0:
        return (
            pd.DataFrame(index=pd.Index(sym.columns, name="sample"), columns=list(signatures), dtype=float),
            mappings,
        )

    values = sym.to_numpy(dtype=float, copy=False)
    if not np.isfinite(values).all():
        raise ValueError("expression matrix contains NaN/inf; clean it before scoring.")

    pos_lookup = {g: i for i, g in enumerate(sym.index)}
    member_positions = {
        name: np.array([pos_lookup[g] for g in m.found], dtype=int)
        for name, m in mappings.items()
    }

    if method == "ssgsea":
        raw = _ssgsea_scores(values, member_positions)
    else:
        raw = _zmean_scores(values, member_positions)

    scores = pd.DataFrame(
        {name: raw[name] for name in signatures},
        index=pd.Index(sym.columns, name="sample"),
    )
    return scores, mappings


def _softmax_over_cohort_z(scores: pd.DataFrame, temperature: float) -> np.ndarray:
    """Softmax over per-signature cohort z-scores.

    Raw enrichment scores for different signatures are not on a common scale
    (a 20-gene astrocytic panel and a 23-gene NF-kB panel have different
    baseline enrichment), so the scores are first standardised *per signature
    across samples* and only then compared. Consequences, stated plainly:

    * the output is a within-cohort relative ranking, not an absolute posterior;
    * with a single sample there is nothing to standardise against, so the
      scores are centred within the sample instead and the result is far
      weaker;
    * a signature with zero variance across the cohort contributes a flat 0.
    """
    if temperature <= 0:
        raise ValueError("temperature must be > 0")
    x = scores.to_numpy(dtype=float, copy=True)
    n_samples = x.shape[0]
    if n_samples == 0:
        return x
    if n_samples >= 2:
        mu = x.mean(axis=0, keepdims=True)
        sd = x.std(axis=0, ddof=0, keepdims=True)
        z = np.divide(x - mu, sd, out=np.zeros_like(x), where=sd > 0)
    else:
        row_mu = x.mean(axis=1, keepdims=True)
        row_sd = x.std(axis=1, ddof=0, keepdims=True)
        z = np.divide(x - row_mu, row_sd, out=np.zeros_like(x), where=row_sd > 0)
    z = z / float(temperature)
    z = z - z.max(axis=1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)


@dataclass(frozen=True)
class SubtypeScores:
    """Per-subtype score matrix plus its mapping report, for inspection."""

    scores: pd.DataFrame
    mappings: dict[str, SignatureMapping]
    method: str

    def report(self) -> pd.DataFrame:
        return mapping_report(self.mappings)

    @property
    def overall_mapping_rate(self) -> float:
        req = sum(m.n_requested for m in self.mappings.values())
        fnd = sum(m.n_found for m in self.mappings.values())
        return (fnd / req) if req else 0.0


def verhaak_score_matrix(
    expression: pd.DataFrame,
    annotation: pd.DataFrame | None = None,
    *,
    method: str = "ssgsea",
    signatures: dict[str, tuple[str, ...]] | None = None,
) -> SubtypeScores:
    """Return the raw ``samples x 4`` subtype score matrix for inspection.

    This is the object to look at when asking "why was this sample called
    Classical" -- the probabilities emitted by :func:`classify_verhaak` are a
    monotone transform of these scores and hide their magnitude.
    """
    sigs = VERHAAK_SIGNATURES if signatures is None else signatures
    scores, mappings = score_signatures(expression, sigs, annotation=annotation, method=method)
    return SubtypeScores(scores=scores, mappings=mappings, method=method)


def classify_verhaak(
    expression: pd.DataFrame,
    annotation: pd.DataFrame | None = None,
    *,
    method: str = "ssgsea",
    temperature: float = 1.0,
    seed: int = 42,
) -> pd.DataFrame:
    """Score and label each sample's Verhaak subtype from expression.

    This is signature-based unsupervised scoring, NOT the supervised RF + 1D-CNN
    ensemble the protocol specifies -- no subtype labels exist in this project,
    so no supervised model can be trained and no accuracy can be measured. The
    returned probabilities are a within-cohort confidence ordering calibrated
    against nothing. See the module docstring.

    Parameters
    ----------
    expression:
        ``gene x sample`` TPM matrix (Ensembl-id indexed if ``annotation`` is
        given, HGNC-symbol indexed otherwise).
    annotation:
        gene_annotation table with ``gene_id`` / ``gene_name``.
    method:
        ``"ssgsea"`` or ``"zmean"``; see :func:`score_signatures`.
    temperature:
        Softmax temperature on the z-scored signature scores. Larger values
        flatten the probability vector. It does not change the argmax label.
    seed:
        Ignored. Accepted only for signature parity with
        ``gliostrata.profiling.subtype.classify_subtypes``; this function is
        fully deterministic and uses no RNG.

    Returns
    -------
    DataFrame with ``patient_id``, ``subtype``, ``subtype_prob_<Subtype>`` for
    each of the four subtypes (summing to 1 per row), the raw
    ``subtype_score_<Subtype>`` columns, and method/version metadata. Column
    names match what ``gliostrata.profiling.subtype`` emits so this is drop-in
    for downstream consumers.
    """
    del seed  # deterministic by construction; see docstring

    sc = verhaak_score_matrix(expression, annotation, method=method)
    scores = sc.scores
    probs = _softmax_over_cohort_z(scores, temperature)

    out = pd.DataFrame({"patient_id": pd.Index(scores.index).astype(str)})
    if len(out):
        labels = [SUBTYPES[i] for i in np.argmax(probs, axis=1)]
    else:
        labels = []
    out["subtype"] = pd.Series(labels, dtype=object)
    for j, st in enumerate(SUBTYPES):
        out[f"subtype_prob_{st}"] = probs[:, j] if len(out) else pd.Series(dtype=float)
    for j, st in enumerate(SUBTYPES):
        out[f"subtype_score_{st}"] = scores.iloc[:, j].to_numpy() if len(out) else pd.Series(dtype=float)
    out["subtype_call_confidence"] = probs.max(axis=1) if len(out) else pd.Series(dtype=float)
    out["subtype_scoring_method"] = f"signature_{method}"
    out["subtype_classifier_version"] = SCORING_VERSION
    out["subtype_is_supervised"] = False
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Label-free validation
# ---------------------------------------------------------------------------
def subtype_distribution(calls: pd.DataFrame) -> pd.DataFrame:
    """Called-subtype distribution next to published TCGA-GBM proportions.

    The only thing this can show is that the calls are **not degenerate** (all
    one class, or an empty class). Agreement with published proportions is weak
    evidence -- the published proportions come from a different labelling
    procedure on an overlapping cohort -- but strong disagreement would be a
    red flag worth acting on.
    """
    n = len(calls)
    counts = calls["subtype"].value_counts() if n else pd.Series(dtype=int)
    rows = []
    for st in SUBTYPES:
        c = int(counts.get(st, 0))
        lo, hi = PUBLISHED_TCGA_GBM_PROPORTIONS[st]
        frac = (c / n) if n else 0.0
        rows.append(
            {
                "subtype": st,
                "n": c,
                "fraction": frac,
                "published_low": lo,
                "published_high": hi,
                "within_published_range": bool(lo <= frac <= hi) if n else False,
            }
        )
    return pd.DataFrame(rows)


def marker_concordance(
    expression: pd.DataFrame,
    calls: pd.DataFrame,
    annotation: pd.DataFrame | None = None,
    *,
    markers: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Check canonical markers are highest in the samples called to their subtype.

    For each marker gene, report its mean and median expression within each
    called subtype, the fold change of the expected group over the rest, and
    whether the expected group is in fact the highest-expressing group.

    This is an **internal-consistency** check: the marker contributed to the
    score that produced the call, so passing it does not prove the calls are
    biologically correct. Failing it would prove something is wrong. Absent
    labels it is the strongest check available, and it is not circular in the
    strong sense -- each marker is 1 of ~20 genes in its signature, so a marker
    that behaved randomly would not track the call.
    """
    mk = CONCORDANCE_MARKERS if markers is None else markers
    sym = to_symbol_matrix(expression, annotation)
    if len(calls) == 0:
        return pd.DataFrame(
            columns=[
                "marker", "expected_subtype", "n_expected", "mean_in_expected",
                "median_in_expected", "mean_in_others", "fold_change",
                "argmax_subtype", "concordant",
            ]
        )

    pid = calls["patient_id"].astype(str).tolist()
    labels = pd.Series(calls["subtype"].to_numpy(), index=pid)
    cols = [c for c in sym.columns.astype(str) if c in set(pid)]
    sym = sym.loc[:, [c in set(cols) for c in sym.columns.astype(str)]]
    sym.columns = sym.columns.astype(str)

    rows = []
    for marker, expected in mk.items():
        if marker not in sym.index:
            rows.append(
                {
                    "marker": marker, "expected_subtype": expected, "n_expected": 0,
                    "mean_in_expected": np.nan, "median_in_expected": np.nan,
                    "mean_in_others": np.nan, "fold_change": np.nan,
                    "argmax_subtype": None, "concordant": False,
                }
            )
            continue
        vals = sym.loc[marker, cols].astype(float)
        grp = labels.loc[cols]
        by = vals.groupby(grp.to_numpy()).mean()
        in_exp = vals[grp.to_numpy() == expected]
        others = vals[grp.to_numpy() != expected]
        rows.append(
            {
                "marker": marker,
                "expected_subtype": expected,
                "n_expected": int(len(in_exp)),
                "mean_in_expected": float(in_exp.mean()) if len(in_exp) else np.nan,
                "median_in_expected": float(in_exp.median()) if len(in_exp) else np.nan,
                "mean_in_others": float(others.mean()) if len(others) else np.nan,
                "fold_change": (
                    float(in_exp.mean() / others.mean())
                    if len(in_exp) and len(others) and others.mean() > 0
                    else np.nan
                ),
                "argmax_subtype": str(by.idxmax()) if len(by) else None,
                "concordant": bool(len(by) and by.idxmax() == expected),
            }
        )
    return pd.DataFrame(rows)


def score_silhouette(scores: pd.DataFrame, labels) -> float:
    """Silhouette coefficient of the calls in the per-subtype score space.

    Computed on per-signature z-scored columns of the score matrix (the same
    space the argmax is taken in) with Euclidean distance. Interpretation:
    softmax-argmax partitions this 4-d space with linear-ish boundaries, so a
    positive silhouette confirms the partition is geometrically coherent; it is
    NOT evidence that the clusters correspond to real biological subtypes.
    Returns NaN when fewer than two classes are present.
    """
    x = scores.to_numpy(dtype=float, copy=True)
    lab = np.asarray(labels)
    if x.shape[0] < 3 or len(set(lab.tolist())) < 2:
        return float("nan")
    mu = x.mean(axis=0, keepdims=True)
    sd = x.std(axis=0, ddof=0, keepdims=True)
    z = np.divide(x - mu, sd, out=np.zeros_like(x), where=sd > 0)
    return float(silhouette_score(z, lab, metric="euclidean"))
