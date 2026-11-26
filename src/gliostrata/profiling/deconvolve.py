"""L2 -- Neftel cell-state and TME composition scoring from bulk expression.

READ THIS FIRST: THESE ARE NOT CELL FRACTIONS
=============================================
This module produces, for each patient, four numbers that sum to 1 for the
Neftel states and three numbers that sum to 1 for the tumour
microenvironment. **They are not cell-type proportions.** They are
marker-signature enrichment scores that have been normalised to sum to 1.

Normalising scores to sum to 1 makes a vector look like a composition. It does
not make it one. There is no reference profile, no non-negative least squares
or support-vector regression against one, no deconvolution model, and
therefore nothing that estimates *how much* of a sample is a given cell type.
What these columns rank is **relative enrichment of a marker programme, within
the supplied cohort** -- "this tumour is more MES-like than most tumours in
this cohort", not "this tumour is 31% MES-like cells".

Three specific consequences, so nobody has to rediscover them:

1. A cohort of 283 uniformly MES-high tumours would still produce a spread of
   MES fractions, because the default fraction rule is cohort-relative.
2. The fractions are not comparable across cohorts, and not comparable to
   CIBERSORTx output.
3. Multiplying a fraction by a cell count, a volume, or anything else physical
   is invalid.

WHAT THE PROTOCOL ASKS FOR, AND WHY IT IS NOT HERE
--------------------------------------------------
The GlioStrata protocol specifies CIBERSORTx with the Neftel 2019 atlas as the
signature source, plus MCP-counter for stromal/endothelial content.
CIBERSORTx is a **web service requiring registration and manual
upload/download of mixture files**; it has no redistributable local
implementation and cannot be placed inside an automated, reproducible
pipeline. ``states.py``'s ``cibersortx_signature_matrix()`` raises
``NotImplementedError`` for that reason and this module does not change that.
MCP-counter is marker-based and *is* reimplementable locally -- the TME part of
this module is close in spirit to it -- but MCP-counter itself reports
arbitrary-unit abundance scores, explicitly NOT proportions, and the published
MCP-counter marker lists are not reproduced here.

PROTOCOL TEST T1 IS NOT ADDRESSED BY THIS MODULE.
T1 requires state-deconvolution Pearson r >= 0.75 against ground truth
(``t1_state_deconv_r`` in ``configs``). **No ground truth exists in this
project**: there is no matched single-cell dataset, no CIBERSORTx output, and
no simulated mixture with known proportions for the TCGA-GBM samples. This
method has therefore been validated against nothing, its correlation with
truth is unknown, and nothing in this module should be read as evidence for or
against that threshold. Validating T1 would require either a CIBERSORTx run on
these same samples or pseudo-bulk mixtures built from the Neftel atlas with
known mixing weights; both are real options and both are future work.

Relationship to the existing stubs
----------------------------------
``profiling/states.py`` and ``profiling/tme.py`` are passthrough stubs: they
copy pre-existing ``state_*`` / ``tme_*`` columns through and raise on any
cohort that does not already contain them (i.e. on every real cohort). This
module is built ALONGSIDE them, matching their output contract exactly
(``patient_id`` + ``state_npc/state_opc/state_ac/state_mes``; ``patient_id`` +
``tme_immune/tme_stromal/tme_endothelial``), so the lead can integrate it as a
separate reviewed step.

Literature basis
----------------
Neftel et al., *Cell* 178:835-849 (2019) derived four malignant-cell
meta-modules in GBM -- NPC-like, OPC-like, AC-like and MES-like -- from
single-cell RNA-seq of 28 tumours. The marker panels below are the canonical,
individually citable members of each meta-module (Neftel 2019 Fig. 2 /
Table S2 anchor genes), not the full ~50-gene published modules, which are not
redistributed here. TME panels are standard lineage markers (see per-panel
comments).

Scoring machinery (ssGSEA, gene-id resolution, duplicate-symbol collapsing) is
shared with ``profiling/verhaak.py`` and imported from there rather than
duplicated; both files belong to the same workstream.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from gliostrata.profiling.verhaak import (
    SignatureMapping,
    mapping_report,
    score_signatures,
    to_symbol_matrix,
)

__all__ = [
    "STATES",
    "TME_FRACTIONS",
    "NEFTEL_SIGNATURES",
    "TME_SIGNATURES",
    "BROAD_IMMUNE_PANEL",
    "STATE_CONCORDANCE_MARKERS",
    "TME_CONCORDANCE_MARKERS",
    "EnrichmentResult",
    "state_score_matrix",
    "tme_score_matrix",
    "deconvolve_states_expression",
    "deconvolve_tme_expression",
    "scores_to_fractions",
    "check_fractions",
    "fraction_distribution",
    "marker_fraction_concordance",
    "immune_signature_agreement",
    "METHOD_NAME",
    "METHOD_VERSION",
]

METHOD_NAME = "marker_enrichment_normalised"
METHOD_VERSION = "0.1.0-not-deconvolution"

#: Column contract of ``gliostrata.profiling.states``.
STATES: tuple[str, ...] = ("state_npc", "state_opc", "state_ac", "state_mes")
#: Column contract of ``gliostrata.profiling.tme``.
TME_FRACTIONS: tuple[str, ...] = ("tme_immune", "tme_stromal", "tme_endothelial")

_STATE_TO_COLUMN = {
    "NPC-like": "state_npc",
    "OPC-like": "state_opc",
    "AC-like": "state_ac",
    "MES-like": "state_mes",
}
_TME_TO_COLUMN = {
    "immune": "tme_immune",
    "stromal": "tme_stromal",
    "endothelial": "tme_endothelial",
}

# ---------------------------------------------------------------------------
# Neftel 2019 meta-module marker panels.
#
# Panels are kept mutually disjoint so no gene votes for two states. Where a
# gene is canonical for two modules in the original paper (ASCL1 and DLL3 span
# NPC/OPC boundaries; CD44 and VIM are MES but also broadly mesenchymal), it is
# assigned to the module it discriminates best and the choice is noted.
# ---------------------------------------------------------------------------
NEFTEL_SIGNATURES: dict[str, tuple[str, ...]] = {
    # NPC-like (neural-progenitor-like). Neftel 2019 anchors: DLL3, DCX, CD24,
    # SOX4, SOX11, STMN1/STMN2, TUBB3, ELAVL4, NNAT, HES6, ASCL1, CCND2 --
    # a neurodevelopmental progenitor programme (NPC1 and NPC2 sub-modules
    # combined here; the paper's NPC1/NPC2 split is not reproduced).
    # MLLT11 and TAGLN3 are neuronal-differentiation genes in the same module.
    "NPC-like": (
        "DLL3", "DCX", "CD24", "SOX4", "SOX11", "STMN1", "STMN2",
        "TUBB3", "ELAVL4", "NNAT", "HES6", "ASCL1", "CCND2",
        "MLLT11", "TAGLN3",
    ),
    # OPC-like (oligodendrocyte-progenitor-like). Neftel 2019 anchors: PDGFRA,
    # OLIG1, OLIG2, APOD, CSPG4 (NG2), SOX10, BCAN, SCRG1, OMG, LHFPL3, PLLP,
    # ALCAM, PTPRZ1, NEU4 -- the canonical OPC lineage programme.
    "OPC-like": (
        "PDGFRA", "OLIG1", "OLIG2", "APOD", "CSPG4", "SOX10",
        "BCAN", "SCRG1", "OMG", "LHFPL3", "PLLP", "ALCAM",
        "PTPRZ1", "NEU4",
    ),
    # AC-like (astrocyte-like). Neftel 2019 anchors: GFAP, AQP4, S100B, ALDOC,
    # MLC1, SPARCL1, CLU, ATP1A2, SLC1A2 (EAAT2), SLC1A3 (EAAT1/GLAST), GJA1
    # (connexin-43), MT3, AGT, CST3, HEPACAM, GATM -- mature astrocytic
    # structural, water/glutamate transport and gap-junction genes.
    "AC-like": (
        "GFAP", "AQP4", "S100B", "ALDOC", "MLC1", "SPARCL1", "CLU",
        "ATP1A2", "SLC1A2", "SLC1A3", "GJA1", "MT3", "AGT", "CST3",
        "HEPACAM", "GATM",
    ),
    # MES-like (mesenchymal-like). Neftel 2019 anchors: CHI3L1 (YKL-40), CD44,
    # VIM, ANXA1, ANXA2, LGALS1, LGALS3, S100A11, TIMP1, CAV1, EMP1, THBS1,
    # CLIC1, PLAUR -- the hypoxia/glycolysis-associated MES programme (MES1 and
    # MES2 sub-modules combined; the paper's hypoxia-specific MES2 split is not
    # reproduced). Note the overlap in biology with the Verhaak Mesenchymal
    # subtype: that is real, not an artefact.
    "MES-like": (
        "CHI3L1", "CD44", "VIM", "ANXA1", "ANXA2", "LGALS1", "LGALS3",
        "S100A11", "TIMP1", "CAV1", "EMP1", "THBS1", "CLIC1", "PLAUR",
    ),
}

# ---------------------------------------------------------------------------
# TME lineage marker panels. Standard lineage markers; the same genes used by
# marker-based abundance estimators (MCP-counter, xCell, ESTIMATE immune/
# stromal scores) without reproducing any of their published panels verbatim.
# ---------------------------------------------------------------------------
TME_SIGNATURES: dict[str, tuple[str, ...]] = {
    # Immune. PTPRC (CD45) is the pan-leukocyte marker. In GBM the immune
    # infiltrate is dominated by microglia and monocyte-derived macrophages,
    # so the panel is weighted to myeloid identity -- CD68, AIF1 (Iba1),
    # ITGAM (CD11b), CSF1R, TYROBP, FCER1G, C1QA/C1QB, CD14, LAPTM5, CD53,
    # HLA-DRA -- with CD3E/CD2 for the (sparse) T-cell compartment. A high
    # value here means "myeloid-rich", which is the dominant real signal.
    "immune": (
        "PTPRC", "CD68", "AIF1", "ITGAM", "CSF1R", "TYROBP", "FCER1G",
        "C1QA", "C1QB", "CD14", "LAPTM5", "CD53", "HLA-DRA", "CD3E", "CD2",
    ),
    # Stromal / fibroblast-pericyte. Fibrillar collagens (COL1A1, COL1A2,
    # COL3A1, COL5A1, COL6A2), fibronectin (FN1), THY1 (CD90), smooth-muscle /
    # pericyte contractile genes (ACTA2, TAGLN, MYL9, PDGFRB) and the
    # small-leucine-rich proteoglycans DCN / LUM / POSTN.
    #
    # CAVEAT specific to GBM: the brain has very little true fibroblast
    # stroma. Most of this signal in a GBM biopsy comes from vascular mural
    # cells (pericytes, smooth muscle) and from leptomeningeal/dural tissue
    # caught in the specimen, so "stromal" here overlaps with "perivascular"
    # and is partially collinear with the endothelial panel.
    "stromal": (
        "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A2", "FN1", "THY1",
        "ACTA2", "TAGLN", "MYL9", "PDGFRB", "DCN", "LUM", "POSTN",
    ),
    # Endothelial. PECAM1 (CD31), VWF, CDH5 (VE-cadherin), CLDN5 (brain
    # endothelial tight junction), ESAM, ROBO4, TIE1, TEK (TIE2), KDR (VEGFR2),
    # FLT1 (VEGFR1), EMCN, PLVAP, ERG and SOX17 (endothelial transcription
    # factors). Relevant to GBM because microvascular proliferation is a
    # diagnostic feature.
    "endothelial": (
        "PECAM1", "VWF", "CDH5", "CLDN5", "ESAM", "ROBO4", "TIE1", "TEK",
        "KDR", "FLT1", "EMCN", "PLVAP", "ERG", "SOX17",
    ),
}

#: A deliberately broader and almost non-overlapping immune panel, used only to
#: cross-check ``tme_immune`` (see :func:`immune_signature_agreement`). Only
#: PTPRC is shared with the scoring panel, so a high correlation is informative
#: rather than tautological.
BROAD_IMMUNE_PANEL: tuple[str, ...] = (
    "PTPRC", "CD3D", "CD8A", "LCK", "CD48", "GZMB", "NKG7", "CCL5",
    "MS4A1", "CD79A", "ITGB2", "FCGR3A", "MNDA", "MS4A4A", "MS4A7",
    "CYBB", "CTSS", "FCGR2A", "IL10RA", "TLR1",
)

#: Marker -> state whose fraction it should be highest in.
#:
#: MEASURED CAVEAT (TCGA-GBM, n=284): GFAP is the textbook astrocytic marker
#: but it is a POOR discriminator here, because it is saturated -- median
#: 3290 TPM, among the top-ranked genes in essentially every GBM sample -- so
#: its within-cohort variance carries little AC-specific information. Measured
#: Spearman against the raw AC-like score is only 0.29, and against the
#: competitive ``state_ac`` fraction 0.11, so the top-quartile concordance
#: check FAILS for GFAP. The non-saturated astrocytic transporters do what is
#: expected: SLC1A2 r=0.79 and AQP4 r=0.58 against the AC-like score. GFAP is
#: kept in this dict deliberately so the failure stays visible rather than
#: being quietly dropped.
STATE_CONCORDANCE_MARKERS: dict[str, str] = {
    "GFAP": "state_ac",
    "AQP4": "state_ac",
    "SLC1A2": "state_ac",
    "PDGFRA": "state_opc",
    "DCX": "state_npc",
    "CHI3L1": "state_mes",
}
#: Marker -> TME compartment whose fraction it should be highest in.
TME_CONCORDANCE_MARKERS: dict[str, str] = {
    "PTPRC": "tme_immune",
    "COL1A1": "tme_stromal",
    "PECAM1": "tme_endothelial",
}


# ---------------------------------------------------------------------------
# Score -> "fraction" conversion
# ---------------------------------------------------------------------------
def scores_to_fractions(
    scores: pd.DataFrame,
    *,
    method: str = "softmax_z",
    temperature: float = 1.0,
) -> pd.DataFrame:
    """Normalise a signature-score matrix so each row sums to 1.

    **This step is where the output stops being interpretable as biology and
    starts being interpretable only as a ranking.** Both available rules are
    arbitrary monotone transforms chosen to satisfy the downstream column
    contract (fractions in [0, 1] summing to 1); neither estimates a
    proportion.

    Parameters
    ----------
    scores:
        ``samples x signatures``.
    method:
        ``"softmax_z"`` (default) standardises each signature column across the
        cohort and applies a softmax. This removes the fact that different
        marker panels have different baseline enrichment -- without it, a panel
        containing VIM and CD44 (extremely highly expressed) would dominate
        every sample. Price: the result is cohort-relative, so a sample's
        fraction changes if the cohort changes, and the cohort mean of every
        fraction is pinned near 1/n_signatures.

        ``"proportional"`` instead clips the incoming scores at zero and
        divides by the row sum, with no cohort standardisation at all. Note
        precisely what that means in composition with ``score_method``: it
        normalises *the score matrix it is handed*, not raw expression. With
        the default ``score_method="ssgsea"`` the inputs are per-sample
        enrichment scores, so the result is cohort-independent (a sample's
        fractions do not change when the cohort changes) at the cost of
        discrimination -- panels with systematically higher baseline
        enrichment keep a larger share in every sample, so the fractions are
        much more rank-constant across samples than ``softmax_z``. Rows whose
        scores are all <= 0 fall back to a uniform ``1/n_signatures``.

    Information-free samples (both methods): a sample whose every score is
    exactly 0 -- which for ssGSEA happens if and only if the sample is
    internally constant (every gene tied), e.g. an all-zero column -- has no
    ranking information at all. It is given uniform fractions. Without this,
    ``"softmax_z"`` would cohort-standardise those zeros into a confident,
    cohort-dependent "preference" for whichever signature happens to have the
    lowest cohort mean, which is an artefact and not a biological call.
    temperature:
        Softmax temperature for ``"softmax_z"``. Larger flattens.

    Returns
    -------
    DataFrame with the same index and columns, rows summing to 1.
    """
    if method not in {"softmax_z", "proportional"}:
        raise ValueError(f"unknown method {method!r}; expected 'softmax_z' or 'proportional'")
    if temperature <= 0:
        raise ValueError("temperature must be > 0")

    x = scores.to_numpy(dtype=float, copy=True)
    n_samples, n_sig = x.shape
    if n_samples == 0 or n_sig == 0:
        return pd.DataFrame(x, index=scores.index, columns=scores.columns)

    # A row of exactly-zero scores carries no ranking information (see the
    # docstring). Detected on the raw scores, before any clipping.
    uninformative = np.all(np.abs(x) <= 1e-12, axis=1)

    if method == "proportional":
        x = np.clip(x, 0.0, None)
        tot = x.sum(axis=1, keepdims=True)
        frac = np.divide(
            x, tot, out=np.full_like(x, 1.0 / n_sig), where=tot > 0
        )
        frac[uninformative] = 1.0 / n_sig
        return pd.DataFrame(frac, index=scores.index, columns=scores.columns)

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
    frac = e / e.sum(axis=1, keepdims=True)
    frac[uninformative] = 1.0 / n_sig
    return pd.DataFrame(frac, index=scores.index, columns=scores.columns)


class EnrichmentResult:
    """Score matrix plus mapping report for a panel set."""

    def __init__(
        self,
        scores: pd.DataFrame,
        mappings: dict[str, SignatureMapping],
        method: str,
    ) -> None:
        self.scores = scores
        self.mappings = mappings
        self.method = method

    def report(self) -> pd.DataFrame:
        return mapping_report(self.mappings)

    @property
    def overall_mapping_rate(self) -> float:
        req = sum(m.n_requested for m in self.mappings.values())
        fnd = sum(m.n_found for m in self.mappings.values())
        return (fnd / req) if req else 0.0


def _score(
    expression: pd.DataFrame,
    annotation: pd.DataFrame | None,
    signatures: dict[str, tuple[str, ...]],
    score_method: str,
) -> EnrichmentResult:
    scores, mappings = score_signatures(
        expression, signatures, annotation=annotation, method=score_method
    )
    return EnrichmentResult(scores=scores, mappings=mappings, method=score_method)


def state_score_matrix(
    expression: pd.DataFrame,
    annotation: pd.DataFrame | None = None,
    *,
    score_method: str = "ssgsea",
) -> EnrichmentResult:
    """Raw Neftel meta-module enrichment scores, before fraction normalisation.

    Inspect this rather than the fractions when asking whether a sample is
    genuinely AC-high or merely AC-highest.
    """
    return _score(expression, annotation, NEFTEL_SIGNATURES, score_method)


def tme_score_matrix(
    expression: pd.DataFrame,
    annotation: pd.DataFrame | None = None,
    *,
    score_method: str = "ssgsea",
) -> EnrichmentResult:
    """Raw immune / stromal / endothelial enrichment scores, before normalisation."""
    return _score(expression, annotation, TME_SIGNATURES, score_method)


def _assemble(
    frac: pd.DataFrame,
    column_map: dict[str, str],
    *,
    method_col: str,
    version_col: str,
) -> pd.DataFrame:
    out = pd.DataFrame({"patient_id": pd.Index(frac.index).astype(str)})
    for sig_name, col in column_map.items():
        out[col] = frac[sig_name].to_numpy() if len(frac) else pd.Series(dtype=float)
    out[method_col] = METHOD_NAME
    out[version_col] = METHOD_VERSION
    out["is_true_deconvolution"] = False
    return out.reset_index(drop=True)


def deconvolve_states_expression(
    expression: pd.DataFrame,
    annotation: pd.DataFrame | None = None,
    *,
    score_method: str = "ssgsea",
    fraction_method: str = "softmax_z",
    temperature: float = 1.0,
    seed: int = 42,
) -> pd.DataFrame:
    """Neftel state "fractions" per patient from bulk expression.

    **The returned ``state_*`` values are marker-enrichment scores normalised to
    sum to 1, not cell-state proportions.** No reference profile is fitted and
    no mixing model is solved, so these numbers do not estimate what fraction
    of the sample is NPC-like, OPC-like, AC-like or MES-like. They rank states
    within a sample and samples within the supplied cohort. This method has
    never been compared to ground truth (none exists in this project), so its
    agreement with true state fractions is unknown, and protocol test T1
    (state-deconvolution Pearson r >= 0.75) is NOT addressed.

    Parameters
    ----------
    expression:
        ``gene x sample`` TPM matrix.
    annotation:
        gene_annotation table; ``None`` means the matrix is symbol-indexed.
    score_method:
        ``"ssgsea"`` or ``"zmean"`` (see
        :func:`gliostrata.profiling.verhaak.score_signatures`).
    fraction_method:
        ``"softmax_z"`` or ``"proportional"`` (see :func:`scores_to_fractions`).
    temperature:
        Softmax temperature.
    seed:
        Ignored; this function is fully deterministic. Accepted only for
        signature parity with ``gliostrata.profiling.states.deconvolve_states``.

    Returns
    -------
    DataFrame with ``patient_id``, ``state_npc``, ``state_opc``, ``state_ac``,
    ``state_mes`` (each in [0, 1], summing to 1 per row) plus method metadata.
    """
    del seed
    res = state_score_matrix(expression, annotation, score_method=score_method)
    frac = scores_to_fractions(res.scores, method=fraction_method, temperature=temperature)
    out = _assemble(
        frac,
        _STATE_TO_COLUMN,
        method_col="state_deconvolution_method",
        version_col="state_deconvolution_version",
    )
    return out


def deconvolve_tme_expression(
    expression: pd.DataFrame,
    annotation: pd.DataFrame | None = None,
    *,
    score_method: str = "ssgsea",
    fraction_method: str = "softmax_z",
    temperature: float = 1.0,
    seed: int = 42,
) -> pd.DataFrame:
    """Immune / stromal / endothelial "fractions" per patient from bulk expression.

    **Same caveat as :func:`deconvolve_states_expression`: these are normalised
    marker-enrichment scores, not cell-type proportions.** They are not
    CIBERSORTx or MCP-counter output and have not been compared to either.
    In particular ``tme_stromal`` in brain tissue is substantially vascular
    mural / meningeal signal and is partly collinear with ``tme_endothelial``
    (see ``TME_SIGNATURES``).

    Returns
    -------
    DataFrame with ``patient_id``, ``tme_immune``, ``tme_stromal``,
    ``tme_endothelial`` (each in [0, 1], summing to 1 per row) plus metadata.
    """
    del seed
    res = tme_score_matrix(expression, annotation, score_method=score_method)
    frac = scores_to_fractions(res.scores, method=fraction_method, temperature=temperature)
    return _assemble(
        frac,
        _TME_TO_COLUMN,
        method_col="tme_deconvolution_method",
        version_col="tme_deconvolution_version",
    )


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------
def check_fractions(
    frame: pd.DataFrame,
    columns: tuple[str, ...],
    *,
    tol: float = 1e-9,
) -> dict:
    """Assert the compositional contract and return a summary dict.

    Checks every value is finite and in [0, 1] and every row sums to 1 within
    ``tol``. Raises ``ValueError`` on violation -- a broken composition is a
    bug, not a warning.
    """
    missing = [c for c in columns if c not in frame.columns]
    if missing:
        raise ValueError(f"frame is missing fraction columns: {missing}")
    x = frame[list(columns)].to_numpy(dtype=float)
    if x.size == 0:
        return {"n": 0, "min": None, "max": None, "max_abs_sum_error": 0.0}
    if not np.isfinite(x).all():
        raise ValueError("fraction columns contain NaN/inf")
    if x.min() < -tol or x.max() > 1.0 + tol:
        raise ValueError(f"fractions out of [0, 1]: min={x.min()}, max={x.max()}")
    sums = x.sum(axis=1)
    err = float(np.abs(sums - 1.0).max())
    if err > max(tol, 1e-9):
        raise ValueError(f"fraction rows do not sum to 1 (max abs error {err})")
    return {
        "n": int(x.shape[0]),
        "min": float(x.min()),
        "max": float(x.max()),
        "max_abs_sum_error": err,
    }


def fraction_distribution(frame: pd.DataFrame, columns: tuple[str, ...]) -> pd.DataFrame:
    """Per-column distribution summary (mean, sd, quartiles, argmax counts)."""
    rows = []
    if len(frame) == 0:
        return pd.DataFrame(
            columns=["column", "mean", "sd", "min", "q25", "median", "q75", "max", "n_argmax"]
        )
    sub = frame[list(columns)].astype(float)
    argmax = sub.to_numpy().argmax(axis=1)
    for j, c in enumerate(columns):
        v = sub[c]
        rows.append(
            {
                "column": c,
                "mean": float(v.mean()),
                "sd": float(v.std(ddof=0)),
                "min": float(v.min()),
                "q25": float(v.quantile(0.25)),
                "median": float(v.median()),
                "q75": float(v.quantile(0.75)),
                "max": float(v.max()),
                "n_argmax": int((argmax == j).sum()),
            }
        )
    return pd.DataFrame(rows)


def marker_fraction_concordance(
    expression: pd.DataFrame,
    fractions: pd.DataFrame,
    markers: dict[str, str],
    annotation: pd.DataFrame | None = None,
    *,
    top_frac: float = 0.25,
) -> pd.DataFrame:
    """Do marker-high samples get the expected fraction?

    For each marker gene, take the top ``top_frac`` of samples by that marker's
    expression and report the mean of the expected fraction column in that group
    versus the rest, plus the Spearman correlation between marker expression and
    the expected fraction across all samples.

    Internal consistency only: the marker is a member of the panel that
    produced the fraction. It can falsify the implementation (a bug in gene
    mapping or normalisation would break it) but cannot validate the biology.
    """
    sym = to_symbol_matrix(expression, annotation)
    sym.columns = sym.columns.astype(str)
    if len(fractions) == 0:
        return pd.DataFrame(
            columns=[
                "marker", "expected_column", "n_top", "mean_in_marker_high",
                "mean_in_rest", "spearman_r", "concordant",
            ]
        )

    pid = fractions["patient_id"].astype(str).tolist()
    keep = [c for c in sym.columns if c in set(pid)]
    f = fractions.set_index(fractions["patient_id"].astype(str))

    rows = []
    for marker, col in markers.items():
        if marker not in sym.index or not keep:
            rows.append(
                {
                    "marker": marker, "expected_column": col, "n_top": 0,
                    "mean_in_marker_high": np.nan, "mean_in_rest": np.nan,
                    "spearman_r": np.nan, "concordant": False,
                }
            )
            continue
        expr = sym.loc[marker, keep].astype(float)
        frac_col = f.loc[keep, col].astype(float)
        n_top = max(1, int(round(len(keep) * top_frac)))
        top_ids = expr.sort_values(ascending=False).index[:n_top]
        mask = frac_col.index.isin(top_ids)
        hi = float(frac_col[mask].mean())
        rest = float(frac_col[~mask].mean()) if (~mask).any() else np.nan
        r = (
            float(pd.Series(expr.to_numpy()).corr(pd.Series(frac_col.to_numpy()), method="spearman"))
            if len(keep) > 2
            else np.nan
        )
        rows.append(
            {
                "marker": marker,
                "expected_column": col,
                "n_top": int(n_top),
                "mean_in_marker_high": hi,
                "mean_in_rest": rest,
                "spearman_r": r,
                "concordant": bool(np.isnan(rest) or hi > rest),
            }
        )
    return pd.DataFrame(rows)


def immune_signature_agreement(
    expression: pd.DataFrame,
    tme: pd.DataFrame,
    annotation: pd.DataFrame | None = None,
    *,
    score_method: str = "ssgsea",
) -> dict:
    """Correlate ``tme_immune`` with an independent, broader immune panel.

    ``BROAD_IMMUNE_PANEL`` shares only ``PTPRC`` with the scoring panel, so a
    strong Spearman/Pearson correlation is real convergent evidence that the
    immune axis tracks leukocyte content -- as far as marker expression can
    show it. It is still not a measurement of immune cell fraction.
    """
    scores, mappings = score_signatures(
        expression,
        {"broad_immune": BROAD_IMMUNE_PANEL},
        annotation=annotation,
        method=score_method,
    )
    if len(tme) == 0 or len(scores) == 0:
        return {"n": 0, "pearson_r": float("nan"), "spearman_r": float("nan"),
                "mapping_rate": mappings["broad_immune"].mapping_rate}
    s = scores["broad_immune"]
    s.index = s.index.astype(str)
    t = tme.set_index(tme["patient_id"].astype(str))["tme_immune"].astype(float)
    common = [i for i in s.index if i in set(t.index)]
    a = pd.Series(s.loc[common].to_numpy())
    b = pd.Series(t.loc[common].to_numpy())
    return {
        "n": len(common),
        "pearson_r": float(a.corr(b, method="pearson")) if len(common) > 2 else float("nan"),
        "spearman_r": float(a.corr(b, method="spearman")) if len(common) > 2 else float("nan"),
        "mapping_rate": mappings["broad_immune"].mapping_rate,
    }
