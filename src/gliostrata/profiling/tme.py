"""L2 — Tumor microenvironment deconvolution.

Stub: real implementation uses CIBERSORTx with LM22 immune signature + MCP-counter
for stromal/endothelial fractions.
"""
from __future__ import annotations

import pandas as pd

TME_FRACTIONS = ["tme_immune", "tme_stromal", "tme_endothelial"]


def deconvolve_tme(cohort: pd.DataFrame, *, seed: int = 42) -> pd.DataFrame:
    """Deconvolve bulk transcriptomes into TME fractions (immune, stromal, endothelial)."""
    missing = [c for c in TME_FRACTIONS if c not in cohort.columns]
    if missing:
        raise ValueError(f"Input cohort missing TME columns: {missing}.")

    out = cohort[["patient_id", *TME_FRACTIONS]].copy()
    out["tme_deconvolution_method"] = "stub_passthrough"
    out["tme_deconvolution_version"] = "0.1.0-stub"
    return out
