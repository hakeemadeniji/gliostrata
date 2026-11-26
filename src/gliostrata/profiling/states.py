"""L2 — Neftel single-cell state deconvolution.

Stub: real implementation uses CIBERSORTx (web tool, import results) or MuSiC2 (R)
to deconvolve bulk transcriptomes into NPC-like / OPC-like / AC-like / MES-like
state fractions.

The stub reads existing state columns from the cohort (which `generate_synthetic_cohort`
populates) and adds a reconciliation metadata column.
"""
from __future__ import annotations

import pandas as pd

STATES = ["state_npc", "state_opc", "state_ac", "state_mes"]


def deconvolve_states(cohort: pd.DataFrame, *, seed: int = 42) -> pd.DataFrame:
    """Deconvolve bulk transcriptomes into Neftel single-cell state fractions.

    Returns a DataFrame with patient_id and the four state fractions.
    """
    missing = [c for c in STATES if c not in cohort.columns]
    if missing:
        raise ValueError(f"Input cohort missing state columns: {missing}. "
                         "Use generate_synthetic_cohort() or implement real deconvolution.")

    out = cohort[["patient_id", *STATES]].copy()
    out["state_deconvolution_method"] = "stub_passthrough"
    out["state_deconvolution_version"] = "0.1.0-stub"
    return out


def cibersortx_signature_matrix() -> pd.DataFrame:
    """STUB — build the Neftel signature matrix for CIBERSORTx.

    Real implementation loads the Neftel 2019 atlas from GEO, normalizes,
    and writes a signature matrix file uploadable to the CIBERSORTx web portal.
    """
    raise NotImplementedError("CIBERSORTx signature matrix builder not implemented. "
                              "See docs/roadmap.md.")
