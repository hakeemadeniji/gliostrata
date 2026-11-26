"""L4 — Tumor growth simulation (PhysiCell wrapper).

Stub: real implementation wraps the PhysiCell C++ agent-based simulator
(Ghaffarizadeh et al. 2018) via pyMCDS for output parsing, calibrated against
≥10 published PDX growth curves.

The stub `simulate_tumor_growth` below returns a synthetic growth curve
(Gompertzian kinetics) parameterized by the patient's tumor profile.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def simulate_tumor_growth(
    patient_features: pd.Series | dict,
    *,
    hours: int = 720,
    dt_hours: float = 1.0,
    seed: int = 42,
) -> pd.DataFrame:
    """Simulate tumor growth over time for a single patient.

    Returns a DataFrame with columns: time_hours, tumor_volume_mm3, n_live_cells.

    Stub uses Gompertzian growth: V(t) = V0 * exp(K/b * (1 - exp(-b*t)))
    where K and b depend on subtype (MES grows fastest) and IDH status (IDH-mut slower).
    """
    rng = np.random.default_rng(seed)

    # Growth rate depends on subtype
    subtype = patient_features.get("subtype", "Classical")
    k_growth = {
        "Classical": 0.020,
        "Mesenchymal": 0.028,  # fastest
        "Proneural": 0.015,
        "Neural": 0.018,
    }.get(subtype, 0.020)

    # IDH-mutant tumors grow slower
    if patient_features.get("idh_mutated", False):
        k_growth *= 0.6

    b = 0.005  # retardation
    v0 = 50.0  # mm^3 initial

    t = np.arange(0, hours + 1, dt_hours)
    # Gompertz: V(t) = V0 * exp((K/b) * (1 - exp(-b*t)))
    vol = v0 * np.exp((k_growth / b) * (1 - np.exp(-b * t)))
    # Add observation noise
    vol = vol * rng.normal(1.0, 0.05, size=len(t))
    n_cells = vol * 1e6  # ~1e6 cells per mm^3

    return pd.DataFrame({
        "time_hours": t,
        "tumor_volume_mm3": vol,
        "n_live_cells": n_cells,
    })


def calibrate_against_pdx(pdx_curves: list[pd.DataFrame]) -> dict:
    """STUB — calibrate PhysiCell growth parameters against PDX growth curves.

    Real implementation:
      1. Load ≥10 PDX series from published supplementary data.
      2. Define PhysiCell parameter priors.
      3. Use ABC or Bayesian inference to fit parameters per PDX series.
      4. Cross-validate: train on half the PDX series, test on the other half.
      5. Return posterior distributions over key parameters.

    Returns stub metadata.
    """
    return {
        "calibrated": False,
        "n_pdx_series_used": len(pdx_curves),
        "rmse_pct": None,
        "posterior_summary": {},
        "note": "Stub — real calibration uses PhysiCell + ABC. See docs/roadmap.md.",
    }
