"""L5 — Survival prediction model.

Stub: real implementation trains a DeepSurv (Cox-PH neural network) on the L3
feature store with TCGA as the training cohort, validated on frozen CGGA + GLASS
test partitions. Calibration via Brier score + decision-curve analysis.

The stub `predict_survival` below returns a per-patient predicted survival
(months) under standard Stupp protocol, using a simple Cox-PH-like linear model
on the log-hazard. This makes the optimizer end-to-end runnable.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# Baseline log-hazard coefficients (illustrative — calibrated to make the stub
# produce survival in the right range. Real implementation learns these from data.)
HAZARD_COEFFS = {
    "age": 0.015,
    "kps": -0.015,
    "sex_male": 0.05,
    "extent_of_resection_STR": 0.20,
    "extent_of_resection_GTR": -0.10,
    "mgmt_methylated": -0.30,
    "idh_mutated": -1.00,
    "subtype_prob_Classical": 0.0,
    "subtype_prob_Mesenchymal": 0.25,
    "subtype_prob_Proneural": -0.30,
    "subtype_prob_Neural": 0.05,
    "state_npc": -0.20,
    "state_opc": -0.10,
    "state_ac": 0.10,
    "state_mes": 0.30,
    "tme_immune": -0.20,
    "tme_stromal": 0.10,
    "tme_endothelial": 0.05,
}

BASELINE_LOG_HAZARD = np.log(1 / 14.0)  # base hazard ~ 1/14 months median


def predict_survival(feature_store: pd.DataFrame, *, seed: int = 42) -> pd.DataFrame:
    """Predict per-patient overall survival under standard Stupp.

    Returns a DataFrame with: patient_id, predicted_survival_months,
    predicted_log_hazard, baseline_stupp_survival_months.

    The `predicted_survival_months` is the model's prediction; the
    `baseline_stupp_survival_months` is a literature-based Stupp baseline
    (currently the same value, but separated so the optimizer can compute
    "survival gain vs Stupp" cleanly).
    """
    rng = np.random.default_rng(seed)
    feat_cols = list(HAZARD_COEFFS.keys())

    # Ensure all feature columns are present
    missing = [c for c in feat_cols if c not in feature_store.columns]
    if missing:
        raise ValueError(f"Feature store missing columns: {missing}")

    X = feature_store[feat_cols].fillna(0).values.astype(float)
    coeffs = np.array([HAZARD_COEFFS[c] for c in feat_cols])

    # Log-hazard: linear predictor + baseline + noise
    log_hazard = BASELINE_LOG_HAZARD + X @ coeffs + rng.normal(0, 0.1, size=len(X))
    # Convert to predicted survival (median) under exponential model: S(t) = exp(-h*t)
    # → median survival = ln(2) / h
    predicted_survival = np.log(2) / np.exp(log_hazard)

    out = pd.DataFrame({
        "patient_id": feature_store["patient_id"].values,
        "predicted_survival_months": predicted_survival,
        "predicted_log_hazard": log_hazard,
        "baseline_stupp_survival_months": predicted_survival,  # stub: same as predicted
    })
    out["survival_model_version"] = "0.1.0-stub-cox-linear"
    return out


def predict_survival_with_config(
    feature_store: pd.DataFrame,
    survival_pred: pd.DataFrame,
    config_df: pd.DataFrame,
) -> pd.DataFrame:
    """Predict survival under a specific therapeutic configuration (used by L6).

    For each (patient, configuration), predicts survival given the configuration's
    simulated tumor kill and the patient's baseline.

    Stub model: survival = baseline_survival * (1 + tumor_kill_fraction * response_modifier)
    where response_modifier depends on patient features and payload match.
    """
    df = config_df.merge(feature_store, on="patient_id", how="left") \
                  .merge(survival_pred, on="patient_id", how="left", suffixes=("", "_pred"))

    # Compute response modifier — depends on (payload, patient subtype) match
    payload_subtype_affinity = {
        "CRAd":            {"Classical": 0.6, "Mesenchymal": 0.5, "Proneural": 0.7, "Neural": 0.6},
        "prodrug_enzyme":  {"Classical": 0.5, "Mesenchymal": 0.7, "Proneural": 0.5, "Neural": 0.5},
        "CAR":             {"Classical": 0.4, "Mesenchymal": 0.6, "Proneural": 0.4, "Neural": 0.4},
        "combo":           {"Classical": 0.8, "Mesenchymal": 0.8, "Proneural": 0.8, "Neural": 0.8},
    }
    subtype = df["subtype"].values
    payload = df["payload"].values
    affinity = np.array([
        payload_subtype_affinity[p][s] for p, s in zip(payload, subtype)
    ])

    # Tumor kill fraction from L4 (synthetic): scales with dose × sites × cell_tropism
    cell_tropism = df["cell_type"].map({
        "iPSC-NSC": 0.85, "MSC": 0.70, "HB1_F3_line": 0.78
    }).values

    # Normalize dose to [0, 1] over [1e6, 5e7]
    dose_norm = np.clip((df["dose_cells"].values - 1e6) / (5e7 - 1e6), 0, 1)
    # Diminishing returns on injection sites
    sites_eff = 1 - np.exp(-0.4 * df["n_injection_sites"].values)
    # Schedule multiplier
    sched_mult = df["schedule"].map({
        "single": 1.0, "fractionated_x2": 1.15, "fractionated_x3": 1.25
    }).values

    tumor_kill = np.clip(
        dose_norm * sites_eff * cell_tropism * sched_mult * affinity,
        0, 0.95
    )

    # Survival = baseline * (1 + 2 * tumor_kill) — capped so huge tumor kills don't explode
    predicted_with_config = df["predicted_survival_months"].values * (1 + 2 * tumor_kill)
    # Toxicity penalty: high doses + multiple fractions increase harm risk
    tox_penalty = (
        0.1 * dose_norm * df["schedule"].map({
            "single": 1.0, "fractionated_x2": 1.5, "fractionated_x3": 2.0
        }).values
    )
    predicted_with_config = predicted_with_config * (1 - tox_penalty)

    out_cols = ["patient_id", "config_id", "rank",
                "cell_type", "payload", "dose_cells", "n_injection_sites", "schedule"]
    out_cols = [c for c in out_cols if c in df.columns]
    out = df[out_cols].copy()
    out["predicted_survival_with_config"] = predicted_with_config
    out["predicted_tumor_kill_fraction"] = tumor_kill
    out["predicted_toxicity_penalty"] = tox_penalty
    out["survival_gain_vs_stupp_months"] = (
        predicted_with_config - df["baseline_stupp_survival_months"].values
    )
    return out
