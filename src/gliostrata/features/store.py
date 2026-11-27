"""L3 — Unified patient feature store.

Assembles per-patient feature vectors from L2 profiling outputs (subtype, states, TME)
plus L0 clinical fields. This module is **working** (not a stub) — it joins the inputs
on patient_id and produces a feature matrix suitable for downstream survival modeling
and optimization.
"""
from __future__ import annotations

import pandas as pd

# Feature columns used by L5 (survival model) and L6 (optimizer) as input.
# Anything not in this list is treated as metadata, not a model feature.
FEATURE_COLUMNS = [
    # Demographics
    "age", "kps",
    "sex_male",  # one-hot encoded
    # Treatment
    "extent_of_resection_STR", "extent_of_resection_GTR",  # biopsy is reference
    # Molecular
    "mgmt_methylated", "idh_mutated",
    # Subtype probabilities (continuous, not argmax — preserves uncertainty)
    "subtype_prob_Classical", "subtype_prob_Mesenchymal",
    "subtype_prob_Proneural", "subtype_prob_Neural",
    # Neftel state fractions
    "state_npc", "state_opc", "state_ac", "state_mes",
    # TME fractions
    "tme_immune", "tme_stromal", "tme_endothelial",
]


def build_feature_store(
    cohort: pd.DataFrame,
    subtypes: pd.DataFrame,
    states: pd.DataFrame,
    tme: pd.DataFrame,
) -> pd.DataFrame:
    """Assemble the unified patient feature store.

    Joins L0 cohort + L2 subtype + L2 states + L2 TME outputs on patient_id,
    one-hot encodes categorical columns, and selects only FEATURE_COLUMNS.

    The output is a per-patient feature matrix with the same row order as
    `cohort`. This is the input to L5 (survival) and L6 (optimizer).
    """
    # Merge everything on patient_id
    df = cohort.merge(subtypes, on="patient_id", how="inner", suffixes=("", "_sub"))
    df = df.merge(states, on="patient_id", how="inner", suffixes=("", "_states"))
    df = df.merge(tme, on="patient_id", how="inner", suffixes=("", "_tme"))

    # Resolve any duplicate column suffixes (preferring the original cohort values)
    for col in ["state_npc", "state_opc", "state_ac", "state_mes",
                "tme_immune", "tme_stromal", "tme_endothelial"]:
        # If the merge created a _states or _tme suffixed duplicate, keep the original
        # cohort column (already merged in from the L2 outputs as well).
        if col not in df.columns and f"{col}_states" in df.columns:
            df[col] = df[f"{col}_states"]
        if col not in df.columns and f"{col}_tme" in df.columns:
            df[col] = df[f"{col}_tme"]

    # One-hot encode categoricals
    df["sex_male"] = (df["sex"] == "M").astype(int)
    for level in ["STR", "GTR"]:
        df[f"extent_of_resection_{level}"] = (df["extent_of_resection"] == level).astype(int)

    # Ensure booleans are int
    for c in ["mgmt_methylated", "idh_mutated"]:
        if c in df.columns:
            df[c] = df[c].astype(int)

    # Select feature columns + keep patient_id + outcome columns for downstream use
    keep = ["patient_id"] + FEATURE_COLUMNS + [
        "overall_survival_months", "vital_status", "treatment_arm", "subtype",
    ]
    keep = [c for c in keep if c in df.columns]
    out = df[keep].copy()
    out.attrs["feature_columns"] = FEATURE_COLUMNS
    out.attrs["version"] = "0.1.0"
    return out


def feature_matrix(feature_store: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split the feature store into (X, y_meta) for modeling.

    X = feature matrix (FEATURE_COLUMNS only, no patient_id).
    y_meta = DataFrame with patient_id, overall_survival_months, vital_status.
    """
    feat_cols = feature_store.attrs.get("feature_columns", FEATURE_COLUMNS)
    X = feature_store[feat_cols].copy()
    y_meta = feature_store[["patient_id", "overall_survival_months", "vital_status"]].copy()
    return X, y_meta
