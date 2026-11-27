"""L6 — Bayesian therapeutic configuration optimizer.

Currently implements random search over the configuration space. The real
implementation will use BoTorch (Bayesian optimization with a Gaussian process
surrogate + expected improvement acquisition function) — but random search
is sufficient for the scaffold to be runnable end-to-end and is a valid
baseline for the T4 "optimizer beats random search" test (where BoTorch
should beat THIS random search).

Public API:
    optimize_configurations(feature_store, survival_pred, n_top=5, seed=42)
        → DataFrame with one row per (patient, ranked config)
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from gliostrata.optimization.config_space import (
    Configuration,
    sample_random_configuration,
)
from gliostrata.prediction.survival import predict_survival_with_config


def _config_to_row(config: Configuration, patient_id: str, rank: int, score: float) -> dict:
    return {
        "patient_id": patient_id,
        "config_id": f"{patient_id}_c{rank}",
        "rank": rank,
        "score": score,
        "cell_type": config.cell_type,
        "payload": config.payload,
        "dose_cells": config.dose_cells,
        "n_injection_sites": config.n_injection_sites,
        "schedule": config.schedule,
    }


def _score_configurations(
    feature_store: pd.DataFrame,
    survival_pred: pd.DataFrame,
    configs: list[Configuration],
    seed: int = 42,
) -> pd.DataFrame:
    """Score a list of configurations across all patients.

    Returns a long-format DataFrame: one row per (patient, config), with the
    predicted survival with that configuration + survival gain vs Stupp.
    """
    rows = []
    for cid, cfg in enumerate(configs):
        for pid in feature_store["patient_id"]:
            rows.append({
                "config_id": f"cfg_{cid}",
                "cell_type": cfg.cell_type,
                "payload": cfg.payload,
                "dose_cells": cfg.dose_cells,
                "n_injection_sites": cfg.n_injection_sites,
                "schedule": cfg.schedule,
                "patient_id": pid,
                "rank": 0,  # placeholder — ranking happens after scoring
            })
    config_df = pd.DataFrame(rows)
    scored = predict_survival_with_config(feature_store, survival_pred, config_df)
    scored["objective"] = scored["survival_gain_vs_stupp_months"]
    return scored


def optimize_configurations(
    feature_store: pd.DataFrame,
    survival_pred: pd.DataFrame,
    *,
    n_top: int = 5,
    n_candidates: int = 100,
    seed: int = 42,
) -> pd.DataFrame:
    """Run random-search optimization for each patient.

    For each patient, samples `n_candidates` random configurations from the
    configuration space, scores each using the L5 survival model + the L4
    tumor-kill proxy, and returns the top-`n_top` configurations per patient
    ranked by `objective` (predicted survival gain vs Stupp).

    Returns a DataFrame with columns:
        patient_id, config_id, rank, score, cell_type, payload, dose_cells,
        n_injection_sites, schedule, predicted_survival_with_config,
        predicted_tumor_kill_fraction, survival_gain_vs_stupp_months
    """
    rng = np.random.default_rng(seed)

    # Sample candidate configurations once (shared across patients for efficiency)
    candidates = [sample_random_configuration(rng) for _ in range(n_candidates)]

    # Score every (patient, candidate) pair
    scored = _score_configurations(feature_store, survival_pred, candidates, seed=seed)

    # Rank per patient by objective
    out_rows = []
    for pid, group in scored.groupby("patient_id"):
        top = group.nlargest(n_top, "objective")
        for rank_i, (_, row) in enumerate(top.iterrows(), start=1):
            out_rows.append({
                "patient_id": pid,
                "config_id": f"{pid}_r{rank_i}",
                "rank": rank_i,
                "score": float(row["objective"]),
                "cell_type": row["cell_type"],
                "payload": row["payload"],
                "dose_cells": int(row["dose_cells"]),
                "n_injection_sites": int(row["n_injection_sites"]),
                "schedule": row["schedule"],
                "predicted_survival_with_config": float(row["predicted_survival_with_config"]),
                "predicted_tumor_kill_fraction": float(row["predicted_tumor_kill_fraction"]),
                "survival_gain_vs_stupp_months": float(row["survival_gain_vs_stupp_months"]),
            })

    out = pd.DataFrame(out_rows)
    out.attrs["n_candidates"] = n_candidates
    out.attrs["n_top"] = n_top
    out.attrs["optimizer"] = "random_search"
    out.attrs["optimizer_version"] = "0.1.0"
    return out


def evaluate_personalization(configs: pd.DataFrame) -> dict:
    """Evaluate the personalization-variance check (T4 threshold).

    Returns a dict with: n_patients, n_distinct_top5_payloads, n_distinct_top5_celltypes,
    passes_threshold (bool, ≥ 3 distinct payload/site combinations in top-5 across 50 patients).
    """
    # Take top-5 per patient
    top5 = configs[configs["rank"] <= 5]
    # Distinct (cell_type, payload) combos across all top-5
    distinct_combos = top5[["cell_type", "payload"]].drop_duplicates()
    return {
        "n_patients": configs["patient_id"].nunique(),
        "n_distinct_top5_cell_payload_combos": len(distinct_combos),
        "n_distinct_top5_payloads": top5["payload"].nunique(),
        "n_distinct_top5_celltypes": top5["cell_type"].nunique(),
        "threshold_min_distinct": 3,
        "passes_threshold": len(distinct_combos) >= 3,
    }
