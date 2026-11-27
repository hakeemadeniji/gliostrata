"""L4 — NSC migration simulation.

Stub: real implementation models NSC migration as a chemotaxis problem — NSCs
are agents whose velocity depends on local gradients of SDF-1, VEGF, uPA, MCP-1,
with chemoattractant fields estimated from tumor transcriptomic signatures.
Validated against published HB1.F3 luciferase imaging.

The stub `simulate_nsc_migration` returns synthetic migration statistics that
depend on cell type and patient features.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# Cell-type-specific migration parameters (illustrative, calibrated to literature order-of-magnitude)
CELL_MIGRATION_PARAMS = {
    "iPSC-NSC":      {"velocity_um_per_h": 45, "tropism_score": 0.85, "half_life_h": 48},
    "MSC":           {"velocity_um_per_h": 28, "tropism_score": 0.70, "half_life_h": 72},
    "HB1_F3_line":   {"velocity_um_per_h": 35, "tropism_score": 0.78, "half_life_h": 60},
}


def simulate_nsc_migration(
    patient_features: pd.Series | dict,
    cell_type: str,
    *,
    hours: int = 168,
    n_cells_injected: int = 2_000_000,
    seed: int = 42,
) -> dict:
    """Simulate NSC migration for a single patient + cell type.

    Returns a dict with summary statistics: mean_distance_um, fraction_at_tumor,
    half_life_h, cell_type, n_cells_injected.

    Stub uses a simple model: cell-type-specific migration velocity × time,
    scaled by a patient-specific tropism modifier.
    """
    rng = np.random.default_rng(seed)

    if cell_type not in CELL_MIGRATION_PARAMS:
        raise ValueError(f"Unknown cell_type: {cell_type}. "
                         f"Choose from {list(CELL_MIGRATION_PARAMS)}.")

    params = CELL_MIGRATION_PARAMS[cell_type]
    velocity = params["velocity_um_per_h"]
    tropism = params["tropism_score"]

    # Patient modifier: MES subtype secretes more chemoattractants → higher tropism
    subtype = patient_features.get("subtype", "Classical")
    tropism_modifier = {
        "Mesenchymal": 1.20,
        "Classical": 1.00,
        "Neural": 0.95,
        "Proneural": 0.85,
    }.get(subtype, 1.0)
    tropism_patient = min(1.0, tropism * tropism_modifier)

    # Mean distance migrated (with stochastic spread)
    mean_dist = velocity * hours * rng.normal(0.9, 0.1)
    fraction_at_tumor = min(1.0, tropism_patient * rng.normal(0.95, 0.05))

    return {
        "cell_type": cell_type,
        "mean_distance_um": float(mean_dist),
        "fraction_at_tumor": float(fraction_at_tumor),
        "half_life_h": params["half_life_h"],
        "n_cells_injected": n_cells_injected,
        "hours_simulated": hours,
        "simulator_version": "0.1.0-stub",
    }
