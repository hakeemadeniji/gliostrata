"""L6 — Therapeutic configuration space definition.

Defines the joint configuration space the Bayesian optimizer searches:
    - cell_type:        iPSC-NSC / MSC / HB1_F3_line
    - payload:          CRAd / prodrug_enzyme / CAR / combo
    - dose_cells:       1e6 to 5e7 (continuous)
    - n_injection_sites: 1 to 8 (discrete)
    - schedule:         single / fractionated_x2 / fractionated_x3

Also defines the constraints and feasibility checks.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd


CELL_TYPES = ["iPSC-NSC", "MSC", "HB1_F3_line"]
PAYLOADS = ["CRAd", "prodrug_enzyme", "CAR", "combo"]
SCHEDULES = ["single", "fractionated_x2", "fractionated_x3"]

DOSE_MIN = 1_000_000
DOSE_MAX = 50_000_000
SITES_MIN = 1
SITES_MAX = 8


@dataclass
class Configuration:
    """A single therapeutic configuration."""

    cell_type: str
    payload: str
    dose_cells: int
    n_injection_sites: int
    schedule: str

    def __post_init__(self) -> None:
        validate_configuration(self)

    def to_dict(self) -> dict[str, Any]:
        return {
            "cell_type": self.cell_type,
            "payload": self.payload,
            "dose_cells": self.dose_cells,
            "n_injection_sites": self.n_injection_sites,
            "schedule": self.schedule,
        }


def validate_configuration(config: Configuration) -> None:
    """Raise ValueError if the configuration violates a constraint."""
    if config.cell_type not in CELL_TYPES:
        raise ValueError(f"cell_type {config.cell_type!r} not in {CELL_TYPES}")
    if config.payload not in PAYLOADS:
        raise ValueError(f"payload {config.payload!r} not in {PAYLOADS}")
    if not (DOSE_MIN <= config.dose_cells <= DOSE_MAX):
        raise ValueError(f"dose_cells {config.dose_cells} outside [{DOSE_MIN}, {DOSE_MAX}]")
    if not (SITES_MIN <= config.n_injection_sites <= SITES_MAX):
        raise ValueError(f"n_injection_sites {config.n_injection_sites} outside "
                         f"[{SITES_MIN}, {SITES_MAX}]")
    if config.schedule not in SCHEDULES:
        raise ValueError(f"schedule {config.schedule!r} not in {SCHEDULES}")
    # Cell type × payload compatibility
    incompatible = {
        ("HB1_F3_line", "CAR"),      # HB1.F3 line is for viral payload delivery, not CAR
        ("MSC", "CRAd"),             # MSCs are not typically CRAd carriers
    }
    if (config.cell_type, config.payload) in incompatible:
        raise ValueError(f"Incompatible cell_type × payload: "
                         f"{config.cell_type} × {config.payload}")


def sample_random_configuration(rng: np.random.Generator) -> Configuration:
    """Sample a uniformly random valid configuration."""
    # Retry until we hit a valid (cell_type, payload) combo
    for _ in range(20):
        cell_type = rng.choice(CELL_TYPES)
        payload = rng.choice(PAYLOADS)
        try:
            return Configuration(
                cell_type=str(cell_type),
                payload=str(payload),
                dose_cells=int(rng.integers(DOSE_MIN, DOSE_MAX + 1)),
                n_injection_sites=int(rng.integers(SITES_MIN, SITES_MAX + 1)),
                schedule=str(rng.choice(SCHEDULES)),
            )
        except ValueError:
            continue
    raise RuntimeError("Failed to sample a valid configuration after 20 attempts")


def configuration_space_summary() -> dict[str, Any]:
    """Return a summary of the configuration space (for documentation / debugging)."""
    n_valid_cell_payload = sum(
        1 for ct in CELL_TYPES for p in PAYLOADS
        if (ct, p) not in {("HB1_F3_line", "CAR"), ("MSC", "CRAd")}
    )
    return {
        "n_cell_types": len(CELL_TYPES),
        "n_payloads": len(PAYLOADS),
        "n_valid_cell_payload_combos": n_valid_cell_payload,
        "dose_range": [DOSE_MIN, DOSE_MAX],
        "n_injection_sites_range": [SITES_MIN, SITES_MAX],
        "n_schedules": len(SCHEDULES),
        "n_discrete_combos": n_valid_cell_payload * SITES_MAX * len(SCHEDULES),
        "continuous_dimensions": 1,  # dose_cells
        "constraints": [
            "HB1_F3_line × CAR incompatible",
            "MSC × CRAd incompatible",
            "dose within delivery-vehicle capacity",
            "injection sites along resection cavity (anatomical constraint)",
        ],
    }
