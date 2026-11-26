"""Configuration loading for GlioStrata.

Loads YAML configs from the configs/ directory. Provides a typed Config object
with sensible defaults so the package works without any config files present.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class Paths:
    data_raw: str = "data/raw"
    data_interim: str = "data/interim"
    data_processed: str = "data/processed"
    data_external: str = "data/external"
    data_results: str = "data/results"
    models: str = "models"


@dataclass
class Thresholds:
    """Pre-registered success thresholds.

    These values MUST match what is committed to OSF before any model touches
    the frozen test partition. Do not change without filing a new OSF version.
    """
    t1_subtype_accuracy_cgga: float = 0.88
    t1_state_deconv_r: float = 0.75
    t2_survival_cindex_cgga: float = 0.65
    t2_beat_baseline_by: float = 0.03
    t2_brier_18m: float = 0.25
    t3_tumor_growth_rmse_pct: float = 15.0
    t3_nsc_migration_r: float = 0.65
    t4_optimizer_beats_random_pct: float = 25.0
    t4_personalization_variance: int = 3
    t5_rmst_gain_months: float = 3.0
    t5_success_fraction: float = 0.50
    t5_harm_fraction_max: float = 0.10


@dataclass
class Config:
    project_name: str = "gliostrata"
    version: str = "0.1.0"
    seed: int = 42
    paths: Paths = field(default_factory=Paths)
    thresholds: Thresholds = field(default_factory=Thresholds)
    raw: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Config":
        paths_data = data.get("paths", {})
        thresholds_data = data.get("thresholds", {})
        project = data.get("project", {})
        return cls(
            project_name=project.get("name", "gliostrata"),
            version=project.get("version", "0.1.0"),
            seed=project.get("seed", 42),
            paths=Paths(**{k: v for k, v in paths_data.items() if hasattr(Paths, k)}),
            thresholds=Thresholds(
                **{k: v for k, v in thresholds_data.items() if hasattr(Thresholds, k)}
            ),
            raw=data,
        )


def load_config(path: str | Path | None = None) -> Config:
    """Load a GlioStrata config from a YAML file.

    If path is None, returns default config. If path doesn't exist, raises FileNotFoundError.
    """
    if path is None:
        return Config()
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Config file not found: {p}")
    with p.open() as f:
        data = yaml.safe_load(f) or {}
    return Config.from_dict(data)


def default_config_path() -> Path:
    """Return the path to the default config, relative to the package root."""
    # Walk up from this file to find configs/default.yaml
    here = Path(__file__).resolve()
    for parent in [here.parent, *here.parents]:
        candidate = parent / "configs" / "default.yaml"
        if candidate.exists():
            return candidate
        # Also check one level up (src/gliostrata -> repo root)
        candidate2 = parent.parent / "configs" / "default.yaml"
        if candidate2.exists():
            return candidate2
    # Fallback — assume CWD has configs/default.yaml
    return Path("configs/default.yaml")
