"""I/O helpers for GlioStrata.

Parquet read/write, JSON read/write, path helpers.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def read_parquet(path: str | Path) -> pd.DataFrame:
    """Read a Parquet file as a DataFrame."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Parquet file not found: {p}")
    return pd.read_parquet(p)


def write_parquet(df: pd.DataFrame, path: str | Path) -> Path:
    """Write a DataFrame to Parquet, creating parent directories."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(p, index=False)
    return p


def read_json(path: str | Path) -> Any:
    """Read a JSON file."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"JSON file not found: {p}")
    with p.open() as f:
        return json.load(f)


def write_json(data: Any, path: str | Path, *, indent: int = 2) -> Path:
    """Write data to JSON, creating parent directories."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w") as f:
        json.dump(data, f, indent=indent, default=str)
    return p


def ensure_dir(path: str | Path) -> Path:
    """Create directory if it doesn't exist; return Path."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p
