"""Tests for scripts/build_tcga_feature_store.py, now a thin wrapper.

The script runs ``gliostrata build-feature-store --mode real``; the rules of
the store are tested in tests/test_feature_store_contract.py. Here: that the
wrapper adds nothing of its own, refuses a CGGA path, and that the old
workaround (dropping null placeholders by hand) is gone.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "build_tcga_feature_store.py"


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "gliostrata_build_tcga_feature_store_script", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_the_script_is_a_thin_wrapper_over_the_cli():
    text = SCRIPT_PATH.read_text(encoding="utf-8")
    assert '"build-feature-store", "--mode", "real"' in text
    for gone in ("def build(", "drop_placeholders", "build_feature_store(",
                 "classify_verhaak", "read_parquet"):
        assert gone not in text, f"the wrapper still carries its own {gone}"


def test_the_wrapper_forwards_to_real_mode_and_refuses_cgga(tmp_path, capsys):
    module = _load_module()
    with pytest.raises(SystemExit) as exc:
        module.main(["--expression", str(tmp_path / "absent" / "cgga" / "x.parquet"),
                     "--report", str(tmp_path / "r.json")])
    assert exc.value.code != 0
    err = capsys.readouterr()
    text = err.out + err.err
    assert "refusing" in text and "held-out cohort" in text
    assert not (tmp_path / "r.json").exists()


def test_the_wrapper_refuses_a_mode_flag():
    module = _load_module()
    with pytest.raises(SystemExit, match="always runs --mode real"):
        module.main(["--mode", "synthetic"])


def test_outcome_values_are_carried_and_counted_but_never_used():
    """Permuting every outcome leaves the feature table unchanged.

    The store build uses no outcome value: outcomes go to the label table and
    are counted as missing where they are missing, nothing else.
    """
    from gliostrata.features.store import build_feature_store
    from tests.test_feature_store_contract import _inputs

    a = build_feature_store(mode="real", **_inputs())
    kw = _inputs()
    kw["labels"]["overall_survival_months"] = kw["labels"][
        "overall_survival_months"].to_numpy()[::-1]
    kw["labels"]["vital_status"] = [0, 1, 0, None, 1, 0, 1, 0]
    b = build_feature_store(mode="real", **kw)
    assert a.features.equals(b.features)
    n_missing_outcome = (b.missing["column"] == "vital_status").sum()
    assert n_missing_outcome == 1
