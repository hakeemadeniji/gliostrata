"""Pins the committed TCGA-GBM fit of variant 1, the model T2 was registered against.

The five numbers the record quotes -- 594 patients passed, 592 modelled, 150 with
KPS set to the median 80, 22 with a missing treatment arm -- live in
tests/fixtures/tcga_fit_fingerprint.json, with the two hazard ratios that show
the missingness terms are themselves prognostic. One test pins the fixture and
runs everywhere. Another recomputes the three fits on the real cohort when it
is on disk and requires each design matrix to be byte-identical to the hash
recorded before the fold-wise imputation change.

Why hashes: variant 1's internal C-index is a quoted number, and variants 2 and
3 were scored once on CGGA from these exact matrices. A change to the design
code that altered one of them would silently move a quoted result.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gliostrata.prediction.cox import CoxPHModel, build_design_matrix

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = json.loads(
    (ROOT / "tests" / "fixtures" / "tcga_fit_fingerprint.json").read_text(encoding="utf-8")
)
REAL_COHORT = ROOT / FIXTURE["inputs"]["cohort"]["path"]
REAL_JOINED = ROOT / FIXTURE["inputs"]["joined"]["path"]

#: fixture section -> (feature tuple in the T2 script, impute, complete cases?)
VARIANTS = {
    "variant_1_baseline_committed": ("VARIANT_1", "median_indicator", False),
    "variant_2_transferable": ("VARIANT_2", "drop", True),
    "variant_3_molecular": ("VARIANT_3", "drop", True),
}


def test_the_fixture_pins_the_five_counts_and_the_two_hazard_ratios():
    """Runs without the data: the numbers the record quotes are these, exactly."""
    v1 = FIXTURE["variant_1_baseline_committed"]
    assert v1["n_passed"] == 594
    assert v1["n_modelled"] == 592
    assert v1["n_kps_imputed"] == 150
    assert v1["kps_imputation_value"] == 80.0
    assert v1["n_treatment_arm_missing"] == 22
    assert v1["kps_missing_indicator_hr"] == pytest.approx(1.57, abs=0.005)
    assert v1["treatment_arm_missing_level_hr"] == pytest.approx(1.20, abs=0.005)

    # The counts add up, and the ratios are the fitted ones.
    assert v1["n_passed"] - v1["n_dropped_missing_outcome"] == v1["n_modelled"]
    assert v1["n_kps_imputed"] < v1["n_modelled"]
    assert v1["kps_missing_indicator_hr"] == v1["hazard_ratios"]["kps_missing"]
    assert v1["treatment_arm_missing_level_hr"] == v1["hazard_ratios"]["treatment_arm_missing"]
    assert "kps_missing" in v1["feature_names"]
    assert "missing" in v1["categorical_levels"]["treatment_arm"]
    assert set(v1["n_missing_other_features"].values()) == {0}

    # Variants 2 and 3 impute nothing, anywhere.
    for name in ("variant_2_transferable", "variant_3_molecular"):
        assert FIXTURE[name]["impute"] == "drop"
        assert set(FIXTURE[name]["n_missing"].values()) == {0}
        assert FIXTURE[name]["n_fit"] == FIXTURE[name]["n_modelled"]
    assert FIXTURE["variant_2_transferable"]["n_fit"] == 570
    assert FIXTURE["variant_3_molecular"]["n_fit"] == 298


def _load_script():
    """Load scripts/validate_cox_cgga.py by path: definitions only, nothing runs."""
    path = ROOT / "scripts" / "validate_cox_cgga.py"
    spec = importlib.util.spec_from_file_location("gliostrata_validate_cox_cgga_script", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _sha_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _sha_array(arr) -> str:
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()


@pytest.mark.skipif(
    not (REAL_COHORT.is_file() and REAL_JOINED.is_file()),
    reason="real TCGA cohort not downloaded",
)
def test_the_real_tcga_fits_are_byte_identical_to_the_fingerprint():
    for key in ("cohort", "joined"):
        path = ROOT / FIXTURE["inputs"][key]["path"]
        assert _sha_file(path) == FIXTURE["inputs"][key]["sha256"], (
            f"{path.name} changed since the fixture was written: re-derive the "
            "counts deliberately and update the fixture, never loosen this test"
        )
    vcox = _load_script()
    tcga = vcox.load_tcga(
        pd.read_parquet(REAL_COHORT),
        pd.read_parquet(REAL_JOINED, columns=["patient_id", *vcox.CALL_COLUMNS]),
    )

    def fit_frame_for(features, complete_case):
        """The frame each variant is fitted on, as the T2 script prepares it."""
        if not complete_case:
            return tcga
        fit_frame, _ = vcox.complete_cases(tcga, tuple(features))
        for feat in features:
            if feat not in ("sex", "treatment_arm") and feat in fit_frame.columns:
                if pd.api.types.is_numeric_dtype(fit_frame[feat]) or feat in vcox.CALL_COLUMNS:
                    fit_frame[feat] = vcox._numeric(fit_frame, feat)
        return fit_frame

    for name, (attr, impute, complete_case) in VARIANTS.items():
        pinned = FIXTURE[name]
        features = tuple(getattr(vcox, attr))
        assert list(features) == pinned["features"]
        x, time, event, report = build_design_matrix(
            fit_frame_for(features, complete_case), features, impute=impute
        )
        assert list(report.feature_names) == pinned["feature_names"], name
        assert (
            _sha_array(x),
            _sha_array(time),
            _sha_array(event),
        ) == (pinned["X_sha256"], pinned["time_sha256"], pinned["event_sha256"]), (
            f"{name}: the design matrix is no longer the one on record"
        )
        assert list(x.shape) == pinned["shape"]
        assert report.n_input == pinned.get("n_passed", pinned.get("n_fit"))
        assert report.n_used == pinned["n_modelled"]
        assert report.dropped_missing_outcome == pinned.get("n_dropped_missing_outcome", 0)
        model = CoxPHModel().fit(x, time, event, feature_names=report.feature_names)
        assert model.n_events_ == pinned["n_events"], name
        ratios = model.summary().set_index("feature")["hazard_ratio"]
        assert {str(k): round(float(v), 6) for k, v in ratios.items()} == pinned["hazard_ratios"]
        if name == "variant_1_baseline_committed":
            assert report.n_imputed["kps"] == pinned["n_kps_imputed"] == 150
            assert report.imputed_values["kps"] == pinned["kps_imputation_value"]
            assert report.n_imputed["treatment_arm"] == pinned["n_treatment_arm_missing"]
            assert report.n_imputed["age"] == report.n_imputed["sex"] == 0
        else:
            assert set(int(v) for v in report.n_imputed.values()) == {0}
            assert report.imputed_values == {}
