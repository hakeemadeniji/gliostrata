"""Tests for L4 (simulation) and L5 (survival prediction)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gliostrata.prediction.survival import predict_survival, predict_survival_with_config
from gliostrata.simulation.nsc import CELL_MIGRATION_PARAMS, simulate_nsc_migration
from gliostrata.simulation.tumor import simulate_tumor_growth


class TestTumorGrowth:
    def test_returns_growth_curve(self, small_cohort):
        patient = small_cohort.iloc[0]
        curve = simulate_tumor_growth(patient, hours=168, seed=42)
        assert len(curve) == 169  # 0 to 168 inclusive
        assert set(curve.columns) == {"time_hours", "tumor_volume_mm3", "n_live_cells"}

    def test_volume_monotonic_increase(self, small_cohort):
        """Gompertzian growth should be monotonic increasing (in expectation)."""
        patient = small_cohort.iloc[0]
        curve = simulate_tumor_growth(patient, hours=240, seed=42)
        # Use a large smoothing window to handle noise
        smoothed = curve["tumor_volume_mm3"].rolling(window=20, center=True).mean()
        # First half should be smaller than second half (mostly)
        first_half = smoothed.iloc[20:80].mean()
        second_half = smoothed.iloc[160:220].mean()
        assert second_half > first_half

    def test_mesenchymal_grows_faster_than_proneural(self, small_cohort):
        """Subtype-dependent growth rate."""
        mes = small_cohort[small_cohort["subtype"] == "Mesenchymal"].iloc[0]
        pn = small_cohort[small_cohort["subtype"] == "Proneural"].iloc[0]
        mes_curve = simulate_tumor_growth(mes, hours=240, seed=42)
        pn_curve = simulate_tumor_growth(pn, hours=240, seed=42)
        # MES final volume > PN final volume (averaged over last 24h to handle noise)
        assert mes_curve["tumor_volume_mm3"].iloc[-24:].mean() > pn_curve["tumor_volume_mm3"].iloc[-24:].mean()


class TestNSCMigration:
    def test_returns_dict_with_required_keys(self, small_cohort):
        patient = small_cohort.iloc[0]
        result = simulate_nsc_migration(patient, "iPSC-NSC", seed=42)
        for key in ["cell_type", "mean_distance_um", "fraction_at_tumor", "half_life_h",
                    "n_cells_injected", "hours_simulated"]:
            assert key in result

    def test_unknown_cell_type_raises(self, small_cohort):
        patient = small_cohort.iloc[0]
        with pytest.raises(ValueError, match="Unknown cell_type"):
            simulate_nsc_migration(patient, "imaginary_cells", seed=42)

    def test_all_cell_types_supported(self, small_cohort):
        patient = small_cohort.iloc[0]
        for ct in CELL_MIGRATION_PARAMS:
            result = simulate_nsc_migration(patient, ct, seed=42)
            assert result["cell_type"] == ct

    def test_fraction_at_tumor_in_range(self, small_cohort):
        patient = small_cohort.iloc[0]
        result = simulate_nsc_migration(patient, "iPSC-NSC", seed=42)
        assert 0 <= result["fraction_at_tumor"] <= 1.0


class TestSurvival:
    def test_returns_one_row_per_patient(self, feature_store):
        pred = predict_survival(feature_store, seed=42)
        assert len(pred) == len(feature_store)

    def test_predicted_survival_positive(self, feature_store):
        pred = predict_survival(feature_store, seed=42)
        assert (pred["predicted_survival_months"] > 0).all()

    def test_idh_mutant_has_better_survival(self, feature_store):
        """IDH-mutant GBM has substantially better prognosis (literature)."""
        pred = predict_survival(feature_store, seed=42)
        df = feature_store.merge(pred, on="patient_id")
        # If we have IDH-mutated samples, their predicted survival should be longer
        idh_mut = df[df["idh_mutated"] == 1]
        idh_wt = df[df["idh_mutated"] == 0]
        if len(idh_mut) > 0 and len(idh_wt) > 0:
            assert idh_mut["predicted_survival_months"].median() > idh_wt["predicted_survival_months"].median()

    def test_predict_with_config(self, feature_store, survival_pred):
        """Smoke test for predict_survival_with_config."""
        # Build a tiny config DataFrame for 2 patients × 2 configs
        cfg = pd.DataFrame({
            "patient_id": [feature_store["patient_id"].iloc[0]] * 2,
            "config_id": ["cfg_0", "cfg_1"],
            "cell_type": ["iPSC-NSC", "MSC"],
            "payload": ["CRAd", "prodrug_enzyme"],
            "dose_cells": [2_000_000, 5_000_000],
            "n_injection_sites": [3, 5],
            "schedule": ["single", "fractionated_x2"],
            "rank": [0, 0],
        })
        out = predict_survival_with_config(feature_store, survival_pred, cfg)
        assert len(out) == len(cfg)
        assert "predicted_survival_with_config" in out.columns
        assert "survival_gain_vs_stupp_months" in out.columns
