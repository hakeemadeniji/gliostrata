"""Tests for L0 (ingest) and L1 (harmonize, qc)."""
from __future__ import annotations

import pandas as pd
import pytest

from gliostrata.data.harmonize import harmonize
from gliostrata.data.ingest import COHORT_COLUMNS, generate_synthetic_cohort
from gliostrata.data.qc import qc_report


class TestSyntheticCohort:
    def test_generates_correct_shape(self):
        df = generate_synthetic_cohort(n=20, seed=42)
        assert len(df) == 20
        assert list(df.columns) == COHORT_COLUMNS

    def test_seed_reproducible(self):
        df1 = generate_synthetic_cohort(n=10, seed=123)
        df2 = generate_synthetic_cohort(n=10, seed=123)
        pd.testing.assert_frame_equal(df1, df2)

    def test_different_seeds_different(self):
        df1 = generate_synthetic_cohort(n=10, seed=1)
        df2 = generate_synthetic_cohort(n=10, seed=2)
        assert not df1.equals(df2)

    def test_subtype_distribution(self):
        df = generate_synthetic_cohort(n=500, seed=42)
        subtype_counts = df["subtype"].value_counts(normalize=True)
        # All four subtypes should appear
        assert set(subtype_counts.index) == {"Classical", "Mesenchymal", "Proneural", "Neural"}
        # Mesenchymal should be the most common (~32%)
        assert subtype_counts["Mesenchymal"] > 0.25
        # Neural should be the least common (~11%)
        assert subtype_counts["Neural"] < 0.20

    def test_state_fractions_sum_to_one(self):
        df = generate_synthetic_cohort(n=20, seed=42)
        state_sum = df[["state_npc", "state_opc", "state_ac", "state_mes"]].sum(axis=1)
        assert (state_sum > 0.98).all() and (state_sum < 1.02).all()

    def test_survival_positive(self):
        df = generate_synthetic_cohort(n=20, seed=42)
        assert (df["overall_survival_months"] > 0).all()

    def test_patient_ids_unique(self):
        df = generate_synthetic_cohort(n=20, seed=42)
        assert df["patient_id"].is_unique


class TestHarmonize:
    def test_passes_through_with_metadata(self, small_cohort):
        out = harmonize(small_cohort)
        assert len(out) == len(small_cohort)
        assert "harmonization_version" in out.columns
        assert "batch_corrected" in out.columns
        # Stub: batch correction not applied yet
        assert (out["batch_corrected"] == False).all()

    def test_requires_columns(self):
        bad = pd.DataFrame({"foo": [1, 2]})
        with pytest.raises(ValueError, match="Missing required columns"):
            harmonize(bad)

    def test_rejects_non_dataframe(self):
        with pytest.raises(TypeError):
            harmonize([1, 2, 3])


class TestQC:
    def test_qc_report_shape(self, small_cohort):
        report = qc_report(small_cohort)
        assert len(report) == len(small_cohort)
        assert set(report.columns) >= {"patient_id", "qc_pass", "qc_flags", "estimated_purity"}

    def test_qc_pass_for_valid_cohort(self, small_cohort):
        report = qc_report(small_cohort)
        # Synthetic cohort has positive survival times → all should pass
        assert report["qc_pass"].all()
