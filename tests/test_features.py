"""Tests for L3 (feature store)."""
from __future__ import annotations

import pandas as pd
import pytest

from gliostrata.features.store import FEATURE_COLUMNS, build_feature_store, feature_matrix


class TestFeatureStore:
    def test_shape(self, feature_store):
        n_patients = 30  # from medium_cohort fixture
        assert len(feature_store) == n_patients

    def test_feature_columns_present(self, feature_store):
        for col in FEATURE_COLUMNS:
            assert col in feature_store.columns, f"Missing feature column: {col}"

    def test_patient_id_present(self, feature_store):
        assert "patient_id" in feature_store.columns
        assert feature_store["patient_id"].is_unique

    def test_no_missing_values_in_features(self, feature_store):
        for col in FEATURE_COLUMNS:
            assert not feature_store[col].isna().any(), f"NaN in feature column: {col}"

    def test_one_hot_encoded_categoricals(self, feature_store):
        # sex_male should be 0 or 1
        assert set(feature_store["sex_male"].unique()).issubset({0, 1})
        # extent_of_resection encoded
        assert set(feature_store["extent_of_resection_STR"].unique()).issubset({0, 1})
        assert set(feature_store["extent_of_resection_GTR"].unique()).issubset({0, 1})

    def test_attrs_set(self, feature_store):
        assert "feature_columns" in feature_store.attrs
        assert "version" in feature_store.attrs


class TestFeatureMatrix:
    def test_split(self, feature_store):
        X, y_meta = feature_matrix(feature_store)
        assert X.shape == (len(feature_store), len(FEATURE_COLUMNS))
        assert "patient_id" in y_meta.columns
        assert "overall_survival_months" in y_meta.columns
        assert "vital_status" in y_meta.columns

    def test_no_metadata_in_X(self, feature_store):
        X, _ = feature_matrix(feature_store)
        assert "patient_id" not in X.columns
        assert "overall_survival_months" not in X.columns
