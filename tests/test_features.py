"""Tests for L3 in synthetic mode and for the layout the L5 to L7 stubs consume.

The contract rules themselves (ownership, reasons, views, real mode) are in
tests/test_feature_store_contract.py.
"""
from __future__ import annotations

import pandas as pd
import pytest

from gliostrata.features.contract import FeatureStoreError, LabelInFeatureMatrixError
from gliostrata.features.store import (
    FEATURE_COLUMNS,
    build_synthetic_store,
    feature_matrix,
    stub_layer_frame,
)


@pytest.fixture
def synthetic_store(medium_cohort):
    from gliostrata.data.harmonize import harmonize
    from gliostrata.profiling.states import deconvolve_states
    from gliostrata.profiling.subtype import classify_subtypes
    from gliostrata.profiling.tme import deconvolve_tme

    cohort = harmonize(medium_cohort)
    return build_synthetic_store(cohort, classify_subtypes(cohort),
                                 deconvolve_states(cohort), deconvolve_tme(cohort))


class TestSyntheticStore:
    def test_every_patient_one_row(self, synthetic_store):
        assert len(synthetic_store.features) == 30
        assert synthetic_store.features["patient_id"].is_unique

    def test_labelled_synthetic(self, synthetic_store):
        assert synthetic_store.mode == "synthetic"
        for rec in synthetic_store.column_records():
            assert rec["source_layer"].startswith("synthetic:")

    def test_outcomes_are_in_the_label_table_only(self, synthetic_store):
        assert "vital_status" not in synthetic_store.features.columns
        assert "vital_status" in synthetic_store.labels.columns


class TestStubLayerFrame:
    """The one-hot layout of the L5 to L7 stubs (conftest's feature_store)."""

    def test_shape(self, feature_store):
        assert len(feature_store) == 30

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
        assert set(feature_store["sex_male"].unique()).issubset({0, 1})
        assert set(feature_store["extent_of_resection_STR"].unique()).issubset({0, 1})
        assert set(feature_store["extent_of_resection_GTR"].unique()).issubset({0, 1})

    def test_no_outcome_label_in_the_frame(self, feature_store):
        assert "overall_survival_months" not in feature_store.columns
        assert "vital_status" not in feature_store.columns

    def test_attrs_set(self, feature_store):
        assert "feature_columns" in feature_store.attrs
        assert "version" in feature_store.attrs

    def test_refuses_a_store_with_an_absent_value(self, synthetic_store):
        features = synthetic_store.features.copy()
        features.loc[0, "extent_of_resection"] = pd.NA
        with pytest.raises(FeatureStoreError, match="absent values"):
            stub_layer_frame(features)


class TestFeatureMatrix:
    def test_complete_case_matrix(self, feature_store):
        matrix = feature_matrix(feature_store, FEATURE_COLUMNS)
        assert matrix.shape == (len(feature_store), len(FEATURE_COLUMNS))
        assert "patient_id" not in matrix.columns

    def test_a_label_is_refused(self, synthetic_store):
        with pytest.raises(LabelInFeatureMatrixError):
            feature_matrix(synthetic_store, ["age", "overall_survival_months"])

    def test_a_gap_is_refused(self, synthetic_store):
        features = synthetic_store.features.copy()
        features.loc[0, "age"] = pd.NA
        with pytest.raises(FeatureStoreError, match="missing values"):
            feature_matrix(features, ["age"])
