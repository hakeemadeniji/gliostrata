"""Pytest fixtures shared across the GlioStrata test suite."""
from __future__ import annotations

import pytest


@pytest.fixture
def small_cohort():
    """A small (n=10) synthetic cohort for fast tests."""
    from gliostrata.data.ingest import generate_synthetic_cohort
    return generate_synthetic_cohort(n=10, seed=42)


@pytest.fixture
def medium_cohort():
    """A medium (n=30) synthetic cohort for optimizer/pipeline tests."""
    from gliostrata.data.ingest import generate_synthetic_cohort
    return generate_synthetic_cohort(n=30, seed=42)


@pytest.fixture
def feature_store(medium_cohort):
    """L3 feature store of the medium cohort (synthetic mode), in the one-hot
    layout the L5 to L7 stubs consume (gliostrata.features.store.stub_layer_frame)."""
    from gliostrata.data.harmonize import harmonize
    from gliostrata.features.store import synthetic_feature_frame
    from gliostrata.profiling.states import deconvolve_states
    from gliostrata.profiling.subtype import classify_subtypes
    from gliostrata.profiling.tme import deconvolve_tme

    cohort = harmonize(medium_cohort)
    subtypes = classify_subtypes(cohort)
    states = deconvolve_states(cohort)
    tme = deconvolve_tme(cohort)
    return synthetic_feature_frame(cohort, subtypes, states, tme)


@pytest.fixture
def survival_pred(feature_store):
    """L5 survival predictions for the medium cohort."""
    from gliostrata.prediction.survival import predict_survival
    return predict_survival(feature_store, seed=42)
