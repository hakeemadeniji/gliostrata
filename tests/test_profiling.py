"""Tests for L2 (subtype, states, tme)."""
from __future__ import annotations

import pandas as pd
import pytest

from gliostrata.profiling.states import STATES, deconvolve_states
from gliostrata.profiling.subtype import SUBTYPES, classify_subtypes
from gliostrata.profiling.tme import TME_FRACTIONS, deconvolve_tme


class TestSubtype:
    def test_returns_one_row_per_patient(self, small_cohort):
        out = classify_subtypes(small_cohort, seed=42)
        assert len(out) == len(small_cohort)
        assert "patient_id" in out.columns
        assert "subtype" in out.columns

    def test_class_probabilities_present(self, small_cohort):
        out = classify_subtypes(small_cohort, seed=42)
        for st in SUBTYPES:
            assert f"subtype_prob_{st}" in out.columns

    def test_class_probabilities_sum_to_one(self, small_cohort):
        out = classify_subtypes(small_cohort, seed=42)
        prob_cols = [f"subtype_prob_{st}" for st in SUBTYPES]
        prob_sum = out[prob_cols].sum(axis=1)
        assert (prob_sum > 0.98).all() and (prob_sum < 1.02).all()

    def test_argmax_matches_input_subtype(self, small_cohort):
        """Stub: classifier trusts the input subtype. Highest prob should match it."""
        out = classify_subtypes(small_cohort, seed=42)
        prob_cols = {f"subtype_prob_{st}": st for st in SUBTYPES}
        # Allow for noise in the stub: the argmax should match within ~80% of cases
        # (since we add gaussian noise to probabilities)
        n_match = sum(
            prob_cols[max(prob_cols, key=lambda c: row[c])] == row["subtype"]
            for _, row in out.iterrows()
        )
        assert n_match / len(out) >= 0.80
        # Just check all subtypes are in the expected set
        assert set(out["subtype"]).issubset(set(SUBTYPES))

    def test_missing_subtype_column_raises(self):
        bad = pd.DataFrame({"patient_id": ["a"], "overall_survival_months": [10]})
        with pytest.raises(ValueError, match="subtype"):
            classify_subtypes(bad)


class TestStates:
    def test_returns_one_row_per_patient(self, small_cohort):
        out = deconvolve_states(small_cohort, seed=42)
        assert len(out) == len(small_cohort)

    def test_state_columns_present(self, small_cohort):
        out = deconvolve_states(small_cohort, seed=42)
        for s in STATES:
            assert s in out.columns

    def test_state_fractions_sum_to_one(self, small_cohort):
        out = deconvolve_states(small_cohort, seed=42)
        state_sum = out[STATES].sum(axis=1)
        assert (state_sum > 0.98).all() and (state_sum < 1.02).all()


class TestTME:
    def test_returns_one_row_per_patient(self, small_cohort):
        out = deconvolve_tme(small_cohort, seed=42)
        assert len(out) == len(small_cohort)

    def test_tme_columns_present(self, small_cohort):
        out = deconvolve_tme(small_cohort, seed=42)
        for f in TME_FRACTIONS:
            assert f in out.columns
