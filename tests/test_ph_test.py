"""Proportional-hazards test. Expected values were cross-checked on 2026-11-27 against
lifelines.statistics.proportional_hazard_test(time_transform="rank"), which agreed to six
decimals on every term; lifelines is not a test dependency."""

import numpy as np
import pytest

from gliostrata.prediction.cox import CoxPHModel
from gliostrata.prediction.ph_test import proportional_hazards_test


def _data():
    rng = np.random.default_rng(1)
    n = 300
    x1 = rng.normal(size=n)
    x2 = rng.binomial(1, .4, n).astype(float)
    x3 = rng.normal(size=n)
    base = rng.exponential(10, n)
    t = base * np.exp(-(0.5 * x1 + 0.7 * x2))
    t = np.where(t > 5, 5 + (t - 5) * np.exp(0.8 * x3), t * np.exp(-0.8 * x3))
    c = rng.uniform(2, 40, n)
    return np.column_stack([x1, x2, x3]), np.minimum(t, c), (t <= c).astype(int)


def test_matches_the_reference_implementation():
    x_mat, time, event = _data()
    model = CoxPHModel(l2=0.0).fit(x_mat, time, event)
    res = proportional_hazards_test(x_mat, time, event, model.coef_, ["x1", "x2", "x3"])
    got = {r["term"]: r["chisq"] for r in res["terms"]}
    assert got["x1"] == pytest.approx(0.553233, abs=2e-5)
    assert got["x2"] == pytest.approx(10.720084, abs=2e-4)
    assert got["x3"] == pytest.approx(48.027691, abs=2e-4)
    assert res["global"]["df"] == 3 and res["n_deaths"] == int(event.sum())


def test_proportional_hazards_data_are_not_flagged():
    rng = np.random.default_rng(7)
    n = 400
    x = rng.normal(size=(n, 2))
    t = rng.exponential(np.exp(-(x @ np.array([0.6, -0.4])))) * 10
    c = rng.uniform(1, 30, n)
    time, event = np.minimum(t, c), (t <= c).astype(int)
    model = CoxPHModel(l2=0.0).fit(x, time, event)
    res = proportional_hazards_test(x, time, event, model.coef_)
    assert res["global"]["p"] > 0.01


def test_rejects_bad_input():
    x_mat, time, event = _data()
    with pytest.raises(ValueError):
        proportional_hazards_test(x_mat, time, event, np.zeros(2))
    with pytest.raises(ValueError):
        proportional_hazards_test(x_mat, time, np.zeros_like(event), np.zeros(3))
