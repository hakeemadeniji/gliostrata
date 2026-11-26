"""Tests for L6b GP-EI Bayesian optimization (``gliostrata.optimization.bayes_opt``).

The load-bearing tests here run against SYNTHETIC objectives whose global
optimum is known analytically. That is deliberate: the production objective is
computed from L4 (simulation) and L5 (survival), both of which are stubs today,
so the optimizer cannot be validated against them. Nothing in this file
exercises the stub objective, and nothing here is protocol test T4 — T4 must be
run on the real objective once L4/L5 are real.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from gliostrata.optimization.bayes_opt import (
    ConfigEncoder,
    SearchSpace,
    SklearnGPSurrogate,
    benchmark_regret,
    expected_improvement,
    gp_ei_optimize,
    random_search_optimize,
    simple_regret,
)
from gliostrata.optimization.config_space import (
    DOSE_MAX,
    DOSE_MIN,
    SITES_MAX,
    SITES_MIN,
    Configuration,
    sample_random_configuration,
    validate_configuration,
)

ENC = ConfigEncoder()

# Keep test budgets small: one GP fit per iteration is the expensive part, and
# it dominates; enlarging the acquisition candidate pool is comparatively cheap,
# so the pool stays generous while the evaluation budget stays modest.
TEST_BUDGET = 26
TEST_N_INITIAL = 10
TEST_N_CANDIDATES = 192
TEST_N_SEEDS = 5


# --------------------------------------------------------------------------
# Synthetic benchmark objectives with analytically known optima
# --------------------------------------------------------------------------
def obj_unimodal(config: Configuration) -> float:
    """Additive, single-peaked. Known optimum = 6.0.

    Maximized by cell_type=iPSC-NSC (+3), payload=combo (+2, a valid pairing),
    schedule=fractionated_x3 (+1), n_injection_sites=6 (zero penalty) and a dose
    whose normalized log position is 0.65. Because dose is integer-valued over
    1e6..5e7 the dose penalty can be driven arbitrarily close to 0, so the
    supremum 3+2+1 = 6.0 is the true optimum (attainable to ~1e-12).
    """
    u = ENC.encode_dose(config.dose_cells)
    value = 0.0
    if config.cell_type == "iPSC-NSC":
        value += 3.0
    if config.payload == "combo":
        value += 2.0
    if config.schedule == "fractionated_x3":
        value += 1.0
    value -= 4.0 * (u - 0.65) ** 2
    value -= 0.15 * (config.n_injection_sites - 6) ** 2
    return value


def _bump(u: float, center: float, height: float, width: float) -> float:
    """Compactly-supported quadratic (Epanechnikov) bump.

    Exactly ``height`` at ``u == center``, exactly 0 for ``|u - center| >= width``,
    and never greater than ``height`` anywhere. The compact support is the point:
    two such bumps placed further than ``width1 + width2`` apart do not interact
    AT ALL, so the maximum of their sum is exactly the taller height — no
    overlap term to account for.
    """
    t = (u - center) / width
    return height * max(0.0, 1.0 - t * t)


def obj_bimodal(config: Configuration) -> float:
    """Two dose peaks plus a categorical interaction. Known optimum = 6.0 EXACTLY.

    The dose axis has a decoy peak of height 2.0 centred at u=0.2 and the true
    peak of height 3.0 centred at u=0.8. The categorical part rewards the
    *interaction* (MSC, CAR) (+2.5) rather than either factor alone, plus
    schedule=single (+0.5), with n_injection_sites=3 incurring no penalty.

    WHY COMPACT BUMPS AND NOT GAUSSIANS: an earlier version of this objective
    used two Gaussians, ``2*exp(-20(u-0.2)^2) + 3*exp(-20(u-0.8)^2)``. Their
    tails overlap, so the true maximum of the sum is NOT 3.0 but ~3.0015 (and
    has no closed form — the stationarity condition is transcendental). The
    declared optimum of 6.0 was therefore wrong by 1.5e-3, which
    ``simple_regret``'s own guard caught at runtime. Since regret is the
    quantity protocol test T4 is measured in, "approximately right" is not good
    enough. Rather than declare a numerically-estimated constant, the objective
    is redesigned so its optimum is exact BY CONSTRUCTION: the bumps have
    disjoint support (|0.8-0.2| = 0.6 > 0.15 + 0.15), each is bounded above by
    its own height, and the categorical terms are exact constants. Hence
    sup = 2.5 + 0.5 + 3.0 = 6.0, attained at u = 0.8 (approached to within
    ~1e-16 because ``dose_cells`` is integral). The bimodality that makes this a
    real exploration test is unchanged. ``test_declared_optima_are_verified_by_
    brute_force`` re-verifies this independently.
    """
    u = ENC.encode_dose(config.dose_cells)
    value = 0.0
    if (config.cell_type, config.payload) == ("MSC", "CAR"):
        value += 2.5
    if config.schedule == "single":
        value += 0.5
    value += _bump(u, 0.2, 2.0, 0.15)
    value += _bump(u, 0.8, 3.0, 0.15)
    value -= 0.1 * abs(config.n_injection_sites - 3)
    return value


def obj_needle(config: Configuration) -> float:
    """Narrow ridge in log-dose at u=0.5, rewarding a single valid combo.
    Known optimum = 4.0 (HB1_F3_line + prodrug_enzyme, fractionated_x2, 8 sites)."""
    u = ENC.encode_dose(config.dose_cells)
    value = 2.0 * math.exp(-60.0 * (u - 0.5) ** 2)
    if (config.cell_type, config.payload) == ("HB1_F3_line", "prodrug_enzyme"):
        value += 1.0
    if config.schedule == "fractionated_x2":
        value += 0.5
    value += 0.5 * (config.n_injection_sites - SITES_MIN) / (SITES_MAX - SITES_MIN)
    return value


def obj_constant(config: Configuration) -> float:
    """No signal at all. Known optimum = 2.5 everywhere."""
    return 2.5


OBJECTIVES = [
    ("unimodal", obj_unimodal, 6.0),
    ("bimodal", obj_bimodal, 6.0),
    ("needle", obj_needle, 4.0),
]


# --------------------------------------------------------------------------
# Independent verification of every declared optimum
# --------------------------------------------------------------------------
def _brute_force_max(objective, *, n_coarse: int = 201, n_refine: int = 1001,
                     refine_rounds: int = 3, n_top: int = 3):
    """Deterministic exhaustive-then-refined search for an objective's true max.

    This is deliberately INDEPENDENT of the optimizer and of the analytic
    reasoning in each objective's docstring: it enumerates every discrete
    configuration (all valid (cell_type, payload) pairs x all schedules x all
    8 site counts) against a log-spaced grid on the dose axis, then refines the
    best few discrete combinations by repeatedly re-gridding a shrinking window
    around the incumbent dose.

    Exhaustive enumeration of the discrete dimensions means the only
    approximation is on the (near-continuous) dose axis. After
    ``refine_rounds`` rounds the dose window is ~(1/n_coarse) * (2/n_refine)**
    ``refine_rounds`` wide, i.e. below 1e-11 in normalized-log-dose units for
    the defaults, so for these objectives (curvature O(10)) the residual error
    is far below the 1e-6 tolerance asserted against.

    Returns ``(best_value, best_config)``.
    """
    combos = SearchSpace().valid_cell_payload_combos()
    schedules = list(ENC.schedules)
    discrete = [(ct, p, sched, sites)
                for ct, p in combos
                for sched in schedules
                for sites in range(SITES_MIN, SITES_MAX + 1)]

    def evaluate(ct, p, sched, sites, u):
        cfg = Configuration(ct, p, ENC.decode_dose(u), sites, sched)
        return objective(cfg), cfg

    coarse_grid = np.linspace(0.0, 1.0, n_coarse)
    scored = []
    for ct, p, sched, sites in discrete:
        best_v, best_u, best_cfg = -np.inf, 0.0, None
        for u in coarse_grid:
            v, cfg = evaluate(ct, p, sched, sites, u)
            if v > best_v:
                best_v, best_u, best_cfg = v, u, cfg
        scored.append((best_v, best_u, (ct, p, sched, sites), best_cfg))

    scored.sort(key=lambda t: -t[0])
    overall_v, overall_cfg = scored[0][0], scored[0][3]
    half_width = 1.0 / (n_coarse - 1)
    for best_v, best_u, key, _cfg in scored[:n_top]:
        lo, hi = best_u - half_width, best_u + half_width
        v_here, cfg_here = best_v, None
        for _ in range(refine_rounds):
            grid = np.linspace(max(lo, 0.0), min(hi, 1.0), n_refine)
            u_here = best_u
            for u in grid:
                v, cfg = evaluate(*key, u)
                if v > v_here:
                    v_here, u_here, cfg_here = v, u, cfg
            span = (min(hi, 1.0) - max(lo, 0.0)) / (n_refine - 1)
            lo, hi = u_here - span, u_here + span
            best_u = u_here
        if v_here > overall_v:
            overall_v, overall_cfg = v_here, cfg_here or _cfg
    return overall_v, overall_cfg


@pytest.mark.parametrize("name,objective,optimum", OBJECTIVES + [("constant", obj_constant, 2.5)])
def test_declared_optima_are_verified_by_brute_force(name, objective, optimum):
    """Every declared optimum is re-derived by independent exhaustive search.

    One declared optimum (bimodal) was wrong, so all of them are verified the
    same way rather than assumed correct. Two separate assertions:
      1. the declared value is a genuine UPPER bound (no configuration beats
         it), which is what keeps simple regret non-negative; and
      2. it is TIGHT (the search gets within 1e-6 of it), which is what keeps
         regret meaningful rather than inflated by a loose bound.
    """
    found, cfg = _brute_force_max(objective)
    assert found <= optimum + 1e-9, (
        f"{name}: brute force found {found!r} > declared optimum {optimum!r} "
        f"at {cfg!r} — the declared optimum is not an upper bound"
    )
    assert optimum - found < 1e-6, (
        f"{name}: declared optimum {optimum!r} is loose — brute force only "
        f"reached {found!r} (gap {optimum - found:.3e})"
    )


def test_gaussian_overlap_regression_is_fixed():
    """Guard the specific bug: the two dose bumps must not interact at all.

    At the true peak (u=0.8) the decoy bump must contribute EXACTLY zero, so
    the maximum of the sum is exactly the taller height. With the previous
    Gaussian formulation the decoy contributed 1.5e-3 there, which is what made
    the declared optimum wrong.
    """
    assert _bump(0.8, 0.2, 2.0, 0.15) == 0.0
    assert _bump(0.2, 0.8, 3.0, 0.15) == 0.0
    assert _bump(0.8, 0.8, 3.0, 0.15) == 3.0
    # The bump is bounded above by its height everywhere, so no overshoot.
    for u in np.linspace(-0.2, 1.2, 141):
        assert _bump(u, 0.2, 2.0, 0.15) <= 2.0
        assert _bump(u, 0.8, 3.0, 0.15) <= 3.0
    # And the decoy really is a decoy: a strictly lower local optimum.
    assert _bump(0.2, 0.2, 2.0, 0.15) < _bump(0.8, 0.8, 3.0, 0.15)


def test_no_optimum_warning_is_emitted_for_the_declared_optima():
    """The suite must be free of the 'declared optimum is wrong' RuntimeWarning."""
    import warnings as _warnings

    for _name, objective, optimum in OBJECTIVES:
        with _warnings.catch_warnings():
            _warnings.simplefilter("error", RuntimeWarning)
            for seed in range(3):
                res = random_search_optimize(objective, budget=200, seed=seed)
                simple_regret(res, optimum)  # raises if it warns


# --------------------------------------------------------------------------
# Core claim: GP-EI gets closer to a known optimum than random search
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name,objective,optimum", OBJECTIVES)
def test_gp_ei_beats_random_search_on_known_optimum(name, objective, optimum):
    """Lower mean simple regret than random search at an identical budget,
    across multiple seeds, on an objective whose optimum we know."""
    res = benchmark_regret(
        objective,
        optimum,
        budget=TEST_BUDGET,
        n_seeds=TEST_N_SEEDS,
        n_initial=TEST_N_INITIAL,
        n_candidates=TEST_N_CANDIDATES,
    )
    assert res["mean_regret_gp_ei"] < res["mean_regret_random"], (
        f"{name}: GP-EI mean regret {res['mean_regret_gp_ei']:.4f} not better than "
        f"random search {res['mean_regret_random']:.4f}"
    )
    # Won on a majority of seeds, not just on the mean (one lucky seed).
    wins = sum(g < r for g, r in zip(res["gp_regrets"], res["random_regrets"]))
    assert wins > TEST_N_SEEDS / 2, f"{name}: GP-EI won only {wins}/{TEST_N_SEEDS} seeds"
    assert res["regret_reduction_pct"] > 0.0


@pytest.mark.parametrize("name,objective,optimum", OBJECTIVES)
def test_regret_is_nonnegative_and_bounded(name, objective, optimum):
    """Declared optima really are upper bounds (regret never negative)."""
    for seed in range(2):
        gp = gp_ei_optimize(objective, budget=14, n_initial=8, seed=seed,
                            n_candidates=96)
        rs = random_search_optimize(objective, budget=14, seed=seed)
        assert simple_regret(gp, optimum) >= 0.0
        assert simple_regret(rs, optimum) >= 0.0
        assert gp.best_value <= optimum + 1e-9


# --------------------------------------------------------------------------
# Determinism
# --------------------------------------------------------------------------
def test_gp_ei_is_deterministic_under_fixed_seed():
    """Identical seed => identical proposals, identical values, identical best."""
    kwargs = dict(budget=18, n_initial=8, seed=1234, n_candidates=96)
    a = gp_ei_optimize(obj_bimodal, **kwargs)
    b = gp_ei_optimize(obj_bimodal, **kwargs)
    assert [ENC.key(c) for c in a.configs] == [ENC.key(c) for c in b.configs]
    assert a.values == b.values
    assert a.best_config == b.best_config
    assert a.best_value == b.best_value


def test_different_seeds_give_different_trajectories():
    """Determinism must not be degenerate — the seed has to actually matter."""
    a = gp_ei_optimize(obj_bimodal, budget=16, n_initial=8, seed=0, n_candidates=64)
    b = gp_ei_optimize(obj_bimodal, budget=16, n_initial=8, seed=1, n_candidates=64)
    assert [ENC.key(c) for c in a.configs] != [ENC.key(c) for c in b.configs]


def test_random_search_is_deterministic_under_fixed_seed():
    a = random_search_optimize(obj_unimodal, budget=40, seed=99)
    b = random_search_optimize(obj_unimodal, budget=40, seed=99)
    assert [ENC.key(c) for c in a.configs] == [ENC.key(c) for c in b.configs]
    assert a.values == b.values


def test_benchmark_harness_is_deterministic():
    kwargs = dict(budget=16, n_seeds=3, n_initial=8, n_candidates=48)
    a = benchmark_regret(obj_unimodal, 6.0, **kwargs)
    b = benchmark_regret(obj_unimodal, 6.0, **kwargs)
    assert a["gp_regrets"] == b["gp_regrets"]
    assert a["random_regrets"] == b["random_regrets"]
    assert a["wilcoxon_pvalue"] == b["wilcoxon_pvalue"]


# --------------------------------------------------------------------------
# Validity: zero invalid proposals, ever
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name,objective,optimum", OBJECTIVES)
def test_every_proposal_is_valid_per_config_space(name, objective, optimum):
    """Every configuration the optimizer evaluates passes config_space's own validator."""
    for seed in range(2):
        for res in (
            gp_ei_optimize(objective, budget=14, n_initial=6, seed=seed, n_candidates=96),
            random_search_optimize(objective, budget=40, seed=seed),
        ):
            assert len(res.configs) == res.n_evaluations
            for cfg in res.configs:
                validate_configuration(cfg)  # raises ValueError if invalid


def test_decode_never_produces_an_incompatible_pairing():
    """Even a vector that maximally favours a forbidden pairing decodes to a valid one."""
    forbidden = {("HB1_F3_line", "CAR"), ("MSC", "CRAd")}
    for cell_type, payload in forbidden:
        vec = np.zeros(ENC.dim)
        vec[ENC.cell_types.index(cell_type)] = 10.0
        vec[len(ENC.cell_types) + ENC.payloads.index(payload)] = 10.0
        cfg = ENC.decode(vec)
        validate_configuration(cfg)
        assert (cfg.cell_type, cfg.payload) not in forbidden


def test_random_vectors_always_decode_to_valid_configs():
    rng = np.random.default_rng(0)
    for _ in range(300):
        cfg = ENC.decode(rng.random(ENC.dim))
        validate_configuration(cfg)


def test_search_space_enumerates_exactly_the_valid_combos():
    combos = SearchSpace().valid_cell_payload_combos()
    assert len(combos) == 10  # 3 cell types x 4 payloads - 2 incompatible
    assert ("HB1_F3_line", "CAR") not in combos
    assert ("MSC", "CRAd") not in combos


def test_search_space_rejects_a_space_with_no_valid_combo():
    with pytest.raises(ValueError, match="no valid"):
        SearchSpace(cell_types=("MSC",), payloads=("CRAd",))


@pytest.mark.parametrize("kwargs", [
    {"cell_types": ("not-a-cell-type",)},
    {"payloads": ("not-a-payload",)},
    {"schedules": ("not-a-schedule",)},
    {"dose_bounds": (DOSE_MIN - 1, DOSE_MAX)},
    {"dose_bounds": (DOSE_MAX, DOSE_MIN)},
    {"sites_bounds": (0, SITES_MAX)},
    {"sites_bounds": (SITES_MIN, SITES_MAX + 1)},
    {"cell_types": ()},
])
def test_search_space_rejects_malformed_bounds(kwargs):
    with pytest.raises(ValueError):
        SearchSpace(**kwargs)


# --------------------------------------------------------------------------
# Log scaling of dose
# --------------------------------------------------------------------------
def test_dose_encoding_is_log_scaled():
    assert ENC.encode_dose(DOSE_MIN) == pytest.approx(0.0)
    assert ENC.encode_dose(DOSE_MAX) == pytest.approx(1.0)
    # u=0.5 is the GEOMETRIC mean, not the arithmetic mean.
    geo = math.sqrt(DOSE_MIN * DOSE_MAX)
    assert ENC.encode_dose(geo) == pytest.approx(0.5, abs=1e-9)
    assert ENC.decode_dose(0.5) == pytest.approx(geo, rel=1e-6)
    # The arithmetic midpoint sits far above 0.5 on a log axis.
    arith = (DOSE_MIN + DOSE_MAX) / 2
    assert ENC.encode_dose(arith) > 0.8


def test_dose_encoding_is_monotone_and_round_trips():
    us = np.linspace(0.0, 1.0, 41)
    doses = [ENC.decode_dose(u) for u in us]
    assert doses == sorted(doses)
    assert all(DOSE_MIN <= d <= DOSE_MAX for d in doses)
    for d in doses:
        assert ENC.decode_dose(ENC.encode_dose(d)) == pytest.approx(d, rel=1e-6)
    # Equal steps in u are equal RATIOS in dose, not equal differences.
    r1 = ENC.decode_dose(0.4) / ENC.decode_dose(0.2)
    r2 = ENC.decode_dose(0.8) / ENC.decode_dose(0.6)
    assert r1 == pytest.approx(r2, rel=1e-3)


def test_dose_encoding_clips_out_of_range_u():
    assert ENC.decode_dose(-5.0) == DOSE_MIN
    assert ENC.decode_dose(5.0) == DOSE_MAX


# --------------------------------------------------------------------------
# Integer rounding of n_injection_sites
# --------------------------------------------------------------------------
def test_sites_decode_returns_integers_in_range():
    for u in np.linspace(-0.5, 1.5, 61):
        s = ENC.decode_sites(u)
        assert isinstance(s, int)
        assert SITES_MIN <= s <= SITES_MAX


def test_sites_decode_endpoints_and_rounding():
    assert ENC.decode_sites(0.0) == SITES_MIN
    assert ENC.decode_sites(1.0) == SITES_MAX
    # Nearest-integer rounding: u maps to 1 + 7u, so u=0.5 -> 4.5 -> 4 or 5.
    assert ENC.decode_sites(0.5) in (4, 5)
    # Just inside each integer's cell rounds to that integer.
    for target in range(SITES_MIN, SITES_MAX + 1):
        u = (target - SITES_MIN) / (SITES_MAX - SITES_MIN)
        assert ENC.decode_sites(u) == target
        assert ENC.decode_sites(u + 0.4 / (SITES_MAX - SITES_MIN)) == target
    assert ENC.encode_sites(ENC.decode_sites(0.3)) == pytest.approx(
        (ENC.decode_sites(0.3) - SITES_MIN) / (SITES_MAX - SITES_MIN))


def test_sites_rounding_endpoint_underrepresentation_is_as_documented():
    """The documented bias: uniform u under-represents the endpoints ~2x."""
    rng = np.random.default_rng(0)
    counts = np.zeros(SITES_MAX + 1, dtype=int)
    for u in rng.random(20_000):
        counts[ENC.decode_sites(u)] += 1
    interior = counts[SITES_MIN + 1:SITES_MAX].mean()
    assert counts[SITES_MIN] == pytest.approx(interior / 2, rel=0.15)
    assert counts[SITES_MAX] == pytest.approx(interior / 2, rel=0.15)


def test_all_proposed_sites_are_integral():
    res = gp_ei_optimize(obj_unimodal, budget=14, n_initial=6, seed=3, n_candidates=64)
    for cfg in res.configs:
        assert isinstance(cfg.n_injection_sites, int)
        assert isinstance(cfg.dose_cells, int)


# --------------------------------------------------------------------------
# Encoder round-trip
# --------------------------------------------------------------------------
def test_encode_decode_round_trip_on_random_configs():
    rng = np.random.default_rng(7)
    for _ in range(100):
        cfg = sample_random_configuration(rng)
        back = ENC.decode(ENC.encode(cfg))
        assert back.cell_type == cfg.cell_type
        assert back.payload == cfg.payload
        assert back.schedule == cfg.schedule
        assert back.n_injection_sites == cfg.n_injection_sites
        assert back.dose_cells == pytest.approx(cfg.dose_cells, rel=1e-6)


def test_encoder_dimension_layout():
    assert ENC.dim == 3 + 4 + 2 + 3
    vec = ENC.encode(Configuration("MSC", "CAR", 5_000_000, 4, "single"))
    assert vec[:3].sum() == 1.0          # one-hot cell type
    assert vec[3:7].sum() == 1.0         # one-hot payload
    assert vec[9:12].sum() == 1.0        # one-hot schedule
    assert 0.0 <= vec[7] <= 1.0 and 0.0 <= vec[8] <= 1.0


# --------------------------------------------------------------------------
# Expected improvement
# --------------------------------------------------------------------------
def test_expected_improvement_is_zero_at_zero_variance():
    ei = expected_improvement(np.array([5.0, 0.0]), np.array([0.0, 0.0]), best_y=1.0)
    assert np.all(ei == 0.0)


def test_expected_improvement_is_nonnegative_and_monotone():
    sigma = np.full(5, 0.5)
    mu = np.linspace(0.0, 2.0, 5)
    ei = expected_improvement(mu, sigma, best_y=1.0)
    assert np.all(ei >= 0.0)
    assert np.all(np.diff(ei) > 0)  # rises with predicted mean
    # Rises with uncertainty at a fixed (pessimistic) mean: exploration.
    ei_lo = expected_improvement(np.array([0.5]), np.array([0.1]), best_y=1.0)
    ei_hi = expected_improvement(np.array([0.5]), np.array([1.0]), best_y=1.0)
    assert ei_hi[0] > ei_lo[0]


def test_expected_improvement_xi_controls_exploration():
    """A larger xi demands a larger gain, so EI at a near-incumbent point shrinks."""
    mu, sigma = np.array([1.05]), np.array([0.2])
    greedy = expected_improvement(mu, sigma, best_y=1.0, xi=0.0)[0]
    cautious = expected_improvement(mu, sigma, best_y=1.0, xi=0.5)[0]
    assert greedy > cautious


def test_expected_improvement_handles_nan_free_constant_predictions():
    ei = expected_improvement(np.full(4, 2.5), np.zeros(4), best_y=2.5)
    assert np.all(np.isfinite(ei))
    assert np.all(ei == 0.0)


# --------------------------------------------------------------------------
# Wilcoxon harness behaviour
# --------------------------------------------------------------------------
def _crippled_baseline(objective, *, budget, seed, space):
    """A baseline that is clearly worse: it spends only 1 of its budget."""
    return random_search_optimize(objective, budget=1, seed=seed, space=space)


def _full_random(objective, *, budget, seed, space):
    return random_search_optimize(objective, budget=budget, seed=seed, space=space)


def test_wilcoxon_detects_a_clearly_better_method():
    res = benchmark_regret(
        obj_unimodal, 6.0, budget=40, n_seeds=12,
        method=_full_random, baseline=_crippled_baseline,
    )
    assert res["mean_regret_gp_ei"] < res["mean_regret_random"]
    assert res["wilcoxon_pvalue"] < 0.01
    assert res["significant"] is True
    assert res["wilcoxon_statistic"] >= 0.0


def test_wilcoxon_claims_nothing_when_methods_are_identical():
    """Identical methods => all paired differences zero => no claim of significance."""
    res = benchmark_regret(
        obj_unimodal, 6.0, budget=30, n_seeds=10,
        method=_full_random, baseline=_full_random,
    )
    assert res["gp_regrets"] == res["random_regrets"]
    assert res["wilcoxon_pvalue"] == 1.0
    assert res["significant"] is False
    assert res["meets_t4_bar"] is False
    assert res["regret_reduction_pct"] == pytest.approx(0.0)
    assert "zero" in res["test_note"]


def test_benchmark_reports_the_t4_bar_explicitly():
    res = benchmark_regret(obj_unimodal, 6.0, budget=16, n_seeds=3, n_initial=8,
                           n_candidates=48)
    assert res["required_reduction_pct"] == 25.0
    assert res["alpha"] == 0.01
    assert res["budget"] == 16 and res["n_seeds"] == 3
    assert isinstance(res["meets_t4_bar"], bool)
    # The harness must carry the caveat that this is not protocol test T4.
    assert "NOT protocol test T4" in res["caveat"]
    assert res["meets_t4_bar"] == (res["regret_reduction_pct"] >= 25.0
                                   and res["wilcoxon_pvalue"] < 0.01)


def test_benchmark_rejects_zero_seeds():
    with pytest.raises(ValueError, match="n_seeds"):
        benchmark_regret(obj_unimodal, 6.0, n_seeds=0)


def test_simple_regret_clips_and_warns_on_a_wrong_declared_optimum():
    res = random_search_optimize(obj_unimodal, budget=10, seed=0)
    with pytest.warns(RuntimeWarning, match="declared optimum is wrong"):
        regret = simple_regret(res, -100.0)
    assert regret == 0.0


# --------------------------------------------------------------------------
# Edge cases
# --------------------------------------------------------------------------
def test_budget_smaller_than_initial_design():
    """The initial design is truncated to the budget; no GP iterations run."""
    res = gp_ei_optimize(obj_unimodal, budget=4, n_initial=20, seed=0)
    assert res.n_evaluations == 4
    assert res.n_initial == 4
    assert res.n_iterations == 0
    validate_configuration(res.best_config)


def test_budget_of_one():
    res = gp_ei_optimize(obj_unimodal, budget=1, n_initial=10, seed=0)
    assert res.n_evaluations == 1
    assert res.best_value == pytest.approx(obj_unimodal(res.best_config))


@pytest.mark.parametrize("kwargs", [
    {"budget": 0}, {"budget": -3}, {"n_initial": 0}, {"n_initial": -1},
])
def test_invalid_budget_parameters_raise(kwargs):
    with pytest.raises(ValueError):
        gp_ei_optimize(obj_unimodal, **kwargs)


def test_random_search_rejects_zero_budget():
    with pytest.raises(ValueError, match="budget"):
        random_search_optimize(obj_unimodal, budget=0)


def test_constant_objective_does_not_crash_or_claim_improvement():
    """No signal: the optimizer must run, stay valid, and claim nothing."""
    res = gp_ei_optimize(obj_constant, budget=16, n_initial=8, seed=0, n_candidates=64)
    assert res.n_evaluations == 16
    assert res.best_value == pytest.approx(2.5)
    assert all(v == pytest.approx(2.5) for v in res.values)
    for cfg in res.configs:
        validate_configuration(cfg)
    assert simple_regret(res, 2.5) == pytest.approx(0.0)

    bench = benchmark_regret(obj_constant, 2.5, budget=14, n_seeds=4, n_initial=6,
                             n_candidates=32)
    assert bench["mean_regret_gp_ei"] == pytest.approx(0.0)
    assert bench["mean_regret_random"] == pytest.approx(0.0)
    assert bench["regret_reduction_pct"] == pytest.approx(0.0)
    assert bench["significant"] is False
    assert bench["meets_t4_bar"] is False


def test_constant_objective_trace_is_flat():
    res = gp_ei_optimize(obj_constant, budget=12, n_initial=5, seed=2, n_candidates=32)
    assert res.best_trace() == [pytest.approx(2.5)] * 12


def test_single_candidate_space():
    """A space containing exactly one configuration: no crash, no invalid proposal."""
    space = SearchSpace(
        cell_types=("MSC",), payloads=("CAR",), schedules=("single",),
        dose_bounds=(5_000_000, 5_000_000), sites_bounds=(3, 3),
    )
    assert space.n_discrete_points() == 1
    res = gp_ei_optimize(obj_bimodal, budget=12, n_initial=4, seed=0, space=space,
                         n_candidates=16)
    assert res.n_evaluations == 12
    keys = {ENC.key(c) for c in res.configs}
    assert len(keys) == 1
    only = res.configs[0]
    assert (only.cell_type, only.payload, only.dose_cells,
            only.n_injection_sites, only.schedule) == ("MSC", "CAR", 5_000_000, 3, "single")
    assert res.best_value == pytest.approx(obj_bimodal(only))
    for cfg in res.configs:
        validate_configuration(cfg)
    assert len(res.top_k(5)) == 1


def test_restricted_space_is_respected():
    space = SearchSpace(cell_types=("iPSC-NSC",), schedules=("fractionated_x2",),
                        dose_bounds=(2_000_000, 4_000_000), sites_bounds=(2, 5))
    res = gp_ei_optimize(obj_unimodal, budget=14, n_initial=6, seed=0, space=space,
                         n_candidates=32)
    for cfg in res.configs:
        validate_configuration(cfg)
        assert cfg.cell_type == "iPSC-NSC"
        assert cfg.schedule == "fractionated_x2"
        assert 2_000_000 <= cfg.dose_cells <= 4_000_000
        assert 2 <= cfg.n_injection_sites <= 5


# --------------------------------------------------------------------------
# Result object
# --------------------------------------------------------------------------
def test_best_trace_is_monotone_non_decreasing():
    res = gp_ei_optimize(obj_bimodal, budget=16, n_initial=8, seed=0, n_candidates=64)
    trace = res.best_trace()
    assert len(trace) == res.n_evaluations
    assert all(b >= a for a, b in zip(trace, trace[1:]))
    assert trace[-1] == pytest.approx(res.best_value)


def test_top_k_is_distinct_and_descending():
    res = gp_ei_optimize(obj_bimodal, budget=18, n_initial=10, seed=0, n_candidates=64)
    top = res.top_k(5)
    assert len(top) <= 5
    values = [v for _, v in top]
    assert values == sorted(values, reverse=True)
    keys = [ENC.key(c) for c, _ in top]
    assert len(set(keys)) == len(keys)
    assert top[0][1] == pytest.approx(res.best_value)


def test_result_metadata():
    res = gp_ei_optimize(obj_unimodal, budget=15, n_initial=5, seed=11, n_candidates=32)
    assert res.method == "gp_ei_sklearn"
    assert res.seed == 11
    assert res.n_initial == 5
    assert res.n_iterations == 10
    assert res.wall_clock_s > 0.0
    assert random_search_optimize(obj_unimodal, budget=5).method == "random_search"


# --------------------------------------------------------------------------
# Surrogate swap seam (for a future BoTorch backend)
# --------------------------------------------------------------------------
def test_surrogate_factory_is_honoured():
    """A custom surrogate is used, proving the BoTorch swap seam works."""
    calls = {"n": 0}

    class DummySurrogate:
        def __init__(self, n_dims, seed):
            self.n_dims = n_dims

        def fit(self, X, y):
            calls["n"] += 1
            self._n = len(y)
            return self

        def predict(self, X):
            # Constant mean, unit std: EI is positive everywhere, so the loop
            # must still produce valid proposals.
            return np.zeros(len(X)), np.ones(len(X))

    res = gp_ei_optimize(obj_unimodal, budget=12, n_initial=4, seed=0,
                         n_candidates=16,
                         surrogate_factory=lambda d, s: DummySurrogate(d, s))
    assert calls["n"] == 8
    assert res.n_evaluations == 12
    for cfg in res.configs:
        validate_configuration(cfg)


def test_sklearn_surrogate_predicts_mean_and_std():
    enc = ConfigEncoder()
    rng = np.random.default_rng(0)
    configs = [sample_random_configuration(rng) for _ in range(12)]
    X = enc.encode_many(configs)
    y = np.array([obj_unimodal(c) for c in configs])
    gp = SklearnGPSurrogate(enc.dim, seed=0).fit(X, y)
    mu, sd = gp.predict(X)
    assert mu.shape == (12,) and sd.shape == (12,)
    assert np.all(sd >= 0.0)
    # Interpolates its own training data closely.
    assert np.max(np.abs(mu - y)) < 1.0


def test_surrogate_predict_before_fit_raises():
    with pytest.raises(RuntimeError, match="before fit"):
        SklearnGPSurrogate(12).predict(np.zeros((1, 12)))


# --------------------------------------------------------------------------
# Personalization variance (the second half of the T4 criteria)
# --------------------------------------------------------------------------
def test_patient_specific_objectives_yield_distinct_best_configs():
    """Different per-patient objectives must not collapse to one configuration.

    This exercises the mechanism behind T4's "at least 3 distinct configs in
    top-5" criterion on synthetic patient objectives; it is NOT that test.
    """
    targets = [
        ("iPSC-NSC", "combo", "single", 2),
        ("MSC", "CAR", "fractionated_x2", 5),
        ("HB1_F3_line", "prodrug_enzyme", "fractionated_x3", 7),
        ("iPSC-NSC", "CRAd", "fractionated_x3", 8),
    ]

    def make_objective(cell_type, payload, schedule, sites):
        def objective(config):
            value = 0.0
            value += 2.0 if config.cell_type == cell_type else 0.0
            value += 2.0 if config.payload == payload else 0.0
            value += 1.0 if config.schedule == schedule else 0.0
            value -= 0.3 * abs(config.n_injection_sites - sites)
            return value
        return objective

    bests = []
    for i, t in enumerate(targets):
        res = gp_ei_optimize(make_objective(*t), budget=18, n_initial=8, seed=i,
                             n_candidates=96)
        validate_configuration(res.best_config)
        bests.append((res.best_config.cell_type, res.best_config.payload))
    assert len(set(bests)) >= 3, f"collapsed to {set(bests)}"
