"""L6b — Gaussian-process Bayesian optimization over the therapeutic configuration space.

This module is built ALONGSIDE the random-search optimizer in
``gliostrata.optimization.bayesian`` (which stays the baseline the protocol's
T4 test requires an optimizer to beat). Nothing here replaces it; integration
into the pipeline is a separate, reviewed step.

Backend choice
--------------
BoTorch/GPyTorch are declared only in the optional "heavy" extra and are NOT
installed in this environment (``import botorch`` and ``import gpytorch`` both
fail; ``torch`` is present but is not sufficient on its own). We therefore use
scikit-learn's :class:`~sklearn.gaussian_process.GaussianProcessRegressor`
(sklearn is a core dependency) with a Matern(5/2) kernel and an
expected-improvement acquisition function. The surrogate is isolated behind the
:class:`Surrogate` protocol and :class:`SklearnGPSurrogate`, so a
``BoTorchSurrogate`` can be dropped in behind the same interface later without
touching the optimizer loop, the encoder, or the benchmark harness.

Objective pluggability (important caveat)
-----------------------------------------
The production objective would be computed from L4 (simulation) and L5
(survival), **both of which are currently stubs**. The optimizer therefore takes
the objective as a plain callable ``Callable[[Configuration], float]`` (to be
MAXIMIZED), and is validated in the test suite against synthetic objectives
whose global optimum is known analytically. No number produced by running this
optimizer on the current stub objective may be reported as a scientific finding.

Mixed search space handling
---------------------------
``cell_type``, ``payload`` and ``schedule`` are categorical and are one-hot
encoded for the GP kernel. ``dose_cells`` spans 1e6..5e7 (~1.7 orders of
magnitude) and is encoded as *normalized log10 dose*. ``n_injection_sites`` is
an integer; it is encoded as a normalized continuous value and rounded on
decode (see :meth:`ConfigEncoder.decode_sites` for an honest note on the bias
this introduces).

Validity is delegated entirely to
:func:`gliostrata.optimization.config_space.validate_configuration` — the
constraints are not reimplemented here. Candidate configurations are generated
so that only valid (cell_type, payload) combinations are ever enumerated, and
every proposal is additionally re-validated before it is returned.
"""
from __future__ import annotations

import math
import time
import warnings
from dataclasses import dataclass, field
from typing import Callable, Protocol, Sequence

import numpy as np

from gliostrata.optimization.config_space import (
    CELL_TYPES,
    DOSE_MAX,
    DOSE_MIN,
    PAYLOADS,
    SCHEDULES,
    SITES_MAX,
    SITES_MIN,
    Configuration,
    validate_configuration,
)

__all__ = [
    "SearchSpace",
    "ConfigEncoder",
    "Surrogate",
    "SklearnGPSurrogate",
    "expected_improvement",
    "OptimizationResult",
    "gp_ei_optimize",
    "random_search_optimize",
    "simple_regret",
    "benchmark_regret",
    "ObjectiveFn",
]

ObjectiveFn = Callable[[Configuration], float]


def _raw_config(cell_type: str, payload: str, dose: int, sites: int,
                schedule: str) -> Configuration:
    """Build a Configuration. ``Configuration.__post_init__`` validates, so this
    raises ValueError for an invalid combination — we never bypass the validator."""
    return Configuration(
        cell_type=cell_type,
        payload=payload,
        dose_cells=int(dose),
        n_injection_sites=int(sites),
        schedule=schedule,
    )


# --------------------------------------------------------------------------
# Search space
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class SearchSpace:
    """The (optionally restricted) region of the configuration space to search.

    Defaults to the full space defined in ``config_space``. Restricting it is
    useful for tests (e.g. a single-candidate space) and for clinically
    constrained per-patient searches. Validity constraints are *always* applied
    on top of this via ``config_space.validate_configuration``.
    """

    cell_types: tuple[str, ...] = tuple(CELL_TYPES)
    payloads: tuple[str, ...] = tuple(PAYLOADS)
    schedules: tuple[str, ...] = tuple(SCHEDULES)
    dose_bounds: tuple[int, int] = (DOSE_MIN, DOSE_MAX)
    sites_bounds: tuple[int, int] = (SITES_MIN, SITES_MAX)

    def __post_init__(self) -> None:
        if not self.cell_types or not self.payloads or not self.schedules:
            raise ValueError("SearchSpace categorical dimensions must be non-empty")
        for name, allowed, given in (
            ("cell_types", CELL_TYPES, self.cell_types),
            ("payloads", PAYLOADS, self.payloads),
            ("schedules", SCHEDULES, self.schedules),
        ):
            bad = [v for v in given if v not in allowed]
            if bad:
                raise ValueError(f"SearchSpace.{name} contains unknown values: {bad}")
        lo, hi = self.dose_bounds
        if not (DOSE_MIN <= lo <= hi <= DOSE_MAX):
            raise ValueError(f"dose_bounds {self.dose_bounds} outside [{DOSE_MIN}, {DOSE_MAX}]")
        slo, shi = self.sites_bounds
        if not (SITES_MIN <= slo <= shi <= SITES_MAX):
            raise ValueError(f"sites_bounds {self.sites_bounds} outside "
                             f"[{SITES_MIN}, {SITES_MAX}]")
        if not self.valid_cell_payload_combos():
            raise ValueError("SearchSpace contains no valid (cell_type, payload) combination")

    def valid_cell_payload_combos(self) -> list[tuple[str, str]]:
        """Enumerate (cell_type, payload) pairs that pass the config_space validator."""
        combos: list[tuple[str, str]] = []
        dose = self.dose_bounds[0]
        sites = self.sites_bounds[0]
        schedule = self.schedules[0]
        for ct in self.cell_types:
            for p in self.payloads:
                try:
                    _raw_config(ct, p, dose, sites, schedule)
                except ValueError:
                    continue
                combos.append((ct, p))
        return combos

    def n_discrete_points(self) -> int:
        """Size of the space if dose were discrete at integer-cell resolution."""
        n_dose = self.dose_bounds[1] - self.dose_bounds[0] + 1
        n_sites = self.sites_bounds[1] - self.sites_bounds[0] + 1
        return len(self.valid_cell_payload_combos()) * n_dose * n_sites * len(self.schedules)


# --------------------------------------------------------------------------
# Encoding
# --------------------------------------------------------------------------
class ConfigEncoder:
    """Encode/decode configurations to a fixed-length numeric vector for the GP.

    Layout (all entries in [0, 1]):
        [one-hot cell_type | one-hot payload | log10-dose | sites | one-hot schedule]

    ``dose_cells`` is encoded as
    ``(log10(dose) - log10(dose_min)) / (log10(dose_max) - log10(dose_min))``.
    Log scaling is used because the dose range spans ~1.7 orders of magnitude;
    on a linear scale a stationary kernel would treat the difference between 1e6
    and 2e6 (a doubling) as negligible next to 4e7 vs 5e7 (a 25% change), which
    is biologically backwards.
    """

    def __init__(self, space: SearchSpace | None = None) -> None:
        self.space = space or SearchSpace()
        self.cell_types = list(self.space.cell_types)
        self.payloads = list(self.space.payloads)
        self.schedules = list(self.space.schedules)
        self.valid_combos = self.space.valid_cell_payload_combos()
        self.dim = len(self.cell_types) + len(self.payloads) + 2 + len(self.schedules)

    # -- continuous sub-encodings ------------------------------------------
    def encode_dose(self, dose_cells: int | float) -> float:
        lo = math.log10(self.space.dose_bounds[0])
        hi = math.log10(self.space.dose_bounds[1])
        if hi <= lo:
            return 0.0
        return float((math.log10(float(dose_cells)) - lo) / (hi - lo))

    def decode_dose(self, u: float) -> int:
        """Map u in [0, 1] back to an integer cell count, geometrically spaced."""
        lo = math.log10(self.space.dose_bounds[0])
        hi = math.log10(self.space.dose_bounds[1])
        u = float(np.clip(u, 0.0, 1.0))
        dose = int(round(10.0 ** (lo + u * (hi - lo))))
        return int(np.clip(dose, self.space.dose_bounds[0], self.space.dose_bounds[1]))

    def encode_sites(self, sites: int | float) -> float:
        lo, hi = self.space.sites_bounds
        if hi == lo:
            return 0.0
        return float((float(sites) - lo) / (hi - lo))

    def decode_sites(self, u: float) -> int:
        """Map u in [0, 1] back to an integer number of injection sites.

        HONEST NOTE ON ROUNDING BIAS: ``n_injection_sites`` is optimized on a
        continuous axis and rounded to the nearest integer. This is the standard
        "relax-then-round" treatment and it is biased in two known ways.

        (1) The GP is fit on the *rounded* integer values that were actually
        evaluated, so the surrogate never observes the fractional points the
        acquisition function could propose; an acquisition optimum landing
        between two integers can be rounded to a point with lower true
        acquisition value than some other integer.

        (2) With round-half-away semantics the interior integers each own an
        interval of width 1/(hi-lo) in u-space while the two endpoints own only
        half that width, so uniform sampling in u under-represents 1 and 8 by
        a factor of two.

        We accept this rather than hide it: the alternative (treating the 8 site
        counts as 8 one-hot categories) throws away the ordering information
        that makes the dimension learnable from few evaluations. The acquisition
        search here proposes candidates that are *already integral* (see
        :func:`_candidate_configs`), which removes effect (1) for the proposals
        themselves and leaves only the endpoint under-representation in (2).
        """
        lo, hi = self.space.sites_bounds
        u = float(np.clip(u, 0.0, 1.0))
        sites = int(round(lo + u * (hi - lo)))
        return int(np.clip(sites, lo, hi))

    # -- full vectors ------------------------------------------------------
    def encode(self, config: Configuration) -> np.ndarray:
        vec = np.zeros(self.dim, dtype=float)
        vec[self.cell_types.index(config.cell_type)] = 1.0
        off = len(self.cell_types)
        vec[off + self.payloads.index(config.payload)] = 1.0
        off += len(self.payloads)
        vec[off] = self.encode_dose(config.dose_cells)
        vec[off + 1] = self.encode_sites(config.n_injection_sites)
        off += 2
        vec[off + self.schedules.index(config.schedule)] = 1.0
        return vec

    def encode_many(self, configs: Sequence[Configuration]) -> np.ndarray:
        if not configs:
            return np.zeros((0, self.dim), dtype=float)
        return np.vstack([self.encode(c) for c in configs])

    def decode(self, vec: np.ndarray) -> Configuration:
        """Decode a vector to the nearest VALID configuration.

        Categoricals are decoded by argmax, except that (cell_type, payload) is
        decoded jointly: we pick the valid pair maximizing the sum of its two
        one-hot scores, so an incompatible pair can never be produced. Validity
        itself is defined by ``config_space.validate_configuration``.
        """
        vec = np.asarray(vec, dtype=float).ravel()
        n_cell, n_pay, n_sch = len(self.cell_types), len(self.payloads), len(self.schedules)
        cell_scores = vec[:n_cell]
        pay_scores = vec[n_cell:n_cell + n_pay]
        best_pair: tuple[str, str] | None = None
        best_score = -np.inf
        for ct, p in self.valid_combos:
            s = cell_scores[self.cell_types.index(ct)] + pay_scores[self.payloads.index(p)]
            if s > best_score:
                best_pair, best_score = (ct, p), s
        assert best_pair is not None  # SearchSpace guarantees >= 1 valid combo
        off = n_cell + n_pay
        dose = self.decode_dose(vec[off])
        sites = self.decode_sites(vec[off + 1])
        sched = self.schedules[int(np.argmax(vec[off + 2:off + 2 + n_sch]))]
        return _raw_config(best_pair[0], best_pair[1], dose, sites, sched)

    def key(self, config: Configuration) -> tuple:
        """A hashable identity for deduplicating evaluated configurations."""
        return (
            config.cell_type,
            config.payload,
            int(config.dose_cells),
            int(config.n_injection_sites),
            config.schedule,
        )


# --------------------------------------------------------------------------
# Surrogate interface + sklearn backend
# --------------------------------------------------------------------------
class Surrogate(Protocol):
    """Minimal surrogate interface, so a BoTorch backend can be swapped in."""

    def fit(self, X: np.ndarray, y: np.ndarray) -> "Surrogate":
        ...

    def predict(self, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return (mean, std) predictions."""
        ...


class SklearnGPSurrogate:
    """Gaussian-process surrogate backed by sklearn's GaussianProcessRegressor.

    Matern(nu=5/2) with per-dimension (ARD) length-scales over the encoded
    space, a constant scale factor and a white-noise term. ``normalize_y=True``
    so the objective's scale need not be tuned per problem.
    """

    def __init__(self, n_dims: int, *, seed: int = 0, n_restarts: int = 2,
                 noise_level: float = 1e-4) -> None:
        self.n_dims = int(n_dims)
        self.seed = int(seed)
        self.n_restarts = int(n_restarts)
        self.noise_level = float(noise_level)
        self._gp = None

    def _build(self):
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

        kernel = ConstantKernel(1.0, (1e-3, 1e3)) * Matern(
            length_scale=np.ones(self.n_dims),
            length_scale_bounds=(1e-2, 1e2),
            nu=2.5,
        ) + WhiteKernel(noise_level=self.noise_level, noise_level_bounds=(1e-10, 1e1))
        return GaussianProcessRegressor(
            kernel=kernel,
            normalize_y=True,
            n_restarts_optimizer=self.n_restarts,
            random_state=self.seed,
        )

    def fit(self, X: np.ndarray, y: np.ndarray) -> "SklearnGPSurrogate":
        gp = self._build()
        with warnings.catch_warnings():
            # A constant / near-constant objective drives kernel hyperparameters
            # to a bound; that is expected here, not an error.
            warnings.simplefilter("ignore")
            gp.fit(np.asarray(X, dtype=float), np.asarray(y, dtype=float))
        self._gp = gp
        return self

    def predict(self, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self._gp is None:
            raise RuntimeError("SklearnGPSurrogate.predict called before fit")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mu, std = self._gp.predict(np.asarray(X, dtype=float), return_std=True)
        return np.asarray(mu, dtype=float), np.asarray(std, dtype=float)


# --------------------------------------------------------------------------
# Acquisition
# --------------------------------------------------------------------------
def expected_improvement(mu: np.ndarray, sigma: np.ndarray, best_y: float,
                         xi: float = 0.01) -> np.ndarray:
    """Expected improvement for MAXIMIZATION.

    ``xi`` is the exploration/exploitation knob: the margin a candidate must
    beat before it counts as an improvement at all. ``xi=0`` is greedy
    improvement weighting; larger ``xi`` demands bigger gains and so pushes mass
    toward high-variance, unexplored regions. It is in units of the objective,
    so callers whose objective is not O(1) should scale it.

    Points with (numerically) zero predictive variance get EI = 0: there is no
    improvement to be had from re-measuring a known noiseless point. This is
    also what keeps a constant objective from producing NaNs.
    """
    from scipy.stats import norm

    mu = np.asarray(mu, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    imp = mu - float(best_y) - float(xi)
    ei = np.zeros_like(mu)
    ok = sigma > 1e-12
    z = np.zeros_like(mu)
    z[ok] = imp[ok] / sigma[ok]
    ei[ok] = imp[ok] * norm.cdf(z[ok]) + sigma[ok] * norm.pdf(z[ok])
    return np.maximum(ei, 0.0)


# --------------------------------------------------------------------------
# Candidate generation
# --------------------------------------------------------------------------
def _sobol_unit(n: int, d: int, seed: int) -> np.ndarray:
    """Scrambled Sobol points in [0,1)^d (deterministic in ``seed``)."""
    from scipy.stats import qmc

    with warnings.catch_warnings():
        # scipy warns when n is not a power of two; we accept the slightly
        # worse balance rather than forcing the budget to a power of two.
        warnings.simplefilter("ignore")
        eng = qmc.Sobol(d=d, scramble=True, seed=seed)
        return np.asarray(eng.random(max(int(n), 1)), dtype=float)


def _config_from_unit(u: np.ndarray, enc: ConfigEncoder) -> Configuration:
    """Map a 4-d unit-cube point to a valid configuration.

    Dimensions: [valid (cell_type, payload) combo, log-dose, sites, schedule].
    Only valid combos are indexed, so this cannot produce an invalid proposal.
    """
    combos = enc.valid_combos
    ci = min(int(float(u[0]) * len(combos)), len(combos) - 1)
    ct, p = combos[ci]
    dose = enc.decode_dose(float(u[1]))
    sites = enc.decode_sites(float(u[2]))
    si = min(int(float(u[3]) * len(enc.schedules)), len(enc.schedules) - 1)
    return _raw_config(ct, p, dose, sites, enc.schedules[si])


def _initial_design(n: int, enc: ConfigEncoder, seed: int) -> list[Configuration]:
    """Sobol (low-discrepancy, Latin-hypercube-like) initial design."""
    if n <= 0:
        return []
    pts = _sobol_unit(n, 4, seed)
    return [_config_from_unit(pts[i], enc) for i in range(n)]


def _perturb(config: Configuration, enc: ConfigEncoder, rng: np.random.Generator,
             scale: float = 0.12) -> Configuration:
    """Local move around ``config``: Gaussian step in log-dose and sites, with an
    occasional categorical flip. Always returns a valid configuration."""
    combos = enc.valid_combos
    ct, p = config.cell_type, config.payload
    if rng.random() < 0.3:
        ct, p = combos[int(rng.integers(len(combos)))]
    sched = config.schedule
    if rng.random() < 0.3:
        sched = enc.schedules[int(rng.integers(len(enc.schedules)))]
    u_dose = float(np.clip(enc.encode_dose(config.dose_cells) + rng.normal(0, scale), 0, 1))
    u_sites = float(np.clip(enc.encode_sites(config.n_injection_sites)
                            + rng.normal(0, scale), 0, 1))
    return _raw_config(ct, p, enc.decode_dose(u_dose), enc.decode_sites(u_sites), sched)


def _candidate_configs(enc: ConfigEncoder, rng: np.random.Generator, n: int,
                       incumbents: Sequence[Configuration]) -> list[Configuration]:
    """Acquisition-search candidate pool of already-integral, already-valid
    configurations: half global (uniform over the valid space), half local moves
    around the best configurations found so far."""
    n = max(int(n), 1)
    n_local = n // 2 if incumbents else 0
    n_global = n - n_local
    out: list[Configuration] = []
    for _ in range(n_global):
        out.append(_config_from_unit(rng.random(4), enc))
    for i in range(n_local):
        out.append(_perturb(incumbents[i % len(incumbents)], enc, rng))
    return out


# --------------------------------------------------------------------------
# Optimizer
# --------------------------------------------------------------------------
@dataclass
class OptimizationResult:
    """Outcome of one optimization run."""

    best_config: Configuration
    best_value: float
    configs: list[Configuration] = field(default_factory=list)
    values: list[float] = field(default_factory=list)
    n_initial: int = 0
    n_iterations: int = 0
    method: str = ""
    seed: int = 0
    wall_clock_s: float = 0.0

    @property
    def n_evaluations(self) -> int:
        return len(self.values)

    def best_trace(self) -> list[float]:
        """Running best-so-far value after each evaluation."""
        trace: list[float] = []
        best = -np.inf
        for v in self.values:
            best = max(best, v)
            trace.append(float(best))
        return trace

    def top_k(self, k: int = 5) -> list[tuple[Configuration, float]]:
        """The k distinct best (configuration, value) pairs, best first."""
        seen: set[tuple] = set()
        pairs = sorted(zip(self.configs, self.values), key=lambda cv: -cv[1])
        out: list[tuple[Configuration, float]] = []
        for cfg, val in pairs:
            ident = (cfg.cell_type, cfg.payload, int(cfg.dose_cells),
                     int(cfg.n_injection_sites), cfg.schedule)
            if ident in seen:
                continue
            seen.add(ident)
            out.append((cfg, float(val)))
            if len(out) >= k:
                break
        return out


def _split_budget(budget: int, n_initial: int) -> tuple[int, int]:
    if budget < 1:
        raise ValueError(f"budget must be >= 1, got {budget}")
    if n_initial < 1:
        raise ValueError(f"n_initial must be >= 1, got {n_initial}")
    # Budget smaller than the initial design: truncate the design and run zero
    # GP iterations rather than overspending the budget.
    n_init = min(n_initial, budget)
    return n_init, budget - n_init


def gp_ei_optimize(
    objective: ObjectiveFn,
    *,
    budget: int = 60,
    n_initial: int = 12,
    seed: int = 0,
    space: SearchSpace | None = None,
    xi: float = 0.01,
    n_candidates: int = 384,
    surrogate_factory: Callable[[int, int], Surrogate] | None = None,
) -> OptimizationResult:
    """Maximize ``objective`` over the configuration space with a GP + EI loop.

    Parameters
    ----------
    objective : callable
        ``Configuration -> float``, to be MAXIMIZED. Pluggable by design: the
        production L4/L5 objective is a stub today, so the optimizer is
        validated against objectives with known optima instead.
    budget : int
        Total number of objective evaluations, initial design INCLUDED. The
        protocol does not specify a budget for T4; it is explicit here.
    n_initial : int
        Size of the Sobol initial design before the GP takes over. If
        ``budget < n_initial`` the design is truncated to ``budget`` and no GP
        iterations run.
    seed : int
        Fixed seed; identical seeds give identical proposals and identical results.
    xi : float
        Expected-improvement exploration margin (see :func:`expected_improvement`).
    n_candidates : int
        Size of the per-iteration acquisition-search candidate pool.
    surrogate_factory : callable, optional
        ``(n_dims, seed) -> Surrogate``. Defaults to :class:`SklearnGPSurrogate`;
        this is the seam for a future BoTorch backend.
    """
    t0 = time.perf_counter()
    enc = ConfigEncoder(space)
    n_init, n_iters = _split_budget(int(budget), int(n_initial))
    rng = np.random.default_rng(seed)
    make_surrogate = surrogate_factory or (lambda d, s: SklearnGPSurrogate(d, seed=s))

    configs: list[Configuration] = []
    values: list[float] = []
    seen: set[tuple] = set()

    for cfg in _initial_design(n_init, enc, seed):
        validate_configuration(cfg)
        configs.append(cfg)
        values.append(float(objective(cfg)))
        seen.add(enc.key(cfg))

    for _ in range(n_iters):
        X = enc.encode_many(configs)
        y = np.asarray(values, dtype=float)
        best_y = float(np.max(y))
        incumbents = [c for c, _ in sorted(zip(configs, values), key=lambda cv: -cv[1])[:4]]
        cands = _candidate_configs(enc, rng, n_candidates, incumbents)
        fresh = [c for c in cands if enc.key(c) not in seen]
        pool = fresh if fresh else cands  # exhausted space: allow a repeat
        try:
            surrogate = make_surrogate(enc.dim, seed).fit(X, y)
            mu, sigma = surrogate.predict(enc.encode_many(pool))
            ei = expected_improvement(mu, sigma, best_y, xi=xi)
        except Exception:  # pragma: no cover - a surrogate failure must not abort
            ei = np.zeros(len(pool))
        if not np.any(ei > 0):
            # No signal (constant objective, or everything informative already
            # known): take the first pool candidate. Deterministic, and we do
            # not pretend it is an informed proposal.
            next_cfg = pool[0]
        else:
            next_cfg = pool[int(np.argmax(ei))]
        validate_configuration(next_cfg)
        configs.append(next_cfg)
        values.append(float(objective(next_cfg)))
        seen.add(enc.key(next_cfg))

    best_i = int(np.argmax(np.asarray(values, dtype=float)))
    return OptimizationResult(
        best_config=configs[best_i],
        best_value=float(values[best_i]),
        configs=configs,
        values=values,
        n_initial=n_init,
        n_iterations=n_iters,
        method="gp_ei_sklearn",
        seed=int(seed),
        wall_clock_s=time.perf_counter() - t0,
    )


def random_search_optimize(
    objective: ObjectiveFn,
    *,
    budget: int = 60,
    seed: int = 0,
    space: SearchSpace | None = None,
) -> OptimizationResult:
    """Uniform random-search baseline.

    Matches the strategy of ``bayesian.optimize_configurations`` (i.i.d. uniform
    draws of valid configurations) behind the same objective/budget interface,
    so regret is compared like for like.
    """
    t0 = time.perf_counter()
    if budget < 1:
        raise ValueError(f"budget must be >= 1, got {budget}")
    enc = ConfigEncoder(space)
    rng = np.random.default_rng(seed)
    configs: list[Configuration] = []
    values: list[float] = []
    for _ in range(int(budget)):
        cfg = _config_from_unit(rng.random(4), enc)
        validate_configuration(cfg)
        configs.append(cfg)
        values.append(float(objective(cfg)))
    best_i = int(np.argmax(np.asarray(values, dtype=float)))
    return OptimizationResult(
        best_config=configs[best_i],
        best_value=float(values[best_i]),
        configs=configs,
        values=values,
        n_initial=int(budget),
        n_iterations=0,
        method="random_search",
        seed=int(seed),
        wall_clock_s=time.perf_counter() - t0,
    )


# --------------------------------------------------------------------------
# Regret benchmark harness (the machinery protocol test T4 needs)
# --------------------------------------------------------------------------
def simple_regret(result: OptimizationResult, true_optimum: float) -> float:
    """Simple regret = true optimum - best value found, clipped at 0.

    Clipping only guards floating-point overshoot; a genuinely negative regret
    means ``true_optimum`` was declared wrong, and we warn in that case.
    """
    regret = float(true_optimum) - float(result.best_value)
    if regret < -1e-9:
        warnings.warn(
            f"best_value {result.best_value!r} exceeds declared true_optimum "
            f"{true_optimum!r}; the declared optimum is wrong.",
            RuntimeWarning,
            stacklevel=2,
        )
    return max(regret, 0.0)


def benchmark_regret(
    objective: ObjectiveFn,
    true_optimum: float,
    *,
    budget: int = 60,
    n_seeds: int = 15,
    n_initial: int = 12,
    seed0: int = 0,
    space: SearchSpace | None = None,
    xi: float = 0.01,
    n_candidates: int = 384,
    baseline: Callable[..., OptimizationResult] | None = None,
    method: Callable[..., OptimizationResult] | None = None,
) -> dict:
    """Compare GP-EI against random search by simple regret at a fixed budget.

    Both methods get the SAME per-seed evaluation budget and the same seed, and
    the comparison is a paired Wilcoxon signed-rank test over seeds
    (``scipy.stats.wilcoxon``, ``alternative="less"``: GP-EI regret lower than
    random-search regret).

    ``method`` / ``baseline`` override the two optimizers (used by the tests to
    feed in deliberately-identical or deliberately-worse methods); each is
    called as ``fn(objective, budget=..., seed=..., space=...)``.

    Returns a dict with per-seed regrets, means, percentage regret reduction,
    the Wilcoxon statistic and p-value, and timing. ``significant`` uses the
    protocol's alpha of 0.01; ``meets_t4_bar`` additionally requires the >= 25%
    regret reduction. Both flags describe THIS objective only — on a synthetic
    objective they do not constitute protocol test T4, which must be run on the
    real L4/L5 objective.
    """
    from scipy.stats import wilcoxon

    if n_seeds < 1:
        raise ValueError(f"n_seeds must be >= 1, got {n_seeds}")

    def _default_method(obj, *, budget, seed, space):
        return gp_ei_optimize(obj, budget=budget, n_initial=n_initial, seed=seed,
                              space=space, xi=xi, n_candidates=n_candidates)

    def _default_baseline(obj, *, budget, seed, space):
        return random_search_optimize(obj, budget=budget, seed=seed, space=space)

    run_method = method or _default_method
    run_baseline = baseline or _default_baseline

    method_regrets: list[float] = []
    baseline_regrets: list[float] = []
    method_times: list[float] = []
    baseline_times: list[float] = []
    for i in range(int(n_seeds)):
        s = int(seed0) + i
        m = run_method(objective, budget=int(budget), seed=s, space=space)
        b = run_baseline(objective, budget=int(budget), seed=s, space=space)
        method_regrets.append(simple_regret(m, true_optimum))
        baseline_regrets.append(simple_regret(b, true_optimum))
        method_times.append(m.wall_clock_s)
        baseline_times.append(b.wall_clock_s)

    m_arr = np.asarray(method_regrets, dtype=float)
    b_arr = np.asarray(baseline_regrets, dtype=float)
    diffs = m_arr - b_arr

    # scipy.stats.wilcoxon raises when every paired difference is zero (the two
    # methods are identical). That is not significance — it is no evidence.
    if n_seeds < 2 or np.allclose(diffs, 0.0):
        stat, pvalue = 0.0, 1.0
        test_note = ("all paired differences are zero (or fewer than 2 seeds): Wilcoxon "
                     "is undefined, reporting p=1.0 — no evidence of a difference")
    else:
        res = wilcoxon(m_arr, b_arr, alternative="less", zero_method="wilcox")
        stat, pvalue = float(res.statistic), float(res.pvalue)
        test_note = "paired Wilcoxon signed-rank, alternative='less' (method regret lower)"

    mean_m, mean_b = float(m_arr.mean()), float(b_arr.mean())
    if mean_b > 0:
        reduction = (mean_b - mean_m) / mean_b * 100.0
    else:
        # The baseline already attained the optimum on every seed: there is no
        # regret left to reduce. Not a 100% improvement, and not a failure.
        reduction = 0.0

    alpha = 0.01
    return {
        "budget": int(budget),
        "n_seeds": int(n_seeds),
        "n_initial": int(n_initial),
        "gp_regrets": method_regrets,
        "random_regrets": baseline_regrets,
        "mean_regret_gp_ei": mean_m,
        "mean_regret_random": mean_b,
        "median_regret_gp_ei": float(np.median(m_arr)),
        "median_regret_random": float(np.median(b_arr)),
        "regret_reduction_pct": float(reduction),
        "wilcoxon_statistic": float(stat),
        "wilcoxon_pvalue": float(pvalue),
        "alpha": alpha,
        "test_note": test_note,
        "significant": bool(pvalue < alpha),
        "required_reduction_pct": 25.0,
        "meets_t4_bar": bool(reduction >= 25.0 and pvalue < alpha),
        "mean_wall_clock_s_gp_ei": float(np.mean(method_times)),
        "mean_wall_clock_s_random": float(np.mean(baseline_times)),
        "caveat": (
            "Synthetic-objective benchmark only. This is NOT protocol test T4: "
            "T4 must be run on the real objective once L4 (simulation) and L5 "
            "(survival) are no longer stubs."
        ),
    }
