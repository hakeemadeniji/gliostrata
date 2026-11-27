"""Real, fitted Cox proportional-hazards survival model (L5).

Why this module exists
----------------------
`gliostrata.prediction.survival` is a STUB: it computes a log-hazard from
hard-coded coefficients plus Gaussian noise and is not fitted to any data.
This module is a genuine maximum-partial-likelihood Cox model. It is built
ALONGSIDE the stub so the existing pipeline keeps working; integration is a
separate, deliberate step.

Implementation choice
---------------------
`lifelines` and `scikit-survival` are declared as OPTIONAL "heavy" extras in
pyproject.toml and are NOT installed in this environment (verified: both raise
ModuleNotFoundError). numpy/scipy/scikit-learn/pandas ARE core dependencies.
Rather than adding a heavy dependency to the critical path of a protocol test,
Cox PH is implemented here directly by maximising the **Breslow** partial
likelihood with analytic gradient and Hessian via `scipy.optimize`. The Hessian
also gives asymptotic standard errors (hence hazard-ratio confidence intervals)
for free. If the project later wants Efron tie handling, penalised splines or
time-varying covariates, adding `lifelines>=0.27` would be the right call and
this module's API is deliberately close to it (`fit`, `predict_log_hazard`,
`predict_survival_function`).

Censoring
---------
`vital_status == 1` is a death (event). `vital_status == 0` means the patient
was ALIVE at last follow-up: a **right-censored** observation, not a death at
that time. Censored observations contribute to risk sets but never to the
event terms of the partial likelihood. Treating them as deaths biases
coefficients toward zero (attenuation) -- `tests/test_cox.py` demonstrates
that bias numerically.

Missingness
-----------
Handled explicitly, never silently. The default strategy is
``median_indicator``: numeric features are imputed with the TRAINING-fold
median and a companion binary ``<feature>_missing`` column is added, so
"KPS was never recorded" is itself available to the model as signal.
Categorical features get an explicit ``missing`` level. Every imputation is
recorded in the returned :class:`DesignReport`.

"Training-fold" is enforced, not assumed. A design is fitted on training rows
(:func:`fit_design_plan`: medians, indicator columns, categorical levels,
reference levels) and applied to whichever rows are being encoded
(:func:`apply_design_plan`). Cross-validation fits one design per fold
(:func:`cross_validate_cohort_cox`), so a validation fold is never encoded
with a value of its own. Before 2026-11-27 the median was computed on the
whole cohort before the folds were drawn, which is not what the paragraph
above said.

Note on future consolidation
----------------------------
`bootstrap_ci` and friends in `gliostrata.utils.stats` are being hardened
concurrently by another workstream, so the bootstrap used here is implemented
locally (``_bootstrap_ci``). Once `utils.stats` is stable these should be
consolidated onto that single implementation.
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import pandas as pd
from scipy import optimize, stats
from sklearn.model_selection import StratifiedKFold

__all__ = [
    "CoxPHModel",
    "DesignReport",
    "DesignPlan",
    "build_design_matrix",
    "fit_design_plan",
    "apply_design_plan",
    "concordance_index",
    "brier_score_ipcw",
    "cross_validate_cox",
    "cross_validate_cohort_cox",
    "fit_clinical_baseline",
    "DEFAULT_CLINICAL_FEATURES",
]

ImputeStrategy = Literal["median_indicator", "median_only", "drop"]

DEFAULT_CLINICAL_FEATURES: tuple[str, ...] = ("age", "sex", "kps", "treatment_arm")


# --------------------------------------------------------------------------
# Design matrix construction (explicit, reported missingness handling)
# --------------------------------------------------------------------------
@dataclass
class DesignReport:
    """Audit record of exactly what was done to build the design matrix."""

    n_input: int
    n_used: int
    dropped_missing_outcome: int
    feature_names: list[str] = field(default_factory=list)
    impute_strategy: str = "median_indicator"
    imputed_values: dict[str, float] = field(default_factory=dict)
    n_imputed: dict[str, int] = field(default_factory=dict)
    categorical_levels: dict[str, list[str]] = field(default_factory=dict)
    dropped_rows_missing_features: int = 0

    def summary(self) -> str:
        lines = [
            f"design matrix: {self.n_used}/{self.n_input} patients used "
            f"({self.dropped_missing_outcome} dropped for missing outcome, "
            f"{self.dropped_rows_missing_features} for missing features)",
            f"impute strategy: {self.impute_strategy}",
            f"features ({len(self.feature_names)}): {', '.join(self.feature_names)}",
        ]
        for col, n in sorted(self.n_imputed.items()):
            if not n:
                continue
            if col in self.imputed_values:
                lines.append(
                    f"  imputed {col}: {n} values -> median {self.imputed_values[col]:.4g}"
                    " (+ missingness indicator)"
                    if self.impute_strategy == "median_indicator"
                    else f"  imputed {col}: {n} values -> median "
                    f"{self.imputed_values[col]:.4g}"
                )
            else:
                lines.append(f"  {col}: {n} missing -> explicit 'missing' level")
        return "\n".join(lines)


def _is_numeric(series: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(series) or pd.api.types.is_bool_dtype(series)


def build_design_matrix(
    df: pd.DataFrame,
    features: Sequence[str] = DEFAULT_CLINICAL_FEATURES,
    *,
    time_col: str = "overall_survival_months",
    event_col: str = "vital_status",
    impute: ImputeStrategy = "median_indicator",
    min_time: float = 1e-8,
    reference_levels: dict[str, str] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, DesignReport]:
    """Build (X, time, event, report) from a clinical cohort frame.

    Patients with a missing outcome (time or event) are dropped -- they cannot
    contribute to a likelihood either way -- and the count is reported.
    Missing *features* are handled by `impute` and never dropped silently.

    The frame is both the training set and the frame being encoded: the
    medians, indicator columns, categorical levels and reference levels are
    fitted on its rows (:func:`fit_design_plan`) and applied to the same rows
    (:func:`apply_design_plan`). To encode one set of rows with a design fitted
    on another, as cross-validation must, use those two functions or
    :func:`cross_validate_cohort_cox`.

    `event` is 1 for a death and 0 for a right-censored observation.
    """
    outcome_ok = df[time_col].notna() & df[event_col].notna()
    plan = fit_design_plan(
        df.loc[outcome_ok],
        features,
        impute=impute,
        reference_levels=reference_levels,
    )
    return apply_design_plan(
        df, plan, time_col=time_col, event_col=event_col, min_time=min_time
    )


@dataclass(frozen=True)
class DesignPlan:
    """What a set of training rows fixes about the design matrix.

    Fitted on training rows only (:func:`fit_design_plan`) and applied to any
    rows (:func:`apply_design_plan`), so that a validation fold is encoded with
    the training fold's medians, missingness indicators, levels and reference
    levels, and with nothing of its own.
    """

    features: tuple[str, ...]
    numeric: tuple[str, ...]
    medians: dict[str, float]
    indicator_for: tuple[str, ...]
    categorical_levels: dict[str, tuple[str, ...]]
    reference_levels: dict[str, str]
    feature_names: tuple[str, ...]
    impute_strategy: str


def fit_design_plan(
    df: pd.DataFrame,
    features: Sequence[str] = DEFAULT_CLINICAL_FEATURES,
    *,
    impute: ImputeStrategy = "median_indicator",
    reference_levels: dict[str, str] | None = None,
) -> DesignPlan:
    """Fix the design on the rows of ``df``: the training rows.

    Every statistic is taken over these rows only: each numeric feature's
    median, whether it carries a missingness indicator column, each categorical
    feature's levels and reference level, and the resulting column order.

    Raises
    ------
    ValueError
        If a feature is absent from the frame, is entirely missing in it, or
        has a reference level that does not occur in it.
    """
    refs = reference_levels or {}
    numeric: list[str] = []
    medians: dict[str, float] = {}
    indicator_for: list[str] = []
    levels_by_feature: dict[str, tuple[str, ...]] = {}
    references: dict[str, str] = {}
    names: list[str] = []

    for feat in features:
        if feat not in df.columns:
            raise ValueError(f"feature {feat!r} not in dataframe")
        series = df[feat]
        if series.isna().all():
            # e.g. the molecular columns, which are not ingested yet.
            raise ValueError(
                f"feature {feat!r} is entirely missing in this cohort; "
                "exclude it rather than imputing a constant"
            )

        if _is_numeric(series):
            values = series.to_numpy(dtype=float)
            numeric.append(feat)
            medians[feat] = float(np.nanmedian(values))
            if np.isnan(values).any() and impute == "median_indicator":
                indicator_for.append(feat)
                names.append(f"{feat}_missing")
            names.append(feat)
        else:
            raw = series.astype("object").where(series.notna(), "missing")
            levels = tuple(sorted({str(v) for v in raw}))
            ref = refs.get(feat)
            if ref is None:
                # Most frequent level is the reference, so HRs read naturally.
                ref = str(raw.astype(str).value_counts().idxmax())
            if ref not in levels:
                raise ValueError(f"reference level {ref!r} not present for {feat}")
            levels_by_feature[feat] = levels
            references[feat] = ref
            names.extend(f"{feat}_{lvl}" for lvl in levels if lvl != ref)

    return DesignPlan(
        features=tuple(features),
        numeric=tuple(numeric),
        medians=medians,
        indicator_for=tuple(indicator_for),
        categorical_levels=levels_by_feature,
        reference_levels=references,
        feature_names=tuple(names),
        impute_strategy=impute,
    )


def apply_design_plan(
    df: pd.DataFrame,
    plan: DesignPlan,
    *,
    time_col: str = "overall_survival_months",
    event_col: str = "vital_status",
    min_time: float = 1e-8,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, DesignReport]:
    """Encode ``df`` with a design fitted on another set of rows.

    A missing numeric value is filled with the training median, whether or not
    the rows being encoded had any. A missingness indicator column exists only
    where the training rows had one. A categorical value the training rows
    never saw, ``missing`` included, is encoded as the reference level. The
    column order is the plan's, so the matrix suits a model fitted on the
    training rows.

    Patients with a missing outcome (time or event) are dropped and counted.
    ``event`` is 1 for a death and 0 for a right-censored observation.
    """
    outcome_ok = df[time_col].notna() & df[event_col].notna()
    dropped_outcome = int((~outcome_ok).sum())
    work = df.loc[outcome_ok].reset_index(drop=True)

    time = work[time_col].to_numpy(dtype=float)
    event = work[event_col].to_numpy(dtype=float)
    if not np.all(np.isin(event, (0.0, 1.0))):
        raise ValueError(f"{event_col} must be 0 (censored/alive) or 1 (death)")
    event = event.astype(int)
    time = np.maximum(time, min_time)  # guard against exact-zero follow-up

    report = DesignReport(
        n_input=len(df),
        n_used=len(work),
        dropped_missing_outcome=dropped_outcome,
        impute_strategy=plan.impute_strategy,
    )
    cols: list[np.ndarray] = [np.zeros(len(work)) for _ in plan.feature_names]

    for feat in plan.features:
        series = work[feat]
        if feat in plan.numeric:
            values = series.to_numpy(dtype=float)
            n_missing = int(np.isnan(values).sum())
            report.n_imputed[feat] = n_missing
            if n_missing:
                if plan.impute_strategy == "drop":
                    raise ValueError(
                        "impute='drop' requires pre-filtered data; "
                        f"{feat} has {n_missing} missing values"
                    )
                fill = plan.medians[feat]
                report.imputed_values[feat] = fill
                values = np.where(np.isnan(values), fill, values)
            if feat in plan.indicator_for:
                cols[plan.feature_names.index(f"{feat}_missing")] = series.isna().to_numpy(
                    dtype=float
                )
            cols[plan.feature_names.index(feat)] = values
        else:
            raw = series.astype("object").where(series.notna(), "missing").astype(str)
            levels = plan.categorical_levels[feat]
            report.categorical_levels[feat] = list(levels)
            report.n_imputed[feat] = int(series.isna().sum())
            for lvl in levels:
                if lvl == plan.reference_levels[feat]:
                    continue
                cols[plan.feature_names.index(f"{feat}_{lvl}")] = (
                    (raw == lvl).to_numpy(dtype=float)
                )

    x = np.column_stack(cols) if cols else np.zeros((len(work), 0))
    report.feature_names = list(plan.feature_names)
    return x, time, event, report


# --------------------------------------------------------------------------
# Breslow partial likelihood
# --------------------------------------------------------------------------
def _risk_set_structure(time: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (order, first_index_of_tied_block) for ascending-sorted time.

    `order` sorts time ascending. For position k in sorted order,
    `block_start[k]` is the first position sharing that same time, so the
    Breslow risk set {j : t_j >= t_k} is positions block_start[k]..n-1.
    """
    order = np.argsort(time, kind="mergesort")
    t_sorted = time[order]
    n = len(t_sorted)
    block_start = np.empty(n, dtype=int)
    k = 0
    while k < n:
        j = k
        while j + 1 < n and t_sorted[j + 1] == t_sorted[k]:
            j += 1
        block_start[k : j + 1] = k
        k = j + 1
    return order, block_start


def _neg_log_partial_likelihood(
    beta: np.ndarray,
    X: np.ndarray,
    event: np.ndarray,
    block_start: np.ndarray,
    l2: float,
    want_hessian: bool = False,
):
    """Negative Breslow partial log-likelihood, gradient, and (optionally) Hessian.

    X, event are assumed already sorted by ascending time.
    """
    n, p = X.shape
    eta = X @ beta
    eta = eta - eta.max()  # numerical stability; cancels in the ratio
    w = np.exp(eta)

    # Reverse cumulative sums give risk-set aggregates in O(n p).
    rev_w = np.cumsum(w[::-1])[::-1]
    rev_wx = np.cumsum((w[:, None] * X)[::-1], axis=0)[::-1]

    ev = event.astype(bool)
    R = rev_w[block_start[ev]]
    A = rev_wx[block_start[ev]] / R[:, None]

    loglik = np.sum(eta[ev] - np.log(R))
    grad = np.sum(X[ev] - A, axis=0)

    nll = -loglik + 0.5 * l2 * float(beta @ beta)
    ngrad = -grad + l2 * beta
    if not want_hessian:
        return nll, ngrad

    wxx = w[:, None, None] * (X[:, :, None] * X[:, None, :])
    rev_wxx = np.cumsum(wxx[::-1], axis=0)[::-1]
    B = rev_wxx[block_start[ev]] / R[:, None, None]
    hess = np.sum(B - A[:, :, None] * A[:, None, :], axis=0)
    nhess = hess + l2 * np.eye(p)
    return nll, ngrad, nhess


@dataclass
class CoxPHModel:
    """Cox proportional-hazards model fitted by Breslow partial likelihood.

    Parameters
    ----------
    l2:
        Ridge penalty on the (standardised) coefficients. A small non-zero
        value keeps the fit well-posed under near-separation / rare dummy
        levels. Default 1e-6 is numerically negligible for inference.
    standardize:
        Centre and scale features before optimisation. Cox PH has no
        intercept and the partial likelihood is invariant to centring, so
        coefficients are exactly rescaled back to original units.
    """

    l2: float = 1e-6
    standardize: bool = True
    max_iter: int = 200
    tol: float = 1e-9

    # fitted state
    coef_: np.ndarray | None = None
    feature_names_: list[str] | None = None
    se_: np.ndarray | None = None
    baseline_times_: np.ndarray | None = None
    baseline_cumhaz_: np.ndarray | None = None
    center_: np.ndarray | None = None
    scale_: np.ndarray | None = None
    log_likelihood_: float | None = None
    n_obs_: int | None = None
    n_events_: int | None = None
    converged_: bool | None = None

    # ---------------- fitting ----------------
    def fit(
        self,
        X: np.ndarray,
        time: np.ndarray,
        event: np.ndarray,
        feature_names: Sequence[str] | None = None,
    ) -> CoxPHModel:
        X = np.asarray(X, dtype=float)
        time = np.asarray(time, dtype=float)
        event = np.asarray(event).astype(int)
        if not np.all(np.isin(event, (0, 1))):
            raise ValueError("event must be 0 (censored) or 1 (death)")
        if event.sum() == 0:
            raise ValueError("cannot fit a Cox model with zero events")
        n, p = X.shape

        self.center_ = X.mean(axis=0) if self.standardize else np.zeros(p)
        if self.standardize:
            sd = X.std(axis=0)
            sd[sd < 1e-12] = 1.0
            self.scale_ = sd
        else:
            self.scale_ = np.ones(p)
        Z = (X - self.center_) / self.scale_

        order, block_start = _risk_set_structure(time)
        Zs, ev_s = Z[order], event[order]

        def fun(b):
            return _neg_log_partial_likelihood(b, Zs, ev_s, block_start, self.l2)

        def hess(b):
            return _neg_log_partial_likelihood(
                b, Zs, ev_s, block_start, self.l2, want_hessian=True
            )[2]

        res = optimize.minimize(
            fun,
            np.zeros(p),
            jac=True,
            hess=hess,
            method="trust-exact",
            options={"maxiter": self.max_iter, "gtol": 1e-8},
        )
        beta_std = res.x

        # Newton-Raphson polish: the trust-region solver stops on its own
        # internal criteria, which on a well-conditioned partial likelihood
        # leaves a small but non-zero gradient. Exact Newton steps with a
        # backtracking line search drive it to machine precision, which is
        # what makes the fit (and therefore every downstream metric) bitwise
        # reproducible.
        nll, ngrad, nhess = _neg_log_partial_likelihood(
            beta_std, Zs, ev_s, block_start, self.l2, want_hessian=True
        )
        for _ in range(self.max_iter):
            if np.max(np.abs(ngrad)) <= self.tol * max(1.0, float(event.sum())):
                break
            try:
                step = np.linalg.solve(nhess, ngrad)
            except np.linalg.LinAlgError:  # pragma: no cover
                break
            improved = False
            scale = 1.0
            for _ls in range(30):
                cand = beta_std - scale * step
                c_nll, c_grad, c_hess = _neg_log_partial_likelihood(
                    cand, Zs, ev_s, block_start, self.l2, want_hessian=True
                )
                if c_nll <= nll:
                    beta_std, nll, ngrad, nhess = cand, c_nll, c_grad, c_hess
                    improved = True
                    break
                scale *= 0.5
            if not improved:
                break
        self.converged_ = bool(
            np.max(np.abs(ngrad)) <= self.tol * max(1.0, float(event.sum()))
        )
        self.log_likelihood_ = float(-nll)

        # Back out original-scale coefficients (centring is a no-op for Cox).
        self.coef_ = beta_std / self.scale_
        try:
            cov_std = np.linalg.inv(nhess)
            self.se_ = np.sqrt(np.clip(np.diag(cov_std), 0, None)) / self.scale_
        except np.linalg.LinAlgError:  # pragma: no cover - degenerate design
            self.se_ = np.full(p, np.nan)

        self.feature_names_ = (
            list(feature_names) if feature_names is not None
            else [f"x{i}" for i in range(p)]
        )
        self.n_obs_, self.n_events_ = int(n), int(event.sum())
        self._fit_baseline(X, time, event)
        return self

    def _fit_baseline(self, X: np.ndarray, time: np.ndarray, event: np.ndarray) -> None:
        """Breslow estimator of the baseline cumulative hazard.

        The baseline is defined at the *mean covariate vector* (i.e. the linear
        predictor used everywhere here is ``(x - center) @ beta``), which keeps
        exp(eta) near 1 and the baseline survival curve interpretable.
        """
        eta = self.predict_log_hazard(X)
        w = np.exp(eta)
        order = np.argsort(time, kind="mergesort")
        t_s, ev_s, w_s = time[order], event[order], w[order]
        rev_w = np.cumsum(w_s[::-1])[::-1]
        _, block_start = _risk_set_structure(time)

        times, increments = [], []
        n = len(t_s)
        k = 0
        while k < n:
            j = k
            while j + 1 < n and t_s[j + 1] == t_s[k]:
                j += 1
            d = int(ev_s[k : j + 1].sum())
            if d > 0:
                times.append(t_s[k])
                increments.append(d / rev_w[block_start[k]])
            k = j + 1
        self.baseline_times_ = np.asarray(times, dtype=float)
        self.baseline_cumhaz_ = np.cumsum(increments) if times else np.zeros(0)

    # ---------------- prediction ----------------
    def _check_fitted(self) -> None:
        if self.coef_ is None:
            raise RuntimeError("model is not fitted")

    def predict_log_hazard(self, X: np.ndarray) -> np.ndarray:
        """Linear predictor (log relative hazard vs the mean patient)."""
        self._check_fitted()
        X = np.asarray(X, dtype=float)
        return (X - self.center_) @ self.coef_

    # risk score == log hazard; higher means worse prognosis
    predict_risk = predict_log_hazard

    def baseline_cumulative_hazard(self, t: np.ndarray | float) -> np.ndarray:
        """Breslow baseline cumulative hazard, right-continuous step function."""
        self._check_fitted()
        t = np.atleast_1d(np.asarray(t, dtype=float))
        if len(self.baseline_times_) == 0:  # pragma: no cover
            return np.zeros_like(t)
        idx = np.searchsorted(self.baseline_times_, t, side="right") - 1
        out = np.where(idx < 0, 0.0, self.baseline_cumhaz_[np.clip(idx, 0, None)])
        return out

    def predict_survival_function(
        self, X: np.ndarray, times: np.ndarray | float
    ) -> np.ndarray:
        """S(t | x) = exp(-H0(t) * exp(eta)); shape (n_samples, n_times)."""
        self._check_fitted()
        times = np.atleast_1d(np.asarray(times, dtype=float))
        H0 = self.baseline_cumulative_hazard(times)
        risk = np.exp(self.predict_log_hazard(X))
        return np.exp(-np.outer(risk, H0))

    def predict_median_survival(self, X: np.ndarray) -> np.ndarray:
        """Median survival time per patient (NaN if the curve never hits 0.5)."""
        self._check_fitted()
        if len(self.baseline_times_) == 0:  # pragma: no cover
            return np.full(len(X), np.nan)
        S = self.predict_survival_function(X, self.baseline_times_)
        out = np.full(S.shape[0], np.nan)
        for i in range(S.shape[0]):
            hit = np.nonzero(S[i] <= 0.5)[0]
            if hit.size:
                out[i] = self.baseline_times_[hit[0]]
        return out

    # ---------------- reporting ----------------
    def summary(self, alpha: float = 0.05) -> pd.DataFrame:
        """Coefficients, hazard ratios, Wald CIs and p-values."""
        self._check_fitted()
        z = stats.norm.ppf(1 - alpha / 2)
        coef, se = self.coef_, self.se_
        with np.errstate(divide="ignore", invalid="ignore"):
            zstat = coef / se
        return pd.DataFrame(
            {
                "feature": self.feature_names_,
                "coef": coef,
                "se": se,
                "hazard_ratio": np.exp(coef),
                "hr_lower": np.exp(coef - z * se),
                "hr_upper": np.exp(coef + z * se),
                "z": zstat,
                "p": 2 * stats.norm.sf(np.abs(zstat)),
            }
        )


# --------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------
def concordance_index(
    time: np.ndarray, event: np.ndarray, risk: np.ndarray
) -> float:
    """Harrell's C-index. `risk` is higher-is-worse (e.g. the log hazard).

    Only *comparable* pairs count: a pair (i, j) is comparable when the one
    with the shorter time had an event, so a censored patient is never assumed
    to have died at censoring. Ties in risk score score 0.5.
    """
    time = np.asarray(time, dtype=float)
    event = np.asarray(event).astype(int)
    risk = np.asarray(risk, dtype=float)
    if not (len(time) == len(event) == len(risk)):
        raise ValueError("time, event and risk must be the same length")

    ti, tj = time[:, None], time[None, :]
    ei = event[:, None]
    # i is the earlier-failing, uncensored member of the pair
    comparable = (ti < tj) & (ei == 1)
    n_comp = comparable.sum()
    if n_comp == 0:
        raise ValueError("no comparable pairs; cannot compute a C-index")
    ri, rj = risk[:, None], risk[None, :]
    concordant = (comparable & (ri > rj)).sum()
    tied = (comparable & (ri == rj)).sum()
    return float((concordant + 0.5 * tied) / n_comp)


def _km_censoring_survival(time: np.ndarray, event: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Kaplan-Meier estimate G(t) of the CENSORING distribution.

    Censoring is the 'event' here (event == 0), which is what IPCW needs.
    """
    order = np.argsort(time, kind="mergesort")
    t_s, cens_s = time[order], (np.asarray(event) == 0).astype(int)[order]
    n = len(t_s)
    times, surv = [], []
    g = 1.0
    k = 0
    while k < n:
        j = k
        while j + 1 < n and t_s[j + 1] == t_s[k]:
            j += 1
        at_risk = n - k
        d = int(cens_s[k : j + 1].sum())
        if d > 0:
            g *= 1.0 - d / at_risk
            times.append(t_s[k])
            surv.append(g)
        k = j + 1
    return np.asarray(times, dtype=float), np.asarray(surv, dtype=float)


def _eval_step(times: np.ndarray, values: np.ndarray, t: np.ndarray) -> np.ndarray:
    t = np.atleast_1d(np.asarray(t, dtype=float))
    if len(times) == 0:
        return np.ones_like(t)
    idx = np.searchsorted(times, t, side="right") - 1
    return np.where(idx < 0, 1.0, values[np.clip(idx, 0, None)])


def brier_score_ipcw(
    time: np.ndarray,
    event: np.ndarray,
    survival_prob: np.ndarray,
    horizon: float,
    *,
    censoring_time: np.ndarray | None = None,
    censoring_event: np.ndarray | None = None,
    eps: float = 1e-8,
) -> float:
    """Inverse-probability-of-censoring-weighted Brier score at `horizon`.

    `survival_prob[i]` is the model's predicted S(horizon | x_i).

    Graf et al. (1999): patients censored before the horizon contribute
    nothing, deaths before the horizon are weighted by 1/G(T_i-) and
    survivors past the horizon by 1/G(horizon). Lower is better; 0.25 is the
    score of an uninformative 0.5 prediction.
    """
    time = np.asarray(time, dtype=float)
    event = np.asarray(event).astype(int)
    S = np.asarray(survival_prob, dtype=float).ravel()
    if len(S) != len(time):
        raise ValueError("survival_prob must have one entry per patient")

    ct = time if censoring_time is None else np.asarray(censoring_time, dtype=float)
    ce = event if censoring_event is None else np.asarray(censoring_event).astype(int)
    g_times, g_vals = _km_censoring_survival(ct, ce)

    g_at_ti = _eval_step(g_times, g_vals, time)
    g_at_h = float(_eval_step(g_times, g_vals, np.array([horizon]))[0])

    died_before = (time <= horizon) & (event == 1)
    alive_after = time > horizon

    contrib = np.zeros_like(S)
    contrib[died_before] = (S[died_before] ** 2) / np.maximum(g_at_ti[died_before], eps)
    contrib[alive_after] = ((1.0 - S[alive_after]) ** 2) / max(g_at_h, eps)
    return float(contrib.sum() / len(S))


def _bootstrap_ci(
    stat_fn,
    n: int,
    *,
    n_boot: int = 1000,
    alpha: float = 0.05,
    seed: int = 42,
) -> tuple[float, float]:
    """Percentile bootstrap CI over patient indices.

    Local implementation on purpose: `gliostrata.utils.stats.bootstrap_ci` is
    being hardened concurrently. TODO(consolidate): switch to utils.stats once
    that API is frozen.
    """
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        try:
            vals.append(stat_fn(idx))
        except ValueError:
            continue
    if not vals:  # pragma: no cover
        return (float("nan"), float("nan"))
    v = np.asarray(vals, dtype=float)
    return (
        float(np.quantile(v, alpha / 2)),
        float(np.quantile(v, 1 - alpha / 2)),
    )


# --------------------------------------------------------------------------
# Cross-validation
# --------------------------------------------------------------------------
@dataclass
class CVResult:
    c_index: float
    c_index_ci: tuple[float, float]
    fold_c_indices: list[float]
    brier_horizon: float
    brier_score: float
    brier_ci: tuple[float, float]
    n_obs: int
    n_events: int
    oof_risk: np.ndarray
    oof_survival_at_horizon: np.ndarray
    n_splits: int
    seed: int

    def summary(self) -> str:
        lo, hi = self.c_index_ci
        blo, bhi = self.brier_ci
        return (
            f"{self.n_splits}-fold CV (seed={self.seed}) on n={self.n_obs} "
            f"({self.n_events} events)\n"
            f"  C-index = {self.c_index:.4f}  95% CI [{lo:.4f}, {hi:.4f}]\n"
            f"  per-fold C = {[round(c, 4) for c in self.fold_c_indices]}\n"
            f"  Brier@{self.brier_horizon:g}mo (IPCW) = {self.brier_score:.4f} "
            f"95% CI [{blo:.4f}, {bhi:.4f}]"
        )


def cross_validate_cox(
    X: np.ndarray,
    time: np.ndarray,
    event: np.ndarray,
    *,
    n_splits: int = 5,
    seed: int = 42,
    horizon: float = 18.0,
    n_boot: int = 1000,
    model_kwargs: dict | None = None,
) -> CVResult:
    """Stratified k-fold CV (stratified on event status, fixed seed).

    Out-of-fold risk scores are pooled and scored once, which avoids the
    optimistic bias of averaging within-fold metrics on small folds; per-fold
    C-indices are also reported. Confidence intervals come from a patient-level
    percentile bootstrap of the pooled out-of-fold predictions.
    """
    X = np.asarray(X, dtype=float)
    time = np.asarray(time, dtype=float)
    event = np.asarray(event).astype(int)
    kwargs = model_kwargs or {}

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    oof_risk = np.full(len(time), np.nan)
    oof_surv = np.full(len(time), np.nan)
    fold_c: list[float] = []

    for train_idx, test_idx in skf.split(X, event):
        m = CoxPHModel(**kwargs).fit(X[train_idx], time[train_idx], event[train_idx])
        r = m.predict_log_hazard(X[test_idx])
        oof_risk[test_idx] = r
        oof_surv[test_idx] = m.predict_survival_function(X[test_idx], horizon)[:, 0]
        try:
            fold_c.append(concordance_index(time[test_idx], event[test_idx], r))
        except ValueError:  # pragma: no cover
            pass

    return _pooled_cv_result(
        time,
        event,
        oof_risk,
        oof_surv,
        fold_c,
        horizon=horizon,
        n_boot=n_boot,
        seed=seed,
        n_splits=n_splits,
    )


def _pooled_cv_result(
    time: np.ndarray,
    event: np.ndarray,
    oof_risk: np.ndarray,
    oof_surv: np.ndarray,
    fold_c: list[float],
    *,
    horizon: float,
    n_boot: int,
    seed: int,
    n_splits: int,
) -> CVResult:
    """Score pooled out-of-fold predictions once, with percentile bootstrap CIs."""
    c = concordance_index(time, event, oof_risk)
    c_ci = _bootstrap_ci(
        lambda idx: concordance_index(time[idx], event[idx], oof_risk[idx]),
        len(time),
        n_boot=n_boot,
        seed=seed,
    )
    b = brier_score_ipcw(time, event, oof_surv, horizon)
    b_ci = _bootstrap_ci(
        lambda idx: brier_score_ipcw(time[idx], event[idx], oof_surv[idx], horizon),
        len(time),
        n_boot=n_boot,
        seed=seed,
    )
    return CVResult(
        c_index=c,
        c_index_ci=c_ci,
        fold_c_indices=fold_c,
        brier_horizon=horizon,
        brier_score=b,
        brier_ci=b_ci,
        n_obs=len(time),
        n_events=int(event.sum()),
        oof_risk=oof_risk,
        oof_survival_at_horizon=oof_surv,
        n_splits=n_splits,
        seed=seed,
    )


def cross_validate_cohort_cox(
    frame: pd.DataFrame,
    features: Sequence[str] = DEFAULT_CLINICAL_FEATURES,
    *,
    impute: ImputeStrategy = "median_indicator",
    reference_levels: dict[str, str] | None = None,
    time_col: str = "overall_survival_months",
    event_col: str = "vital_status",
    n_splits: int = 5,
    seed: int = 42,
    horizon: float = 18.0,
    n_boot: int = 1000,
    model_kwargs: dict | None = None,
) -> tuple[CVResult, list[DesignReport]]:
    """Stratified k-fold CV in which each fold's design is fitted on its training rows.

    The folds are those :func:`cross_validate_cox` draws on the same patients
    (stratified on event status, same seed), and the pooled out-of-fold scores
    are computed the same way. What differs is where the encoding comes from:
    each training fold fixes its own medians, indicator columns, categorical
    levels and reference levels (:func:`fit_design_plan`), and its validation
    fold is encoded with those and nothing of its own. When nothing has to be
    imputed the two functions give identical predictions.

    Returns
    -------
    (CVResult, reports)
        ``reports`` holds, per fold, the :class:`DesignReport` of the validation
        fold. Its ``imputed_values`` are the training fold's medians.
    """
    outcome_ok = frame[time_col].notna() & frame[event_col].notna()
    rows = frame.loc[outcome_ok].reset_index(drop=True)
    # The outcome arrays are exactly those build_design_matrix would return.
    _, time, event, _ = apply_design_plan(
        frame,
        fit_design_plan(rows, features, impute=impute, reference_levels=reference_levels),
        time_col=time_col,
        event_col=event_col,
    )

    kwargs = model_kwargs or {}
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    oof_risk = np.full(len(time), np.nan)
    oof_surv = np.full(len(time), np.nan)
    fold_c: list[float] = []
    fold_reports: list[DesignReport] = []

    for train_idx, test_idx in skf.split(np.zeros(len(time)), event):
        train_rows = rows.iloc[train_idx]
        plan = fit_design_plan(
            train_rows, features, impute=impute, reference_levels=reference_levels
        )
        x_train, time_train, event_train, _ = apply_design_plan(
            train_rows, plan, time_col=time_col, event_col=event_col
        )
        x_test, time_test, event_test, report = apply_design_plan(
            rows.iloc[test_idx], plan, time_col=time_col, event_col=event_col
        )
        m = CoxPHModel(**kwargs).fit(x_train, time_train, event_train)
        r = m.predict_log_hazard(x_test)
        oof_risk[test_idx] = r
        oof_surv[test_idx] = m.predict_survival_function(x_test, horizon)[:, 0]
        try:
            fold_c.append(concordance_index(time_test, event_test, r))
        except ValueError:  # pragma: no cover
            pass
        fold_reports.append(report)

    cv = _pooled_cv_result(
        time,
        event,
        oof_risk,
        oof_surv,
        fold_c,
        horizon=horizon,
        n_boot=n_boot,
        seed=seed,
        n_splits=n_splits,
    )
    return cv, fold_reports


def fit_clinical_baseline(
    cohort: pd.DataFrame,
    *,
    features: Sequence[str] = DEFAULT_CLINICAL_FEATURES,
    impute: ImputeStrategy = "median_indicator",
    n_splits: int = 5,
    seed: int = 42,
    horizon: float = 18.0,
    n_boot: int = 1000,
    model_kwargs: dict | None = None,
) -> tuple[CoxPHModel, CVResult, DesignReport]:
    """Fit + cross-validate the honest clinical-only baseline on a cohort.

    Returns the full-data model (for coefficients / survival curves), the
    cross-validated metrics (for honest performance), and the design report
    (for the missingness audit trail).

    Each cross-validation fold is encoded with a design fitted on that fold's
    training rows (:func:`cross_validate_cohort_cox`): the imputation medians,
    the missingness indicators and the categorical levels come from the
    training fold, never from the cohort as a whole. The returned model and
    report are the full-data fit, in which every modelled patient is training.
    """
    cv, _ = cross_validate_cohort_cox(
        cohort,
        features,
        impute=impute,
        n_splits=n_splits,
        seed=seed,
        horizon=horizon,
        n_boot=n_boot,
        model_kwargs=model_kwargs,
    )
    X, time, event, report = build_design_matrix(cohort, features, impute=impute)
    model = CoxPHModel(**(model_kwargs or {})).fit(
        X, time, event, feature_names=report.feature_names
    )
    return model, cv, report
