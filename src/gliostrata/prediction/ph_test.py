"""Proportional-hazards check for a fitted Cox model (Grambsch and Therneau, 1994).

Pre-registration section 3.8: every Cox model used in a registered read is reported with
this test. A violation is reported as a limitation. It never changes the model, the
verdict or the analysis.

The statistic is the one computed by R's ``survival::cox.zph`` (survival 2.x), on the
scaled Schoenfeld residuals:

* Schoenfeld residual at the k-th death: r_k = x_k - xbar(t_k). xbar is the risk-set mean
  of the covariates, weighted by exp(x beta), over everyone still at risk at t_k
  (Breslow: tied deaths share one risk set).
* V is the inverse of the Breslow information matrix at beta, and d the number of deaths.
  The scaled residuals are s_k = beta + d V r_k.
* g_k is the rank of the k-th death time among the deaths (average ranks for ties), and
  gc_k = g_k - mean(g).
* Per covariate j:  z_j = (sum_k gc_k s_kj)^2 / (d V_jj sum_k gc_k^2), chi-square, 1 df.
* Globally:  z = (gc' R) V (R' gc) d / sum_k gc_k^2, where R is the matrix of residuals
  r_k, chi-square, p df.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy import stats

__all__ = ["proportional_hazards_test"]


def _rank_average(values: np.ndarray) -> np.ndarray:
    return stats.rankdata(values, method="average")


def proportional_hazards_test(x_mat, time, event, beta,
                              feature_names: Sequence[str] | None = None) -> dict:
    """Per-covariate and global Grambsch-Therneau tests, with rank-of-time transform."""
    x_mat = np.asarray(x_mat, dtype=float)
    time = np.asarray(time, dtype=float)
    event = np.asarray(event).astype(int)
    beta = np.asarray(beta, dtype=float)
    n, p = x_mat.shape
    if beta.shape != (p,):
        raise ValueError("beta must have one coefficient per column")
    names = list(feature_names) if feature_names is not None else [f"x{j}" for j in range(p)]

    order = np.argsort(time, kind="mergesort")
    t_sorted = time[order]
    w = np.exp(x_mat @ beta)
    death_idx = np.flatnonzero(event == 1)
    d = len(death_idx)
    if d < 2:
        raise ValueError("fewer than two deaths; the test is undefined")

    resid = np.empty((d, p))
    info = np.zeros((p, p))
    for k, i in enumerate(death_idx):
        # Risk set at t_i: everyone with time >= t_i (Breslow).
        start = np.searchsorted(t_sorted, time[i], side="left")
        at_risk = order[start:]
        wr = w[at_risk]
        xr = x_mat[at_risk]
        s0 = wr.sum()
        xbar = (wr @ xr) / s0
        resid[k] = x_mat[i] - xbar
        centred = xr - xbar
        info += (centred * wr[:, None]).T @ centred / s0
    var = np.linalg.inv(info)

    g = _rank_average(time[death_idx])
    gc = g - g.mean()
    sgc2 = float(gc @ gc)
    scaled = beta[None, :] + d * resid @ var
    per = []
    for j in range(p):
        num = float(gc @ scaled[:, j]) ** 2
        z = num / (d * var[j, j] * sgc2)
        rho = float(np.corrcoef(gc, scaled[:, j])[0, 1]) if np.std(scaled[:, j]) > 0 else float("nan")
        per.append({"term": names[j], "rho": rho, "chisq": z, "df": 1,
                    "p": float(stats.chi2.sf(z, 1))})
    tvec = gc @ resid
    z_glob = float(tvec @ var @ tvec) * d / sgc2
    return {"transform": "rank", "n": n, "n_deaths": d, "terms": per,
            "global": {"chisq": z_glob, "df": p, "p": float(stats.chi2.sf(z_glob, p))}}
