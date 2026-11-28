"""L4-lite: the closed-form delivery model registered by pre-registration item 11(h).

Item 11(h) names this module "L4-lite: closed-form delivery model, uncalibrated,
literature-parameterised". It replaces the stubs of `gliostrata.simulation.nsc`, whose
"affinity" table had no source. It is **declared uncalibrated** and that declaration stands
whatever T3 returns: T3 data would calibrate the acting fraction, the tropism intercept and
slope and the coverage constant, but T4 and T5 do not require them (item 11(h)).

For patient i and configuration c = (carrier k, payload p, dose n_cells, sites n, schedule s):

    z_i,k   = mean over the carrier's tropism genes of the standardised expression
    n_eff   = n_cells * f_k * logistic(a_k + g_k * z_i,k)
    cov(n)  = 1 - exp(-n / n0)
    sigma(s) = the registered schedule factor

Every constant comes from `configs/response_surface_v1.yaml` through
`gliostrata.simulation.parameters`; this module defines none of its own.

**The standardisation reference is a frozen TCGA mean and SD, never a CGGA cohort
statistic.** A per-cohort standardisation would let the sealed external cohort leak into the
transform, which section 3.1 forbids, so `tropism_score` takes the reference explicitly and
refuses to compute one from the data it is given.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from gliostrata.simulation.parameters import ParameterError, ParameterSet


@dataclass(frozen=True)
class TropismReference:
    """A frozen per-gene mean and standard deviation, estimated on TCGA and never on CGGA."""

    mean: dict[str, float]
    sd: dict[str, float]
    source: str

    def __post_init__(self) -> None:
        missing = set(self.mean) ^ set(self.sd)
        if missing:
            raise ParameterError(f"tropism reference mean and sd cover different genes: {sorted(missing)}")
        bad = sorted(g for g, s in self.sd.items() if not np.isfinite(s) or s <= 0)
        if bad:
            raise ParameterError(f"tropism reference has a non-positive standard deviation: {bad}")
        if "cgga" in self.source.lower():
            raise ParameterError(
                "the tropism reference names CGGA as its source. Section 3.1 and item 11 "
                "require a frozen TCGA reference: a cohort statistic taken from the sealed "
                "external cohort would leak it into the transform.")


def tropism_score(
    expression: pd.DataFrame,
    carrier: str,
    params: ParameterSet,
    reference: TropismReference,
) -> pd.Series:
    """Standardised tropism score per patient, for one carrier.

    `expression` is genes x patients on the log2 CPM scale. Genes are matched by exact
    symbol, as elsewhere in the project. A gene absent from the matrix or from the reference
    is dropped for every patient, and the count of genes actually used is attached to the
    result's `attrs` so that it is reported rather than hidden.
    """
    genes = params.tropism_genes(carrier)
    usable = [g for g in genes if g in expression.index and g in reference.mean]
    if not usable:
        raise ParameterError(
            f"none of {carrier}'s tropism genes {genes} is present in both the expression "
            "matrix and the frozen reference; the score cannot be computed")
    standardised = np.vstack([
        (expression.loc[gene].to_numpy(dtype=float) - reference.mean[gene]) / reference.sd[gene]
        for gene in usable
    ])
    # A missing value drops that gene for that patient only, as T1's rule does.
    with np.errstate(invalid="ignore"):
        score = np.nanmean(standardised, axis=0)
    out = pd.Series(score, index=expression.columns, name=f"tropism_{carrier}")
    out.attrs["genes_registered"] = list(genes)
    out.attrs["genes_used"] = usable
    out.attrs["n_genes_used"] = len(usable)
    return out


def _logistic(x: np.ndarray | float) -> np.ndarray | float:
    return 1.0 / (1.0 + np.exp(-np.asarray(x, dtype=float)))


def effective_dose(
    dose_cells: float | np.ndarray,
    carrier: str,
    tropism: float | np.ndarray,
    params: ParameterSet,
) -> np.ndarray:
    """n_eff = dose * f_k * logistic(a_k + g_k * z). The engrafted, acting cell count."""
    dose = np.asarray(dose_cells, dtype=float)
    if np.any(dose < 0):
        raise ValueError("dose_cells cannot be negative")
    fraction = params.cell_fraction_acting(carrier)
    logit = params.tropism_logit_intercept + params.tropism_logit_slope * np.asarray(
        tropism, dtype=float)
    return dose * fraction * _logistic(logit)


def coverage(n_injection_sites: int | np.ndarray, params: ParameterSet) -> np.ndarray:
    """cov(n) = 1 - exp(-n / n0). Saturating spatial coverage in the number of sites."""
    sites = np.asarray(n_injection_sites, dtype=float)
    if np.any(sites < 1):
        raise ValueError("n_injection_sites starts at 1 (item 4's registered range is 1 to 8)")
    return 1.0 - np.exp(-sites / params.coverage_saturation_sites)


def schedule_factor(schedule: str, params: ParameterSet) -> float:
    """The registered multiplier for a single or fractionated schedule."""
    return params.schedule_factor(schedule)


def delivery_terms(
    dose_cells: float | np.ndarray,
    carrier: str,
    tropism: float | np.ndarray,
    n_injection_sites: int | np.ndarray,
    schedule: str,
    params: ParameterSet,
) -> dict[str, np.ndarray | float]:
    """Every delivery quantity the response surface needs, in one call."""
    return {
        "n_eff": effective_dose(dose_cells, carrier, tropism, params),
        "coverage": coverage(n_injection_sites, params),
        "schedule_factor": schedule_factor(schedule, params),
        "calibrated": False,  # item 11(h), and not negotiable by configuration
    }
