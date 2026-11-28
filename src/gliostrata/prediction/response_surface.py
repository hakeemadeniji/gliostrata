"""The T4/T5 response surface registered by pre-registration item 11 (Design A).

The survival model (item 11(g), variant 2: age, sex, arm) supplies the Stupp baseline curve.
This module turns a configuration into a predicted curve by a log-hazard offset:

    E_i(c) = Emax_p * m_p(i) * sigma(s) * cov(n) * n_eff^h / (n_eff^h + ED50^h)
    T_i(c) = tau_p * (dose / dose_ref)^q * (1 + eta * (n - 1)) + delta_delay
    S_c,i(t) = S_stupp,i(t) ^ exp(-E_i(c) + T_i(c))
    RMST_i(c) = integral of S_c,i over [0, tau];  gain = RMST_i(c) - RMST_i(Stupp)

Three registered properties this module enforces rather than documents:

* **Emax is never silently zero.** `ParameterSet.emax` raises while a payload's value waits
  on the search of item 11 (i.4). A zero enters only through `null_payload=True`, which is
  item 11(o)'s control and says so in the returned record.
* **"combo" is per carrier** (item 11(j)). The pair comes from the parameter file, which is
  checked against item 4's exclusions at load time, so a combo can never carry a payload its
  carrier is registered as unable to carry.
* **Toxicity is always present** (item 11(k)). It is the only source of harm, so T5's harm
  row could not fail without it; there is no flag to switch it off.

Constants live in `configs/response_surface_v1.yaml` and nowhere else (item 11(n)).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from gliostrata.simulation.parameters import ParameterError, ParameterSet

SINGLE_PAYLOADS = ("CRAd", "prodrug_enzyme", "CAR")


@dataclass(frozen=True)
class EffectRecord:
    """One configuration's log-hazard effect, with every term kept for reporting."""

    payload: str
    carrier: str
    benefit: float            # E, the log-hazard reduction (positive is good)
    toxicity: float           # T, the excess log hazard (positive is harm)
    log_hazard_offset: float  # -E + T, what the curve is raised to
    terms: dict[str, float] = field(default_factory=dict)
    null_payload: bool = False
    uncalibrated: bool = True  # item 11(h); never set False by this module


def patient_modifier(payload: str, modifiers: dict[str, float], params: ParameterSet) -> float:
    """m_p(i) = exp(sum of w_p * u_i) over the modifiers that apply to this payload.

    A modifier the payload does not take contributes nothing, by a registered slope of zero
    (for example the immune-TME term does not apply to prodrug_enzyme). An unknown modifier
    name is an error rather than a silent zero.
    """
    total = 0.0
    for name, value in modifiers.items():
        slope = params.modifier_slope(name, payload)  # raises on an unknown name
        total += slope * float(value)
    return float(np.exp(total))


def _single_payload_effect(
    payload: str,
    n_eff: float,
    modifier: float,
    sigma: float,
    cov: float,
    params: ParameterSet,
    *,
    null_payload: bool,
) -> float:
    emax = params.null_payload_emax() if null_payload else params.emax(payload)
    hill = params.hill_h
    ed50 = params.ed50_cells
    if n_eff < 0:
        raise ValueError("n_eff cannot be negative")
    occupancy = 0.0 if n_eff == 0 else n_eff**hill / (n_eff**hill + ed50**hill)
    return float(emax * modifier * sigma * cov * occupancy)


def payload_benefit(
    payload: str,
    carrier: str,
    n_eff: float,
    modifiers: dict[str, float],
    sigma: float,
    cov: float,
    params: ParameterSet,
    *,
    null_payload: bool = False,
) -> tuple[float, dict[str, float]]:
    """E, the log-hazard reduction, with its per-payload parts.

    For "combo" the carrier's registered pair combine by Bliss independence on the hazard
    scale, (E_a + E_b) * synergy, with the synergy factor from the parameter file.
    """
    if payload == "combo":
        pair = params.combo_pair(carrier)
        parts = {
            p: _single_payload_effect(
                p, n_eff, patient_modifier(p, modifiers, params), sigma, cov, params,
                null_payload=null_payload)
            for p in pair
        }
        total = sum(parts.values()) * params.synergy
        return float(total), {**parts, "synergy": params.synergy, "pair": pair}
    if payload not in SINGLE_PAYLOADS:
        raise ParameterError(f"unknown payload {payload!r}")
    modifier = patient_modifier(payload, modifiers, params)
    value = _single_payload_effect(
        payload, n_eff, modifier, sigma, cov, params, null_payload=null_payload)
    return value, {payload: value, "modifier": modifier}


def toxicity(
    payload: str,
    carrier: str,
    dose_cells: float,
    n_injection_sites: int,
    params: ParameterSet,
) -> float:
    """T, the excess log hazard. Item 11(k): the only source of harm, and never switched off."""
    if dose_cells < 0:
        raise ValueError("dose_cells cannot be negative")
    if n_injection_sites < 1:
        raise ValueError("n_injection_sites starts at 1")
    tau_p = params.tau_payload(payload, carrier)
    dose_term = (dose_cells / params.dose_reference_cells) ** params.dose_exponent_q
    site_term = 1.0 + params.site_penalty_eta * (n_injection_sites - 1)
    return float(tau_p * dose_term * site_term + params.delay_cost_log_hr)


def configuration_effect(
    payload: str,
    carrier: str,
    dose_cells: float,
    n_injection_sites: int,
    n_eff: float,
    modifiers: dict[str, float],
    sigma: float,
    cov: float,
    params: ParameterSet,
    *,
    null_payload: bool = False,
) -> EffectRecord:
    """Assemble E and T into the offset the baseline curve is raised to."""
    benefit, parts = payload_benefit(
        payload, carrier, n_eff, modifiers, sigma, cov, params, null_payload=null_payload)
    harm = toxicity(payload, carrier, dose_cells, n_injection_sites, params)
    return EffectRecord(
        payload=payload, carrier=carrier, benefit=benefit, toxicity=harm,
        log_hazard_offset=-benefit + harm,
        terms={**parts, "toxicity": harm, "n_eff": n_eff, "coverage": cov,
               "schedule_factor": sigma},
        null_payload=null_payload,
    )


def configuration_curve(baseline_survival: np.ndarray, log_hazard_offset: float) -> np.ndarray:
    """S_c(t) = S_stupp(t) ** exp(offset). A proportional-hazards transform of the curve."""
    surv = np.asarray(baseline_survival, dtype=float)
    if np.any(surv < 0) or np.any(surv > 1):
        raise ValueError("a survival curve must lie in [0, 1]")
    return surv ** float(np.exp(log_hazard_offset))


def rmst(survival: np.ndarray, params: ParameterSet) -> float:
    """Restricted mean survival time to tau, by the trapezoid rule on the registered grid.

    `survival` must be sampled on `time_grid(params)`: the same fixed grid for every
    configuration and every patient, so that two RMSTs are always comparable.
    """
    surv = np.asarray(survival, dtype=float)
    grid = time_grid(params)
    if surv.shape[-1] != grid.size:
        raise ValueError(
            f"survival has {surv.shape[-1]} points but the registered grid has {grid.size}; "
            "RMST is only comparable on the one grid (item 11(b))")
    return float(np.trapezoid(surv, grid)) if hasattr(np, "trapezoid") else float(
        np.trapz(surv, grid))


def time_grid(params: ParameterSet) -> np.ndarray:
    """The fixed 0 to tau grid every curve is evaluated on (item 11(b))."""
    step = params.grid_step_months
    n_points = int(round(params.tau_months / step)) + 1
    return np.linspace(0.0, params.tau_months, n_points)


def gain_months(
    baseline_survival: np.ndarray,
    effect: EffectRecord,
    params: ParameterSet,
) -> float:
    """gain = RMST(configuration) - RMST(Stupp), in months (item 11(b) and (c))."""
    baseline_rmst = rmst(baseline_survival, params)
    treated = configuration_curve(baseline_survival, effect.log_hazard_offset)
    return rmst(treated, params) - baseline_rmst
