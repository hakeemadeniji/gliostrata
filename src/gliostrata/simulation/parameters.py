"""Load and validate the T4/T5 response-surface parameters.

Pre-registration section 2.2 item 11(n) registers `configs/response_surface_v1.yaml` as the
single place every numeric parameter of the response surface, the delivery module and the
toxicity term lives. This module reads that file and holds **no constants of its own**, so
that the registered checksum cannot be broken by an edit that changes no parameter.

Two guards matter more than the loading:

* ``ParameterSet.require_frozen()`` refuses to proceed while the file says ``frozen: false``.
  Item 11(e) and (m) allow a T4 or T5 read only after the freeze amendment is entered, and a
  pending Emax value must not reach a registered number by accident.
* ``ParameterSet.emax(payload)`` raises if that payload's value is still
  ``pending_registered_search``. The only way to get a zero out of it is to ask for the
  null-payload control of item 11(o) explicitly, by name.

Nothing here reads a cohort, so no read-guard access is recorded.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

DEFAULT_PATH = Path("configs/response_surface_v1.yaml")
PENDING = "pending_registered_search"
CARRIERS = ("iPSC-NSC", "MSC", "HB1_F3_line")
SINGLE_PAYLOADS = ("CRAd", "prodrug_enzyme", "CAR")
# Item 4's registered incompatibilities, repeated here only to check the combo pairs
# against them. `gliostrata.optimization.config_space` remains authoritative.
INCOMPATIBLE = frozenset({("HB1_F3_line", "CAR"), ("MSC", "CRAd")})


class ParameterError(RuntimeError):
    """The parameter file is missing, inconsistent, or not usable for what was asked."""


def _value(entry: Any, where: str) -> Any:
    """Return the number an entry carries, whether it is wrapped or bare."""
    if isinstance(entry, dict):
        if "value" not in entry:
            raise ParameterError(f"{where}: a mapping with no 'value' key")
        return entry["value"]
    return entry


@dataclass(frozen=True)
class ParameterSet:
    """The parsed parameter file, with the guards item 11 requires."""

    doc: dict[str, Any]
    path: Path
    sha256: str

    # -- guards ---------------------------------------------------------------------------

    @property
    def frozen(self) -> bool:
        return bool(self.doc.get("frozen", False))

    def require_frozen(self, what: str) -> None:
        """Refuse a registered computation while the parameter file is still a draft."""
        if not self.frozen:
            raise ParameterError(
                f"{what} needs frozen parameters. {self.path} says frozen: false, so the "
                "freeze amendment of pre-registration item 11(e) and (m) has not been "
                "entered. Dry runs and the controls of item 11(o) may proceed; a registered "
                "T4 or T5 read may not.")

    def pending(self) -> list[str]:
        """Payloads whose Emax is still waiting on the registered search of item 11(i.4)."""
        emax = self.doc["response"]["emax"]
        return [p for p in SINGLE_PAYLOADS
                if isinstance(emax.get(p), dict) and emax[p].get("status") == PENDING]

    # -- horizon --------------------------------------------------------------------------

    @property
    def tau_months(self) -> float:
        return float(_value(self.doc["horizon"]["tau_months"], "horizon.tau_months"))

    @property
    def grid_step_months(self) -> float:
        return float(_value(self.doc["horizon"]["grid_step_months"], "horizon.grid_step"))

    # -- delivery -------------------------------------------------------------------------

    def tropism_genes(self, carrier: str) -> tuple[str, ...]:
        genes = self.doc["delivery"]["tropism_genes"]
        if carrier not in genes:
            raise ParameterError(f"no tropism genes registered for carrier {carrier!r}")
        return tuple(genes[carrier])

    def cell_fraction_acting(self, carrier: str) -> float:
        table = self.doc["delivery"]["cell_fraction_acting"]
        if carrier not in table:
            raise ParameterError(f"no acting fraction registered for carrier {carrier!r}")
        return float(_value(table[carrier], f"cell_fraction_acting.{carrier}"))

    @property
    def tropism_logit_intercept(self) -> float:
        return float(_value(self.doc["delivery"]["tropism_logit_intercept"], "a_k"))

    @property
    def tropism_logit_slope(self) -> float:
        return float(_value(self.doc["delivery"]["tropism_logit_slope_per_sd"], "g_k"))

    @property
    def coverage_saturation_sites(self) -> float:
        return float(_value(self.doc["delivery"]["coverage_saturation_sites"], "n0"))

    def schedule_factor(self, schedule: str) -> float:
        table = self.doc["delivery"]["schedule_factor"]
        if schedule not in table:
            raise ParameterError(f"no schedule factor registered for {schedule!r}")
        return float(_value(table[schedule], f"schedule_factor.{schedule}"))

    # -- response -------------------------------------------------------------------------

    @property
    def ed50_cells(self) -> float:
        return float(_value(self.doc["response"]["dose_response"]["ed50_cells"], "ed50"))

    @property
    def hill_h(self) -> float:
        return float(_value(self.doc["response"]["dose_response"]["hill_h"], "hill_h"))

    @property
    def translation_factor(self) -> float:
        return float(_value(
            self.doc["response"]["emax"]["translation_factor"], "translation_factor"))

    def emax(self, payload: str) -> float:
        """The registered Emax for a single payload. Raises while the value is pending."""
        if payload not in SINGLE_PAYLOADS:
            raise ParameterError(
                f"emax() takes a single payload {SINGLE_PAYLOADS}, not {payload!r}; "
                "'combo' is assembled from its carrier's pair by combo_pair()")
        entry = self.doc["response"]["emax"][payload]
        value = _value(entry, f"emax.{payload}")
        if value is None or (isinstance(entry, dict) and entry.get("status") == PENDING):
            raise ParameterError(
                f"Emax for {payload!r} is still {PENDING}: the registered search of "
                "pre-registration item 11 (i.4) has run, but no ratio has been entered for "
                "this payload yet. Enter it from docs/t4t5_emax_search/ and nowhere else. "
                "For the null-payload control of item 11(o), call null_payload_emax().")
        return float(value)

    @staticmethod
    def null_payload_emax() -> float:
        """Item 11(o)'s null-payload control: Emax = 0, which must give a gain of exactly 0."""
        return 0.0

    def combo_pair(self, carrier: str) -> tuple[str, str]:
        """Item 11(j): the two payloads this carrier is registered as able to carry."""
        pairs = self.doc["response"]["combo"]["pairs"]
        if carrier not in pairs:
            raise ParameterError(f"no combo pair registered for carrier {carrier!r}")
        pair = tuple(pairs[carrier])
        if len(pair) != 2:
            raise ParameterError(f"combo pair for {carrier!r} is not a pair: {pair}")
        return pair  # type: ignore[return-value]

    @property
    def synergy(self) -> float:
        return float(_value(self.doc["response"]["combo"]["synergy"], "synergy"))

    def modifier_slope(self, name: str, payload: str) -> float:
        """A patient-modifier slope, zero for a payload the modifier does not apply to."""
        table = self.doc["response"]["patient_modifier"]
        if name not in table:
            raise ParameterError(f"no patient modifier registered named {name!r}")
        entry = table[name]
        if payload not in entry.get("applies_to", []):
            return 0.0
        return float(_value(entry, f"patient_modifier.{name}"))

    @property
    def modifier_names(self) -> tuple[str, ...]:
        """The actual modifiers, not every key under `patient_modifier`.

        A modifier is a mapping carrying a `value`. Annotation keys sitting alongside them
        (a note, a pointer to the reference they standardise against) are not modifiers, and
        returning them would have made any note look like a term in the surface.
        """
        table = self.doc["response"]["patient_modifier"]
        return tuple(name for name, entry in table.items()
                     if isinstance(entry, dict) and "value" in entry)

    # -- toxicity -------------------------------------------------------------------------

    @property
    def dose_reference_cells(self) -> float:
        return float(_value(self.doc["toxicity"]["dose_reference_cells"], "dose_reference"))

    def tau_payload(self, payload: str, carrier: str) -> float:
        """Item 11(k) and (j): combo's toxicity is the sum of its carrier's two payloads."""
        table = self.doc["toxicity"]["tau_payload"]
        if payload == "combo":
            if table.get("combo") != "sum_of_pair":
                raise ParameterError(
                    "toxicity.tau_payload.combo must be 'sum_of_pair' (item 11(j)); "
                    f"found {table.get('combo')!r}")
            return sum(self.tau_payload(p, carrier) for p in self.combo_pair(carrier))
        if payload not in table:
            raise ParameterError(f"no toxicity coefficient registered for {payload!r}")
        return float(_value(table[payload], f"tau_payload.{payload}"))

    @property
    def dose_exponent_q(self) -> float:
        return float(_value(self.doc["toxicity"]["dose_exponent_q"], "q"))

    @property
    def site_penalty_eta(self) -> float:
        return float(_value(self.doc["toxicity"]["site_penalty_eta"], "eta"))

    @property
    def delay_cost_log_hr(self) -> float:
        return float(_value(self.doc["toxicity"]["delay_cost_log_hr"], "delta_delay"))


def _check(params: ParameterSet) -> None:
    """Internal consistency checks that must hold whether or not the file is frozen."""
    if params.tau_months <= 0 or params.grid_step_months <= 0:
        raise ParameterError("tau and the grid step must both be positive")
    if params.grid_step_months > params.tau_months:
        raise ParameterError("the grid step is larger than the horizon")
    if params.doc["survival_model"]["variant"] != 2:
        raise ParameterError(
            "item 11(g) registers variant 2 (age, sex, arm); the file names variant "
            f"{params.doc['survival_model']['variant']}")
    for banned in ("IDH", "MGMT", "idh", "mgmt"):
        if banned in params.doc["survival_model"].get("covariates", []):
            raise ParameterError(
                f"item 11(g) and (l) put IDH and MGMT in no part of the model; found {banned}")
    if params.ed50_cells <= 0 or params.hill_h <= 0:
        raise ParameterError("ED50 and the Hill exponent must both be positive")
    if params.coverage_saturation_sites <= 0:
        raise ParameterError("the coverage saturation constant must be positive")

    # Item 11(j): every carrier has a pair, of two payloads it may actually carry.
    for carrier in CARRIERS:
        pair = params.combo_pair(carrier)
        if len(set(pair)) != 2:
            raise ParameterError(f"{carrier}'s combo pair repeats a payload: {pair}")
        for payload in pair:
            if payload not in SINGLE_PAYLOADS:
                raise ParameterError(f"{carrier}'s combo pair names unknown payload {payload!r}")
            if (carrier, payload) in INCOMPATIBLE:
                raise ParameterError(
                    f"{carrier}'s combo pair contains {payload!r}, which item 4 registers as "
                    "incompatible with it. This is the contradiction item 11(j) exists to "
                    "avoid; a uniform combo definition would reintroduce it.")
    for carrier in params.doc["response"]["combo"]["pairs"]:
        if carrier not in CARRIERS:
            raise ParameterError(f"combo pair registered for unknown carrier {carrier!r}")

    # The tropism reference must never be a CGGA statistic (delivery.tropism_reference).
    source = str(params.doc["delivery"]["tropism_reference"].get("source", "")).lower()
    if "cgga" in source:
        raise ParameterError(
            "the tropism reference names CGGA. Item 11 and section 3.1 require a frozen "
            "TCGA reference, so that the sealed cohort cannot leak into the transform.")


def load_parameters(path: str | Path = DEFAULT_PATH) -> ParameterSet:
    """Read, hash and validate the registered parameter file."""
    path = Path(path)
    try:
        raw = path.read_bytes()
    except FileNotFoundError as exc:
        raise ParameterError(
            f"{path} not found. Item 11(n) registers it as the only place the T4/T5 "
            "parameters live; the module carries no defaults to fall back on.") from exc
    doc = yaml.safe_load(raw.decode("utf-8"))
    if not isinstance(doc, dict):
        raise ParameterError(f"{path} did not parse as a mapping")
    for section in ("horizon", "survival_model", "delivery", "response", "toxicity"):
        if section not in doc:
            raise ParameterError(f"{path} has no {section!r} section")
    params = ParameterSet(doc=doc, path=path, sha256=hashlib.sha256(raw).hexdigest())
    _check(params)
    return params
