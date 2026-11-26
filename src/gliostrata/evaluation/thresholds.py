"""Pre-registered evaluation thresholds (Protocol Table 8-2, tests T1-T5).

This module is the *machine-readable pre-registration* of the GlioStrata
computational success criteria. It exists so that no threshold can be
negotiated after the fact: the numbers live in version control, they are
hashed, and the hash is asserted by the test suite.

Three guarantees the module is designed to provide
---------------------------------------------------
1. **Immutability.** Every threshold is a frozen dataclass and the registry is
   exposed through a read-only mapping, so an accidental in-place edit at
   runtime raises instead of silently changing the bar.
2. **Loud failure on edits.** :data:`THRESHOLD_TABLE_SHA256` is a SHA-256 over
   the canonical serialisation of the whole table. ``tests/test_evaluation.py``
   asserts it. Changing *any* pre-registered number, direction, CI bound or
   metric key breaks the build and forces the author to update the hash as a
   deliberate, reviewable commit.
3. **No evaluation against unregistered metrics.** :func:`evaluate` raises
   :class:`UnregisteredMetricError` for any metric key that was not
   pre-registered, which blocks the classic p-hacking route of inventing a
   friendlier metric once the results are in.

Decision rule (Protocol 8.2)
----------------------------
PASS requires the relevant bound of the bootstrap 95% confidence interval to
clear the threshold -- the *lower* bound for "higher is better" metrics, the
*upper* bound for the harm metric. The point estimate alone is never
sufficient.

This module adds one refinement the protocol leaves implicit, so that honest
reporting is possible: a result whose point estimate clears the threshold but
whose CI bound does not is reported as ``INCONCLUSIVE`` rather than ``FAIL``.
A ``FAIL`` means the point estimate itself missed. Both are threshold misses
for pre-registration purposes (neither may be reported as a PASS); the
distinction exists purely so the write-up can describe *how* it missed.
"""
from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from enum import Enum
from types import MappingProxyType

__all__ = [
    "CIBound",
    "Direction",
    "Status",
    "Threshold",
    "Verdict",
    "UnregisteredMetricError",
    "THRESHOLDS",
    "THRESHOLD_TABLE_SHA256",
    "canonical_table",
    "threshold_table_checksum",
    "get_threshold",
    "registered_metrics",
    "evaluate",
    "evaluate_all",
    "format_threshold_table",
]


class Direction(str, Enum):
    """Comparison direction of a pre-registered threshold."""

    GTE = "gte"  # metric must be >= threshold ("higher is better")
    LTE = "lte"  # metric must be <= threshold ("lower is better")

    @property
    def symbol(self) -> str:
        return ">=" if self is Direction.GTE else "<="


class CIBound(str, Enum):
    """Which bound of the bootstrap 95% CI must clear the threshold."""

    LOWER = "lower"
    UPPER = "upper"


class Status(str, Enum):
    """Outcome of evaluating a metric against its pre-registered threshold."""

    PASS = "PASS"
    FAIL = "FAIL"
    INCONCLUSIVE = "INCONCLUSIVE"


class UnregisteredMetricError(KeyError):
    """Raised when a caller tries to evaluate a metric that was not pre-registered."""


@dataclass(frozen=True, slots=True)
class Threshold:
    """One row of Protocol Table 8-2.

    Attributes
    ----------
    key:
        Stable identifier used by :func:`evaluate`. Part of the hashed table,
        so renaming a key is a pre-registration amendment.
    test_id:
        Protocol test this threshold belongs to (``T1``..``T5``).
    metric:
        Human-readable metric name as written in Table 8-2.
    direction:
        ``GTE`` if higher is better, ``LTE`` if lower is better.
    value:
        The pre-registered threshold value.
    ci_bound:
        Which bootstrap 95% CI bound must clear ``value`` for a PASS.
    justification:
        The protocol's stated reason for this particular number.
    unit:
        Unit of ``value`` (``"fraction"``, ``"correlation"``, ``"count"``,
        ``"percent_of_final_volume"``, ``"percent_reduction"``, ...).
    point_estimate_limit:
        Optional *additional* requirement on the point estimate itself, used by
        the T5 safety floor where the protocol states two numbers
        (point estimate <= 10%, CI upper bound <= 15%).
    auxiliary_requirement:
        Optional non-CI condition that must also hold for a PASS (T4's
        ``Wilcoxon p < 0.01``). Supplied to :func:`evaluate` as
        ``auxiliary_satisfied``; the verdict cannot be PASS while it is unmet
        or unknown.
    notes:
        Implementation caveats recorded at pre-registration time.
    """

    key: str
    test_id: str
    metric: str
    direction: Direction
    value: float
    ci_bound: CIBound
    justification: str
    unit: str
    point_estimate_limit: float | None = None
    auxiliary_requirement: str | None = None
    notes: str = field(default="")

    def as_row(self) -> dict[str, object]:
        """Return a plain-JSON dict of this threshold (used for hashing)."""
        row = asdict(self)
        row["direction"] = self.direction.value
        row["ci_bound"] = self.ci_bound.value
        return row

    def describe(self) -> str:
        """One-line human description, e.g. ``T2 survival C-index ... >= 0.65 (CI lower)``."""
        base = (
            f"{self.test_id} {self.metric}: {self.direction.symbol} {self.value:g} "
            f"(95% CI {self.ci_bound.value} bound)"
        )
        if self.point_estimate_limit is not None:
            base += f"; point estimate {self.direction.symbol} {self.point_estimate_limit:g}"
        if self.auxiliary_requirement:
            base += f"; also requires {self.auxiliary_requirement}"
        return base


# ---------------------------------------------------------------------------
# Protocol Table 8-2 -- PRE-REGISTERED. Do not edit without amending
# THRESHOLD_TABLE_SHA256 below and documenting the amendment in
# docs/preregistration.md. The test suite will fail loudly otherwise.
# ---------------------------------------------------------------------------

_TABLE: tuple[Threshold, ...] = (
    Threshold(
        key="T1_subtype_accuracy_cgga",
        test_id="T1",
        metric="Subtyping accuracy (CGGA external)",
        direction=Direction.GTE,
        value=0.88,
        ci_bound=CIBound.LOWER,
        justification="Best published external subtyping accuracy is ~0.88-0.90.",
        unit="fraction",
    ),
    Threshold(
        key="T1_state_deconvolution_pearson_r",
        test_id="T1",
        metric="State deconvolution Pearson r",
        direction=Direction.GTE,
        value=0.75,
        ci_bound=CIBound.LOWER,
        justification="CIBERSORTx-class deconvolution typically achieves r 0.70-0.85.",
        unit="correlation",
        notes=(
            "Reference truth for Neftel state fractions must be the "
            "single-cell-derived fractions named in the analysis plan; the "
            "protocol does not state whether r is pooled across states or "
            "averaged per state (see docs/preregistration.md, Ambiguities)."
        ),
    ),
    Threshold(
        key="T2_survival_cindex_cgga",
        test_id="T2",
        metric="Survival C-index (CGGA external)",
        direction=Direction.GTE,
        value=0.65,
        ci_bound=CIBound.LOWER,
        justification="Best published external GBM survival C-index is ~0.65-0.68.",
        unit="cindex",
    ),
    Threshold(
        key="T2_cindex_improvement_vs_baseline",
        test_id="T2",
        metric="Improvement over best-published 2024 baseline",
        direction=Direction.GTE,
        value=0.03,
        ci_bound=CIBound.LOWER,
        justification="Demonstrate improvement, not parity, against the 2024 state of the art.",
        unit="cindex_delta",
        notes=(
            "Requires the comparator model to be re-run on the same frozen "
            "partition; a delta against a number quoted from a paper computed "
            "on a different cohort is not a valid bootstrap CI."
        ),
    ),
    Threshold(
        key="T2_brier_18mo",
        test_id="T2",
        metric="Brier score @ 18 months",
        direction=Direction.LTE,
        value=0.25,
        ci_bound=CIBound.UPPER,
        justification="GBM 18-month Brier baseline is ~0.23-0.27.",
        unit="brier",
    ),
    Threshold(
        key="T3_growth_rmse_pct_final_volume",
        test_id="T3",
        metric="Tumor-growth RMSE vs PDX",
        direction=Direction.LTE,
        value=15.0,
        ci_bound=CIBound.UPPER,
        justification="Standard PDX growth-modeling benchmark tolerance.",
        unit="percent_of_final_volume",
    ),
    Threshold(
        key="T3_nsc_migration_r_in_vivo",
        test_id="T3",
        metric="NSC-migration r vs in vivo",
        direction=Direction.GTE,
        value=0.65,
        ci_bound=CIBound.LOWER,
        justification="In vivo NSC tracking data is noisy, so the bar is set below T1's r.",
        unit="correlation",
    ),
    Threshold(
        key="T4_optimizer_regret_reduction_pct",
        test_id="T4",
        metric="Optimizer beats random search (regret reduction)",
        direction=Direction.GTE,
        value=25.0,
        ci_bound=CIBound.LOWER,
        justification="A statistical test, not a point estimate, must carry the claim.",
        unit="percent_reduction",
        auxiliary_requirement="Wilcoxon signed-rank p < 0.01",
    ),
    Threshold(
        key="T4_personalization_distinct_configs",
        test_id="T4",
        metric="Personalization variance (distinct configs in top-5 across 50 patients)",
        direction=Direction.GTE,
        value=3.0,
        ci_bound=CIBound.LOWER,
        justification="Proves the optimizer has not converged to a universal recommendation.",
        unit="count",
        notes=(
            "A count over one fixed 50-patient cohort has no sampling "
            "distribution unless the cohort is itself resampled; if no "
            "bootstrap is performed the CI is reported as degenerate "
            "(lower = point = upper)."
        ),
    ),
    Threshold(
        key="T5_endtoend_rmst_success_fraction",
        test_id="T5",
        metric="End-to-end success (fraction of held-out patients with >= 3-month RMST gain)",
        direction=Direction.GTE,
        value=0.50,
        ci_bound=CIBound.LOWER,
        justification="Pre-registered primary success criterion from the working hypothesis.",
        unit="fraction",
    ),
    Threshold(
        key="T5_safety_predicted_harm_fraction",
        test_id="T5",
        metric="Safety floor (predicted harm rate)",
        direction=Direction.LTE,
        value=0.15,
        ci_bound=CIBound.UPPER,
        justification="Hard safety stop: harm rate <= 10% with 95% CI upper bound <= 15%.",
        unit="fraction",
        point_estimate_limit=0.10,
        notes=(
            "Two numbers are pre-registered. `value` is the CI-upper-bound "
            "requirement (0.15); `point_estimate_limit` is the point-estimate "
            "requirement (0.10). Both must hold for a PASS."
        ),
    ),
)

THRESHOLDS: Mapping[str, Threshold] = MappingProxyType({t.key: t for t in _TABLE})
"""Read-only registry of the pre-registered thresholds, keyed by ``Threshold.key``."""

if len(THRESHOLDS) != len(_TABLE):  # pragma: no cover - guards a typo in the table above
    raise RuntimeError("Duplicate threshold key in the pre-registered table.")


def canonical_table() -> str:
    """Return the canonical JSON serialisation of the threshold table.

    Keys are sorted and separators fixed so the string -- and therefore the
    checksum -- depends only on the registered content, not on dict ordering,
    whitespace or Python version.
    """
    rows = [t.as_row() for t in _TABLE]
    rows.sort(key=lambda r: str(r["key"]))
    return json.dumps(rows, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def threshold_table_checksum() -> str:
    """SHA-256 hex digest over :func:`canonical_table`."""
    return hashlib.sha256(canonical_table().encode("utf-8")).hexdigest()


THRESHOLD_TABLE_SHA256 = "96bbcd49dcda20f8a0bc9c8b71427aaa24cf754dc6d8d7b4f7d6a2df3d2874ab"
"""Frozen checksum of the pre-registered threshold table.

If a test fails on this constant, a pre-registered number changed. That is not
a test to "fix": either revert the edit, or amend the pre-registration
deliberately (update this constant, record the amendment and its date in
``docs/preregistration.md``, and disclose it in the paper).
"""


def registered_metrics() -> tuple[str, ...]:
    """Return all registered metric keys, sorted."""
    return tuple(sorted(THRESHOLDS))


def get_threshold(metric_key: str) -> Threshold:
    """Look up a pre-registered threshold.

    Raises
    ------
    UnregisteredMetricError
        If ``metric_key`` was not pre-registered.
    """
    try:
        return THRESHOLDS[metric_key]
    except KeyError:
        raise UnregisteredMetricError(
            f"'{metric_key}' is not a pre-registered metric. Evaluating against an "
            f"unregistered metric would defeat the pre-registration, so this is refused. "
            f"Registered metrics: {', '.join(registered_metrics())}."
        ) from None


@dataclass(frozen=True, slots=True)
class Verdict:
    """Structured outcome of one threshold evaluation."""

    metric_key: str
    test_id: str
    metric: str
    status: Status
    point_estimate: float
    ci_lower: float
    ci_upper: float
    threshold: float
    direction: Direction
    ci_bound: CIBound
    decisive_value: float
    reason: str
    justification: str
    auxiliary_requirement: str | None = None
    auxiliary_satisfied: bool | None = None

    @property
    def passed(self) -> bool:
        return self.status is Status.PASS

    def to_dict(self) -> dict[str, object]:
        row = asdict(self)
        row["status"] = self.status.value
        row["direction"] = self.direction.value
        row["ci_bound"] = self.ci_bound.value
        return row

    def __str__(self) -> str:  # pragma: no cover - convenience
        return f"[{self.status.value}] {self.test_id} {self.metric} -- {self.reason}"


def _clears(direction: Direction, observed: float, limit: float) -> bool:
    return observed >= limit if direction is Direction.GTE else observed <= limit


def evaluate(
    metric_key: str,
    point_estimate: float,
    ci_lower: float,
    ci_upper: float,
    *,
    auxiliary_satisfied: bool | None = None,
) -> Verdict:
    """Evaluate a ``(point_estimate, ci_lower, ci_upper)`` triple against a threshold.

    Parameters
    ----------
    metric_key:
        A key from :func:`registered_metrics`. Anything else raises
        :class:`UnregisteredMetricError`.
    point_estimate, ci_lower, ci_upper:
        The observed metric and its bootstrap 95% confidence interval. Must
        satisfy ``ci_lower <= point_estimate <= ci_upper``.
    auxiliary_satisfied:
        Whether the threshold's :attr:`Threshold.auxiliary_requirement` holds.
        Required (as ``True``) for a PASS on thresholds that declare one.

    Returns
    -------
    Verdict
        Including the ``reason`` string explaining *why* the verdict is what it
        is -- the decisive comparison, in words.

    Raises
    ------
    UnregisteredMetricError
        For an unregistered ``metric_key``.
    ValueError
        If the CI is malformed, or if any input is not finite.
    """
    spec = get_threshold(metric_key)

    for name, val in (
        ("point_estimate", point_estimate),
        ("ci_lower", ci_lower),
        ("ci_upper", ci_upper),
    ):
        fval = float(val)
        if fval != fval or fval in (float("inf"), float("-inf")):
            raise ValueError(f"{name} must be a finite number, got {val!r}.")

    point_estimate, ci_lower, ci_upper = float(point_estimate), float(ci_lower), float(ci_upper)
    if not (ci_lower <= point_estimate <= ci_upper):
        raise ValueError(
            f"Malformed confidence interval for '{metric_key}': requires "
            f"ci_lower <= point_estimate <= ci_upper, got "
            f"({ci_lower}, {point_estimate}, {ci_upper})."
        )

    decisive = ci_lower if spec.ci_bound is CIBound.LOWER else ci_upper
    bound_clears = _clears(spec.direction, decisive, spec.value)
    point_clears = _clears(spec.direction, point_estimate, spec.value)

    sym = spec.direction.symbol
    bound_phrase = (
        f"95% CI {spec.ci_bound.value} bound {decisive:g} {sym} {spec.value:g}"
        if bound_clears
        else f"95% CI {spec.ci_bound.value} bound {decisive:g} does not clear {sym} {spec.value:g}"
    )

    # Secondary point-estimate limit (T5 safety floor).
    extra_ok = True
    extra_phrase = ""
    if spec.point_estimate_limit is not None:
        extra_ok = _clears(spec.direction, point_estimate, spec.point_estimate_limit)
        extra_phrase = (
            f"; point estimate {point_estimate:g} "
            f"{'satisfies' if extra_ok else 'violates'} the pre-registered limit "
            f"{sym} {spec.point_estimate_limit:g}"
        )

    # Auxiliary non-CI requirement (T4 Wilcoxon test).
    aux_ok = True
    aux_unknown = False
    aux_phrase = ""
    if spec.auxiliary_requirement is not None:
        aux_ok = auxiliary_satisfied is True
        aux_unknown = auxiliary_satisfied is None
        if aux_unknown:
            aux_phrase = (
                f"; auxiliary requirement '{spec.auxiliary_requirement}' was not "
                f"reported, so a PASS cannot be claimed"
            )
        else:
            aux_phrase = (
                f"; auxiliary requirement '{spec.auxiliary_requirement}' "
                f"{'holds' if aux_ok else 'does not hold'}"
            )

    tail = f"{extra_phrase}{aux_phrase}."
    if not extra_ok:
        # A violated point-estimate limit is a hard miss (T5 safety floor).
        status = Status.FAIL
        reason = f"FAIL: pre-registered point-estimate limit violated; {bound_phrase}{tail}"
    elif spec.auxiliary_requirement is not None and not aux_ok and not aux_unknown:
        status = Status.FAIL
        reason = f"FAIL: pre-registered auxiliary requirement not met; {bound_phrase}{tail}"
    elif bound_clears and aux_ok:
        status = Status.PASS
        reason = f"PASS: {bound_phrase}{tail}"
    elif not bound_clears and point_clears:
        # Point estimate on the right side of the line but the CI does not
        # settle it -- an honest "inconclusive", never reportable as a PASS.
        status = Status.INCONCLUSIVE
        reason = (
            f"INCONCLUSIVE: point estimate {point_estimate:g} {sym} {spec.value:g} but "
            f"{bound_phrase}{tail}"
        )
    elif bound_clears and aux_unknown:
        status = Status.INCONCLUSIVE
        reason = f"INCONCLUSIVE: {bound_phrase}{tail}"
    else:
        status = Status.FAIL
        reason = (
            f"FAIL: point estimate {point_estimate:g} does not clear {sym} "
            f"{spec.value:g}; {bound_phrase}{tail}"
        )

    return Verdict(
        metric_key=spec.key,
        test_id=spec.test_id,
        metric=spec.metric,
        status=status,
        point_estimate=point_estimate,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        threshold=spec.value,
        direction=spec.direction,
        ci_bound=spec.ci_bound,
        decisive_value=decisive,
        reason=reason,
        justification=spec.justification,
        auxiliary_requirement=spec.auxiliary_requirement,
        auxiliary_satisfied=auxiliary_satisfied,
    )


def evaluate_all(
    results: Mapping[str, tuple[float, float, float]],
    *,
    auxiliary: Mapping[str, bool] | None = None,
) -> dict[str, Verdict]:
    """Evaluate several metrics at once.

    ``results`` maps metric key -> ``(point_estimate, ci_lower, ci_upper)``.
    ``auxiliary`` optionally maps metric key -> whether its auxiliary
    requirement holds. Unregistered keys raise.
    """
    aux = auxiliary or {}
    return {
        key: evaluate(key, *triple, auxiliary_satisfied=aux.get(key))
        for key, triple in results.items()
    }


def format_threshold_table() -> str:
    """Render the pre-registered table as Markdown (used by docs/preregistration.md)."""
    lines = [
        "| Test | Metric | Requirement | Decisive bound | Justification |",
        "| --- | --- | --- | --- | --- |",
    ]
    for t in _TABLE:
        req = f"{t.direction.symbol} {t.value:g}"
        if t.point_estimate_limit is not None:
            req += f" (CI) and {t.direction.symbol} {t.point_estimate_limit:g} (point)"
        if t.auxiliary_requirement:
            req += f" and {t.auxiliary_requirement}"
        lines.append(
            f"| {t.test_id} | {t.metric} | {req} | 95% CI {t.ci_bound.value} | {t.justification} |"
        )
    return "\n".join(lines)
