"""Model-ready populations: named views over the feature store (decision D2).

The store keeps every patient. A view declares the columns it requires and
excludes, step by step, the patients who lack one, counting each step. Its
rules come from availability alone, fixed before any outcome is looked at; no
view compares outcomes.

Three views are defined:

``profiling``
    Patients with an expression sample. The population of the L2 profiles.
``transferable_model``
    Age, sex and treatment arm present, and an observed outcome (both labels).
    The feature set of T2's transferable model (variant 2).
``molecular_model``
    The transferable requirements plus an IDH call and an MGMT call, each with
    call status ``called``. The feature set of T2's molecular model (variant 3).

A missing outcome label is an exclusion reason of a model view, never a
feature. KPS is not required by any view: it is missing for 152 of the 594
TCGA patients and is not available in CGGA at all; a view that needs it is
declared separately, as complete cases (decision D4).

Each step is counted sequentially: ``n_excluded`` is the number of patients
still in the view who lack that requirement. ``n_missing_overall`` is the
number of store patients who lack it, whatever the earlier steps did. Counts
are also split by the frozen development partition when one is given.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from gliostrata.features.contract import LABEL_COLUMNS, FeatureStore


@dataclass(frozen=True)
class Requirement:
    """One exclusion step: patients lacking ``column`` (as defined) leave."""

    step: str
    table: str          # "features" or "labels"
    column: str
    rule: str           # "present", "true" or "called"
    reason: str


@dataclass(frozen=True)
class View:
    name: str
    description: str
    requirements: tuple[Requirement, ...]
    feature_columns: tuple[str, ...]
    label_columns: tuple[str, ...] = ()

    def required_columns(self) -> list[str]:
        return [r.column for r in self.requirements]


def _present(step, col, reason, table="features"):
    return Requirement(step, table, col, "present", reason)


_TRANSFERABLE_REQ = (
    _present("age recorded", "age", "no age"),
    _present("sex recorded", "sex", "no sex"),
    _present("treatment arm recorded", "treatment_arm", "no treatment arm"),
    _present("survival time recorded", "overall_survival_months",
             "no survival time (n_missing_outcome)", "labels"),
    _present("vital status recorded", "vital_status",
             "no vital status (n_missing_outcome)", "labels"),
)

VIEWS: dict[str, View] = {
    "profiling": View(
        name="profiling",
        description="cohort patients with an expression sample",
        requirements=(
            Requirement("expression sample", "features", "expression_available",
                        "true", "no expression sample"),
        ),
        feature_columns=(),  # filled from the store: the L2 columns present
    ),
    "transferable_model": View(
        name="transferable_model",
        description="age, sex and treatment arm present, and an observed outcome",
        requirements=_TRANSFERABLE_REQ,
        feature_columns=("age", "sex", "treatment_arm"),
        label_columns=LABEL_COLUMNS,
    ),
    "molecular_model": View(
        name="molecular_model",
        description="transferable requirements plus an IDH call and an MGMT call",
        requirements=_TRANSFERABLE_REQ + (
            Requirement("IDH called", "features", "idh_call_status", "called",
                        "no IDH call (not tested or not callable)"),
            Requirement("MGMT called", "features", "mgmt_call_status", "called",
                        "no MGMT call (not tested or not callable)"),
        ),
        feature_columns=("age", "sex", "treatment_arm", "idh_mutated", "mgmt_methylated"),
        label_columns=LABEL_COLUMNS,
    ),
}

L2_PREFIXES = ("subtype", "state_", "tme_")


@dataclass
class ViewResult:
    view: View
    patient_ids: list[str]
    features: pd.DataFrame
    labels: pd.DataFrame | None
    flow: list[dict]

    def flow_table(self) -> pd.DataFrame:
        return pd.DataFrame(self.flow)


def _meets(store: FeatureStore, req: Requirement) -> pd.Series:
    frame = store.labels if req.table == "labels" else store.features
    col = frame[req.column]
    if req.rule == "present":
        ok = col.notna()
    elif req.rule == "true":
        ok = col.fillna(False).astype(bool)
    elif req.rule == "called":
        ok = col.astype("string").eq("called").fillna(False)
    else:
        raise ValueError(f"unknown rule {req.rule!r}")
    return pd.Series(ok.to_numpy(dtype=bool), index=store.features["patient_id"].astype(str))


def _split_counts(ids, split_of: dict[str, str] | None) -> dict[str, int]:
    if split_of is None:
        return {}
    counts = {"train": 0, "test": 0, "not_in_partition": 0}
    for pid in ids:
        counts[split_of.get(pid, "not_in_partition")] += 1
    return counts


def partition_lookup(train_ids, test_ids) -> dict[str, str]:
    lookup = {str(p): "train" for p in train_ids}
    for p in test_ids:
        if str(p) in lookup:
            raise ValueError(f"patient {p} is in both splits")
        lookup[str(p)] = "test"
    return lookup


def build_view(store: FeatureStore, view: View | str,
               split_of: dict[str, str] | None = None) -> ViewResult:
    """Apply a view's requirements to the store and count every step."""
    view = VIEWS[view] if isinstance(view, str) else view
    ids = store.features["patient_id"].astype(str)
    keep = pd.Series(True, index=ids.to_numpy())
    flow = [{"step": "store", "reason": "all store patients",
             "n_excluded": 0, "n_missing_overall": 0, "n_remaining": int(keep.sum()),
             **{f"n_remaining_{k}": v for k, v in _split_counts(ids, split_of).items()},
             **{f"n_excluded_{k}": 0 for k in _split_counts([], split_of)}}]
    for req in view.requirements:
        ok = _meets(store, req)
        excluded = keep & ~ok
        keep = keep & ok
        ex_ids = [p for p, e in excluded.items() if e]
        rem_ids = [p for p, k in keep.items() if k]
        flow.append({
            "step": req.step, "reason": req.reason,
            "n_excluded": len(ex_ids),
            "n_missing_overall": int((~ok).sum()),
            "n_remaining": len(rem_ids),
            **{f"n_remaining_{k}": v for k, v in _split_counts(rem_ids, split_of).items()},
            **{f"n_excluded_{k}": v for k, v in _split_counts(ex_ids, split_of).items()},
        })
    kept = [p for p, k in keep.items() if k]
    if view.name == "profiling":
        cols = [c for c in store.features.columns
                if c.startswith(L2_PREFIXES) and store.columns[c].role == "feature"]
    else:
        cols = list(view.feature_columns)
    feats = store.features.set_index(store.features["patient_id"].astype(str))
    features = feats.loc[kept, ["patient_id", *cols]].reset_index(drop=True)
    if view.name == "profiling":
        gaps = {c: int(n) for c, n in features[cols].isna().sum().items() if n}
        if gaps:
            raise ValueError(f"profiling view: L2 scoring failed for some patients {gaps}")
    labels = None
    if view.label_columns:
        lab = store.labels.set_index(store.labels["patient_id"].astype(str))
        labels = lab.loc[kept, ["patient_id", *view.label_columns]].reset_index(drop=True)
    return ViewResult(view=view, patient_ids=kept, features=features,
                      labels=labels, flow=flow)


def build_all_views(store: FeatureStore,
                    split_of: dict[str, str] | None = None) -> dict[str, ViewResult]:
    return {name: build_view(store, view, split_of) for name, view in VIEWS.items()}
