"""L1 — clinical harmonisation with explicit column ownership.

WHY THIS EXISTS
    The L0 clinical cohort (``clinical_cohort.parquet``) was given its schema
    before any measurement existed, so it carries all-null placeholder columns
    under the names that L2 and the molecular call tables produce (``subtype``,
    ``state_*``, ``tme_*``, ``mgmt_methylated``, ``idh_mutated``). A placeholder
    under such a name can shadow a computed value downstream. This module emits
    the clinical table WITHOUT those names, and it never drops data silently: a
    column under an L2 or molecular name that carries any value stops it.

WHAT IT EMITS
    Two tables, each keyed by ``patient_id``, one row per cohort patient, in
    cohort order:

    clinical  ``patient_id`` and the clinical-record columns in
              :data:`CLINICAL_COLUMNS`, plus one ``reason__<column>`` string
              column per clinical column. The reason is null where the value
              is present and a code from
              :data:`gliostrata.features.contract.MISSING_REASONS` where it is
              not.
    labels    ``patient_id`` and the outcome columns in :data:`LABEL_COLUMNS`,
              with their own ``reason__<column>`` columns. Outcomes are labels,
              not features: they live in this separate table and never in the
              clinical table.

HOW A MISSING VALUE GETS ITS REASON
    With the raw GDC records (``clinical_raw.parquet``) supplied, each missing
    value is traced to its source field:

    ``not_recorded_at_source``     the raw field is empty (null, blank, or an
                                   empty treatment list ``[]``);
    ``recorded_unknown_at_source`` the raw field holds an explicit "unknown"
                                   token ("Unknown", "Not Reported", ...);
    ``lost_in_harmonisation``      the raw field holds some other value and
                                   the harmonised value is still missing;
    ``not_available_in_cohort``    the field does not exist in the raw records
                                   at all (extent of resection in TCGA-GBM);
    ``missing_in_source_table``    no raw records were supplied, so the
                                   missing value cannot be traced further.

    Nothing is imputed and no value is changed: values are copied from the L0
    cohort as they are.
"""
from __future__ import annotations

import json

import pandas as pd

#: Clinical-record columns this layer owns, in output order.
CLINICAL_COLUMNS: tuple[str, ...] = (
    "age", "sex", "kps", "treatment_arm", "extent_of_resection",
)

#: Outcome columns. Labels, not features; they go to the label table only.
LABEL_COLUMNS: tuple[str, ...] = ("overall_survival_months", "vital_status")

#: Names produced by L2 profiling or by the molecular call tables. The clinical
#: table must never carry them.
L2_AND_MOLECULAR_NAMES: tuple[str, ...] = (
    "subtype",
    "state_npc", "state_opc", "state_ac", "state_mes",
    "tme_immune", "tme_stromal", "tme_endothelial",
    "mgmt_methylated", "idh_mutated", "mgmt_score",
)

REASON_PREFIX = "reason__"

#: For each harmonised column, the raw GDC field(s) it is derived from.
#: ``None`` means the raw records have no such field.
RAW_SOURCE_FIELDS: dict[str, tuple[str, ...] | None] = {
    "age": ("age_at_index_raw",),
    "sex": ("sex_at_birth_raw",),
    "kps": ("kps_raw",),
    "treatment_arm": ("treatment_types",),
    "extent_of_resection": None,
    "overall_survival_months": ("days_to_death", "days_to_follow_up"),
    "vital_status": ("vital_status_raw",),
}

#: Raw tokens that record "unknown" rather than leave the field empty.
UNKNOWN_TOKENS = frozenset({
    "unknown", "not reported", "not available", "not applicable", "na", "n/a",
    "[not available]", "[unknown]", "[not reported]",
})


def reason_column(column: str) -> str:
    return f"{REASON_PREFIX}{column}"


def _raw_state(value) -> str:
    """Classify one raw field value: 'empty', 'unknown' or 'value'."""
    if value is None:
        return "empty"
    try:
        if pd.isna(value):
            return "empty"
    except (TypeError, ValueError):
        pass
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return "empty"
        if text.lower() in UNKNOWN_TOKENS:
            return "unknown"
        if text.startswith("["):
            try:
                parsed = json.loads(text)
            except (json.JSONDecodeError, ValueError):
                return "value"
            if isinstance(parsed, list) and not [p for p in parsed if p]:
                return "empty"
    return "value"


def _trace_reason(column: str, raw_row: pd.Series | None, raw_columns) -> str:
    if raw_row is None:
        return "missing_in_source_table"
    fields = RAW_SOURCE_FIELDS.get(column)
    if fields is None or not all(f in raw_columns for f in fields):
        return "not_available_in_cohort"
    states = [_raw_state(raw_row[f]) for f in fields]
    if all(s == "empty" for s in states):
        return "not_recorded_at_source"
    if all(s in ("empty", "unknown") for s in states):
        return "recorded_unknown_at_source"
    return "lost_in_harmonisation"


def harmonise_clinical(
    cohort: pd.DataFrame,
    raw: pd.DataFrame | None = None,
    *,
    raw_id_column: str = "submitter_id",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split the L0 cohort into an owned clinical table and a label table.

    Raises ``ValueError`` if a column under an L2 or molecular name carries a
    value, if the cohort has a column this layer does not know, if a required
    column is missing, or if ``patient_id`` is not unique.
    """
    if "patient_id" not in cohort.columns:
        raise ValueError("cohort: no patient_id column")
    ids = cohort["patient_id"].astype(str)
    if ids.duplicated().any():
        raise ValueError("cohort: duplicate patient_id values")

    known = {"patient_id", *CLINICAL_COLUMNS, *LABEL_COLUMNS, *L2_AND_MOLECULAR_NAMES}
    unknown = [c for c in cohort.columns if c not in known]
    if unknown:
        raise ValueError(
            f"cohort: column(s) {unknown} are not owned by the clinical layer and "
            "are not known L2 or molecular names; refusing to drop them silently"
        )
    carrying = [c for c in cohort.columns
                if c in L2_AND_MOLECULAR_NAMES and cohort[c].notna().any()]
    if carrying:
        raise ValueError(
            f"cohort: column(s) {carrying} use an L2 or molecular name and carry "
            "values. The clinical table must not supply them; remove them at L0 "
            "or route them through the input that owns them"
        )
    absent = [c for c in (*CLINICAL_COLUMNS, *LABEL_COLUMNS) if c not in cohort.columns]
    if absent:
        raise ValueError(f"cohort: required column(s) missing: {absent}")

    raw_by_id = None
    if raw is not None:
        if raw_id_column not in raw.columns:
            raise ValueError(f"raw: no {raw_id_column} column")
        raw_ids = raw[raw_id_column].astype(str)
        if raw_ids.duplicated().any():
            raise ValueError(f"raw: duplicate {raw_id_column} values")
        raw_by_id = raw.set_index(raw_ids)
        not_in_raw = sorted(set(ids) - set(raw_by_id.index))
        if not_in_raw:
            raise ValueError(f"raw: {len(not_in_raw)} cohort patient(s) have no raw record")

    def _table(columns: tuple[str, ...]) -> pd.DataFrame:
        out = pd.DataFrame({"patient_id": ids.to_numpy()})
        for col in columns:
            out[col] = cohort[col].to_numpy()
        for col in columns:
            reasons = []
            for pid, value in zip(ids, cohort[col], strict=True):
                if not pd.isna(value):
                    reasons.append(None)
                    continue
                row = None if raw_by_id is None else raw_by_id.loc[pid]
                reasons.append(_trace_reason(col, row, raw_by_id.columns if raw_by_id is not None else ()))
            out[reason_column(col)] = pd.array(reasons, dtype="string")
        return out

    clinical = _table(CLINICAL_COLUMNS)
    labels = _table(LABEL_COLUMNS)
    # Fix dtypes so that the written files do not depend on how the cohort
    # happened to encode an all-null column.
    clinical["age"] = pd.to_numeric(clinical["age"]).astype("float64")
    clinical["kps"] = pd.to_numeric(clinical["kps"]).astype("float64")
    for col in ("sex", "treatment_arm", "extent_of_resection"):
        clinical[col] = clinical[col].astype("string")
    labels["overall_survival_months"] = pd.to_numeric(
        labels["overall_survival_months"]).astype("float64")
    labels["vital_status"] = pd.to_numeric(labels["vital_status"]).astype("Int8")
    clinical["patient_id"] = clinical["patient_id"].astype("string")
    labels["patient_id"] = labels["patient_id"].astype("string")
    return clinical, labels


def reason_counts(table: pd.DataFrame) -> dict[str, dict[str, int]]:
    """``{column: {reason: count}}`` for every ``reason__`` column with a reason."""
    out: dict[str, dict[str, int]] = {}
    for col in table.columns:
        if not col.startswith(REASON_PREFIX):
            continue
        counts = table[col].dropna().value_counts()
        if len(counts):
            out[col[len(REASON_PREFIX):]] = {str(k): int(v) for k, v in sorted(counts.items())}
    return out
