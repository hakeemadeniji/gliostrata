"""L3 — Unified patient feature store.

``build_feature_store`` is the single code path that assembles the store, in
both modes:

* ``mode="real"`` — TCGA inputs: the L1 clinical and label tables
  (:mod:`gliostrata.data.clinical`), the expression-derived L2 outputs
  (single-sample subtype, state and TME scores) and the molecular calls.
* ``mode="synthetic"`` — the offline simulator and the L2 stubs. Labelled as
  synthetic in every file it writes. This is the default for tests and for the
  offline pipeline.

Each input is taken by explicit name and owns a fixed set of columns. The
build raises, before assembling anything, on:

(a) a column supplied by more than one input (:class:`DuplicateSourceError`);
(b) a null-only column that shadows a computed feature, i.e. an all-null
    column under a name another input computes (:class:`ShadowColumnError`);
(c) a column an input does not own, an outcome label among them
    (:class:`ColumnOwnershipError`);
(d) in real mode, a stub or cohort-relative subtype output
    (:class:`StubInRealModeError`).

The store keeps every patient of the clinical table and never filters;
model-ready populations are the named views of :mod:`gliostrata.features.views`.
See :mod:`gliostrata.features.contract` for the files, the column classes and
the missing-value reason codes.

The L5 to L7 layers are still stubs that expect the earlier one-hot layout.
:func:`stub_layer_frame` derives that layout from a store, and refuses to do so
where it would have to turn an absent value into a value.
"""
from __future__ import annotations

from collections.abc import Iterable, Sequence

import pandas as pd

from gliostrata.data.clinical import CLINICAL_COLUMNS, REASON_PREFIX
from gliostrata.features.contract import (
    CALL_STATUSES,
    LABEL_COLUMNS,
    MODES,
    ColumnOwnershipError,
    ColumnSpec,
    DuplicateSourceError,
    FeatureStore,
    FeatureStoreError,
    LabelInFeatureMatrixError,
    ShadowColumnError,
    StubInRealModeError,
    UnexplainedMissingError,
)

STATE_COLUMNS: tuple[str, ...] = ("state_npc", "state_opc", "state_ac", "state_mes")
TME_COLUMNS: tuple[str, ...] = ("tme_immune", "tme_stromal", "tme_endothelial")
MOLECULAR_COLUMNS: tuple[str, ...] = (
    "mgmt_score", "mgmt_threshold", "mgmt_methylated", "mgmt_call_status",
    "idh_mutated", "idh_call_status",
)
SUBTYPE_FIXED: tuple[str, ...] = ("subtype", "subtype_n_genes_used")
SUBTYPE_SINGLE_SAMPLE_PREFIX = "subtype_corr_"
SUBTYPE_COHORT_RELATIVE_PREFIX = "subtype_prob_"

#: Metadata columns each L2 input may carry. They are read, checked for being
#: constant, recorded in the column metadata and not copied into the store.
INPUT_METADATA: dict[str, tuple[str, ...]] = {
    "subtypes": (
        "subtype_mode", "subtype_centroid_set", "subtype_centroid_sha256",
        "subtype_classifier_version", "subtype_scoring_method",
        "subtype_is_supervised", "subtype_scaling",
    ),
    "states": ("state_deconvolution_method", "state_deconvolution_version",
               "is_true_deconvolution"),
    "tme": ("tme_deconvolution_method", "tme_deconvolution_version",
            "is_true_deconvolution"),
}

#: Columns computed by an input other than the clinical table. A null-only
#: column under one of these names, supplied by another input, is a shadow.
COMPUTED_FEATURES: frozenset[str] = frozenset({
    *SUBTYPE_FIXED, *STATE_COLUMNS, *TME_COLUMNS, *MOLECULAR_COLUMNS,
})

#: Legacy layout consumed by the L5 to L7 stubs (survival, optimiser, handoff).
FEATURE_COLUMNS = [
    "age", "kps",
    "sex_male",
    "extent_of_resection_STR", "extent_of_resection_GTR",
    "mgmt_methylated", "idh_mutated",
    "subtype_prob_Classical", "subtype_prob_Mesenchymal",
    "subtype_prob_Proneural", "subtype_prob_Neural",
    "state_npc", "state_opc", "state_ac", "state_mes",
    "tme_immune", "tme_stromal", "tme_endothelial",
]

_STUB_MARKERS = ("stub_passthrough",)


def _owned(input_name: str, column: str) -> bool:
    if column == "patient_id":
        return True
    if column.startswith(REASON_PREFIX):
        return input_name in ("clinical", "labels")
    if input_name == "clinical":
        return column in CLINICAL_COLUMNS
    if input_name == "labels":
        return column in LABEL_COLUMNS
    if input_name == "subtypes":
        return (column in SUBTYPE_FIXED
                or column.startswith(SUBTYPE_SINGLE_SAMPLE_PREFIX)
                or column.startswith(SUBTYPE_COHORT_RELATIVE_PREFIX))
    if input_name == "states":
        return column in STATE_COLUMNS
    if input_name == "tme":
        return column in TME_COLUMNS
    if input_name == "molecular":
        return column in MOLECULAR_COLUMNS
    return False


def _computed_by_other(input_name: str, column: str) -> bool:
    if column in ("patient_id",) or column.startswith(REASON_PREFIX):
        return False
    if _owned(input_name, column):
        return False
    return (column in COMPUTED_FEATURES
            or column.startswith(SUBTYPE_SINGLE_SAMPLE_PREFIX)
            or column.startswith(SUBTYPE_COHORT_RELATIVE_PREFIX))


def _ids(frame: pd.DataFrame, label: str) -> pd.Series:
    if "patient_id" not in frame.columns:
        raise FeatureStoreError(f"{label}: no patient_id column")
    ids = frame["patient_id"].astype(str)
    if ids.duplicated().any():
        raise FeatureStoreError(f"{label}: duplicate patient_id values")
    return ids


def _split_metadata(name: str, frame: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    meta_cols = [c for c in INPUT_METADATA.get(name, ()) if c in frame.columns]
    meta: dict = {}
    for col in meta_cols:
        values = frame[col].dropna().unique().tolist()
        if len(values) > 1:
            raise FeatureStoreError(f"{name}: metadata column {col} is not constant")
        value = values[0] if values else None
        if hasattr(value, "item"):
            value = value.item()
        meta[col] = value
    return frame.drop(columns=meta_cols), meta


def check_inputs(inputs: dict[str, pd.DataFrame]) -> None:
    """Rules (b) and (a) of the module docstring, in that order."""
    for name, frame in inputs.items():
        for col in frame.columns:
            if _computed_by_other(name, col) and frame[col].isna().all():
                raise ShadowColumnError(
                    f"input {name!r} supplies {col!r} with no value in any row; "
                    f"{col!r} is a computed feature owned by another input, and a "
                    "null-only column under its name would shadow it"
                )
    suppliers: dict[str, list[str]] = {}
    for name, frame in inputs.items():
        for col in frame.columns:
            if col != "patient_id":
                suppliers.setdefault(col, []).append(name)
    dupes = {c: s for c, s in suppliers.items() if len(s) > 1}
    if dupes:
        raise DuplicateSourceError(
            "column(s) supplied by more than one input: "
            + "; ".join(f"{c} by {', '.join(s)}" for c, s in sorted(dupes.items()))
        )


def check_ownership(inputs: dict[str, pd.DataFrame]) -> None:
    """Rule (c): every column belongs to the input that supplies it."""
    for name, frame in inputs.items():
        foreign = [c for c in frame.columns if not _owned(name, c)]
        if foreign:
            labels = [c for c in foreign if c in LABEL_COLUMNS]
            hint = (f"; {labels} are outcome labels and belong only in the label table"
                    if labels else "")
            raise ColumnOwnershipError(
                f"input {name!r} supplies column(s) it does not own: {foreign}{hint}"
            )


def _check_real_mode(subtypes: pd.DataFrame, meta: dict[str, dict]) -> None:
    sub_meta = meta.get("subtypes", {})
    mode = sub_meta.get("subtype_mode")
    if mode != "single_sample":
        raise StubInRealModeError(
            f"real mode needs single-sample subtype calls (subtype_mode "
            f"'single_sample'); got subtype_mode={mode!r}. Cohort-relative and "
            "stub subtype outputs are refused for per-patient features"
        )
    prob = [c for c in subtypes.columns if c.startswith(SUBTYPE_COHORT_RELATIVE_PREFIX)]
    if prob:
        raise StubInRealModeError(
            f"real mode refuses cohort-relative subtype probabilities {prob}"
        )
    for name, values in meta.items():
        for col, value in values.items():
            text = str(value)
            if text in _STUB_MARKERS or text.endswith("-stub"):
                raise StubInRealModeError(
                    f"input {name!r} is a stub output ({col}={value!r}); a stub "
                    "consumed in real mode is refused"
                )


def _reason_series(frame: pd.DataFrame, col: str) -> pd.Series | None:
    rc = f"{REASON_PREFIX}{col}"
    return frame[rc] if rc in frame.columns else None


def _reasons_from_table(table: pd.DataFrame, columns: Sequence[str]) -> list[tuple]:
    rows = []
    for col in columns:
        reasons = _reason_series(table, col)
        null = table[col].isna()
        for pid, is_null, reason in zip(
            table["patient_id"].astype(str), null,
            reasons if reasons is not None else [None] * len(table), strict=True,
        ):
            if not is_null:
                if reasons is not None and not pd.isna(reason):
                    raise FeatureStoreError(
                        f"{col}: patient {pid} has a value and a missing reason")
                continue
            rows.append((pid, col, None if pd.isna(reason) else str(reason)))
    return rows


def build_feature_store(
    *,
    clinical: pd.DataFrame,
    labels: pd.DataFrame,
    subtypes: pd.DataFrame,
    states: pd.DataFrame,
    tme: pd.DataFrame,
    molecular: pd.DataFrame,
    expression_samples: Iterable[str],
    mode: str = "synthetic",
    sources: dict[str, list[dict]] | None = None,
    scoring: dict[str, str] | None = None,
) -> FeatureStore:
    """Assemble the store. Every argument is required and named.

    ``expression_samples`` is the list of sample ids in the expression input
    (in synthetic mode, the simulator's patients). ``sources`` maps an input
    name to its files (``[{"path", "sha256"}]``); ``scoring`` maps an input
    name to a note on how its values were scored.
    """
    if mode not in MODES:
        raise FeatureStoreError(f"mode must be one of {MODES}, got {mode!r}")
    sources = sources or {}
    scoring = dict(scoring or {})

    frames: dict[str, pd.DataFrame] = {}
    meta: dict[str, dict] = {}
    for name, frame in (("clinical", clinical), ("labels", labels),
                        ("subtypes", subtypes), ("states", states), ("tme", tme),
                        ("molecular", molecular)):
        if not isinstance(frame, pd.DataFrame):
            raise FeatureStoreError(f"input {name!r} must be a DataFrame")
        _ids(frame, name)
        frames[name], meta[name] = _split_metadata(name, frame.copy())
    check_inputs(frames)
    if mode == "real":
        _check_real_mode(frames["subtypes"], meta)
    check_ownership(frames)
    if "subtype_mode" in meta["subtypes"]:
        scoring.setdefault("subtypes", str(meta["subtypes"]["subtype_mode"]))

    clin = frames["clinical"].copy()
    clin["patient_id"] = clin["patient_id"].astype(str)
    ids = clin["patient_id"].tolist()
    id_set = set(ids)
    missing_cols = [c for c in CLINICAL_COLUMNS if c not in clin.columns]
    if missing_cols:
        raise FeatureStoreError(f"clinical: required column(s) missing: {missing_cols}")

    lab = frames["labels"].copy()
    lab["patient_id"] = lab["patient_id"].astype(str)
    if set(lab["patient_id"]) != id_set:
        raise FeatureStoreError("labels: patients differ from the clinical table")
    missing_cols = [c for c in LABEL_COLUMNS if c not in lab.columns]
    if missing_cols:
        raise FeatureStoreError(f"labels: required column(s) missing: {missing_cols}")
    lab = lab.set_index("patient_id").loc[ids].reset_index()

    samples = [str(s) for s in expression_samples]
    if len(set(samples)) != len(samples):
        raise FeatureStoreError("expression: duplicate sample ids")
    has_expr = [pid in set(samples) for pid in ids]
    profiled = [pid for pid, h in zip(ids, has_expr, strict=True) if h]

    l2: dict[str, pd.DataFrame] = {}
    for name in ("subtypes", "states", "tme"):
        frame = frames[name].copy()
        frame["patient_id"] = frame["patient_id"].astype(str)
        if set(frame["patient_id"]) != set(profiled):
            raise FeatureStoreError(
                f"{name}: patients differ from the cohort patients with an "
                f"expression sample ({len(frame)} rows against {len(profiled)})"
            )
        l2[name] = frame.set_index("patient_id")

    mol = frames["molecular"].copy()
    mol["patient_id"] = mol["patient_id"].astype(str)
    if set(mol["patient_id"]) != id_set:
        raise FeatureStoreError("molecular: patients differ from the clinical table")
    missing_cols = [c for c in MOLECULAR_COLUMNS if c not in mol.columns]
    if missing_cols:
        raise FeatureStoreError(f"molecular: required column(s) missing: {missing_cols}")
    mol = mol.set_index("patient_id").loc[ids].reset_index()
    for status_col, value_col in (("mgmt_call_status", "mgmt_methylated"),
                                  ("idh_call_status", "idh_mutated")):
        bad = sorted(set(mol[status_col].dropna().astype(str)) - set(CALL_STATUSES))
        if bad or mol[status_col].isna().any():
            raise FeatureStoreError(f"{status_col}: values outside {CALL_STATUSES}: {bad}")
        called = mol[status_col].astype(str).eq("called")
        if (called != mol[value_col].notna()).any():
            raise FeatureStoreError(
                f"{value_col}: a call exists exactly when {status_col} is 'called'; "
                "the input breaks that"
            )

    # --- assemble the features table, in clinical order --------------------
    feats = pd.DataFrame({"patient_id": pd.array(ids, dtype="string")})
    for col in CLINICAL_COLUMNS:
        feats[col] = clin[col].to_numpy()
    feats["expression_available"] = pd.array(has_expr, dtype="boolean")

    sub_cols = [c for c in l2["subtypes"].columns]
    # Score columns keep the classifier's subtype order.
    ordered_sub = (["subtype"]
                   + [c for c in sub_cols if c.startswith(SUBTYPE_SINGLE_SAMPLE_PREFIX)]
                   + [c for c in sub_cols if c.startswith(SUBTYPE_COHORT_RELATIVE_PREFIX)]
                   + [c for c in ("subtype_n_genes_used",) if c in sub_cols])
    if "subtype" not in sub_cols:
        raise FeatureStoreError("subtypes: no subtype column")
    l2_columns: dict[str, list[str]] = {
        "subtypes": ordered_sub,
        "states": [c for c in STATE_COLUMNS],
        "tme": [c for c in TME_COLUMNS],
    }
    for name, cols in l2_columns.items():
        absent = [c for c in cols if c not in l2[name].columns]
        if absent:
            raise FeatureStoreError(f"{name}: required column(s) missing: {absent}")
        block = l2[name].reindex(ids)[cols]
        for col in cols:
            feats[col] = block[col].to_numpy()
    for col in MOLECULAR_COLUMNS:
        feats[col] = mol[col].to_numpy()

    # dtypes: nullable throughout, so absent stays absent
    for col in ("age", "kps", "mgmt_score", "mgmt_threshold", *STATE_COLUMNS, *TME_COLUMNS,
                *[c for c in ordered_sub if c.startswith(("subtype_corr_", "subtype_prob_"))]):
        feats[col] = pd.to_numeric(feats[col]).astype("Float64")
    for col in ("sex", "treatment_arm", "extent_of_resection", "subtype",
                "mgmt_call_status", "idh_call_status"):
        feats[col] = feats[col].astype("string")
    if "subtype_n_genes_used" in feats.columns:
        feats["subtype_n_genes_used"] = pd.to_numeric(
            feats["subtype_n_genes_used"]).astype("Int64")
    for col in ("mgmt_methylated", "idh_mutated"):
        values = pd.to_numeric(feats[col].astype("Float64"))
        bad = values.dropna()
        if (~bad.isin([0.0, 1.0])).any():
            raise FeatureStoreError(f"{col}: values other than 0, 1 or missing")
        feats[col] = values.map(lambda v: pd.NA if pd.isna(v) else bool(v)).astype("boolean")

    labels_out = pd.DataFrame({"patient_id": pd.array(ids, dtype="string")})
    labels_out["overall_survival_months"] = pd.to_numeric(
        lab["overall_survival_months"]).astype("Float64").to_numpy()
    labels_out["vital_status"] = pd.array(
        pd.to_numeric(lab["vital_status"]).astype("Int8"), dtype="Int8")

    # --- one reason per null cell --------------------------------------------
    rows: list[tuple] = []
    rows += _reasons_from_table(clin, CLINICAL_COLUMNS)
    rows += _reasons_from_table(lab, LABEL_COLUMNS)
    expr_by_id = dict(zip(ids, has_expr, strict=True))
    for cols in l2_columns.values():
        for col in cols:
            for pid, is_null in zip(ids, feats[col].isna(), strict=True):
                if is_null:
                    rows.append((pid, col, "scoring_failed" if expr_by_id[pid]
                                 else "no_expression_sample"))
    status_of = {
        "mgmt_score": "mgmt_call_status", "mgmt_threshold": "mgmt_call_status",
        "mgmt_methylated": "mgmt_call_status", "idh_mutated": "idh_call_status",
    }
    for col, status_col in status_of.items():
        for pid, is_null, status in zip(ids, feats[col].isna(), feats[status_col], strict=True):
            if not is_null:
                continue
            rows.append((pid, col, None if status == "called" else str(status)))
    # A column absent for every patient with no traced reason is not available
    # in this cohort; any other untraced null is an error.
    resolved = []
    for pid, col, reason in rows:
        if reason is None:
            frame = labels_out if col in LABEL_COLUMNS else feats
            if frame[col].isna().all():
                reason = "not_available_in_cohort"
            else:
                raise UnexplainedMissingError(
                    f"{col}: patient {pid} has no value and no reason code")
        resolved.append((pid, col, reason))
    missing = pd.DataFrame(resolved, columns=["patient_id", "column", "reason"])
    col_order = {c: i for i, c in enumerate([*feats.columns, *labels_out.columns])}
    pos = {pid: i for i, pid in enumerate(ids)}
    missing = (missing.assign(_c=missing["column"].map(col_order),
                              _p=missing["patient_id"].map(pos))
               .sort_values(["_c", "_p"]).drop(columns=["_c", "_p"])
               .reset_index(drop=True))
    for col in missing.columns:
        missing[col] = missing[col].astype("string")
    n_null = int(feats.drop(columns=["patient_id"]).isna().sum().sum()
                 + labels_out.drop(columns=["patient_id"]).isna().sum().sum())
    if n_null != len(missing):
        raise UnexplainedMissingError(
            f"{n_null} null cells but {len(missing)} reason codes")

    # --- column metadata -----------------------------------------------------
    layer = "synthetic:" if mode == "synthetic" else ""
    specs: dict[str, ColumnSpec] = {}

    def add(name, role, src_layer, src_input, klass, desc, score=None):
        frame = labels_out if role == "label" else feats
        if role != "id" and frame[name].isna().all():
            klass = "not_available"
        specs[name] = ColumnSpec(
            name=name, role=role, source_layer=f"{layer}{src_layer}",
            source_input=src_input, klass=klass, description=desc,
            scoring=score, source_files=tuple(sources.get(src_input, ())),
        )

    add("patient_id", "id", "L0", "clinical", "clinical_record", "patient identifier")
    clinical_desc = {
        "age": "age at index, years", "sex": "sex at birth, M or F",
        "kps": "Karnofsky performance status as recorded; never imputed",
        "treatment_arm": "ChemoRadiation, MonoTherapy or NoAdjuvant, from the "
                         "recorded treatment types",
        "extent_of_resection": "GTR, STR or biopsy as recorded; missing is "
                               "not_available, never biopsy",
    }
    for col in CLINICAL_COLUMNS:
        add(col, "feature", "L1", "clinical", "clinical_record", clinical_desc[col])
    add("expression_available", "status", "L0", "expression", "measured",
        "the patient has an expression sample in the expression input")
    sub_score = scoring.get("subtypes")
    for col in ordered_sub:
        role = "status" if col == "subtype_n_genes_used" else "feature"
        add(col, role, "L2", "subtypes", "derived_from_measurement",
            {"subtype": "Verhaak subtype call",
             "subtype_n_genes_used": "signature genes found in the expression input"}
            .get(col, "per-subtype score of the subtype classifier"), sub_score)
    for col in STATE_COLUMNS:
        add(col, "feature", "L2", "states", "derived_from_measurement",
            "Neftel state marker-enrichment score, normalised to sum to 1; not a "
            "cell proportion", scoring.get("states"))
    for col in TME_COLUMNS:
        add(col, "feature", "L2", "tme", "derived_from_measurement",
            "TME marker-enrichment score, normalised to sum to 1; not a cell "
            "proportion", scoring.get("tme"))
    add("mgmt_score", "feature", "L0", "molecular", "measured",
        "MGMT promoter score: mean beta value of the two probes")
    add("mgmt_threshold", "status", "L0", "molecular", "derived_from_measurement",
        "threshold applied to mgmt_score to make the call")
    add("mgmt_methylated", "feature", "L0", "molecular", "derived_from_measurement",
        "mgmt_score >= mgmt_threshold; missing where there is no call, never 0")
    add("mgmt_call_status", "status", "L0", "molecular", "derived_from_measurement",
        "called, not_tested or tested_not_callable")
    add("idh_mutated", "feature", "L0", "molecular", "measured",
        "IDH1 R132 or IDH2 R172 hotspot variant; missing where not sequenced, "
        "never 0")
    add("idh_call_status", "status", "L0", "molecular", "measured",
        "called, not_tested or tested_not_callable")
    add("overall_survival_months", "label", "L1", "labels", "clinical_record",
        "outcome label; never a feature")
    add("vital_status", "label", "L1", "labels", "clinical_record",
        "outcome label (1 dead, 0 alive); never a feature")

    store_meta = {
        "input_metadata": {k: v for k, v in meta.items() if v},
        "scoring": scoring,
        "n_expression_samples": len(samples),
        "n_expression_samples_outside_cohort": len(set(samples) - id_set),
    }
    return FeatureStore(mode=mode, features=feats, labels=labels_out,
                        missing=missing, columns=specs, meta=store_meta)


# ---------------------------------------------------------------------------
# Feature matrices
# ---------------------------------------------------------------------------
def feature_matrix(features: pd.DataFrame | FeatureStore,
                   columns: Sequence[str]) -> pd.DataFrame:
    """The ``patient x feature`` matrix for ``columns``, complete cases only.

    Raises :class:`LabelInFeatureMatrixError` if a column is an outcome label
    (or a reason column), and :class:`FeatureStoreError` if a requested column
    has a missing value: model-ready matrices are built from a view, never by
    filling gaps.
    """
    frame = features.features if isinstance(features, FeatureStore) else features
    cols = list(columns)
    labels = [c for c in cols if c in LABEL_COLUMNS or c.startswith(REASON_PREFIX)]
    if labels:
        raise LabelInFeatureMatrixError(
            f"{labels} are outcome labels or reason codes and are forbidden in a "
            "feature matrix"
        )
    present_labels = [c for c in frame.columns if c in LABEL_COLUMNS]
    if present_labels and not cols:
        raise LabelInFeatureMatrixError(f"frame carries label column(s) {present_labels}")
    absent = [c for c in cols if c not in frame.columns]
    if absent:
        raise FeatureStoreError(f"feature column(s) not in the frame: {absent}")
    matrix = frame[cols].copy()
    gaps = {c: int(n) for c, n in matrix.isna().sum().items() if n}
    if gaps:
        raise FeatureStoreError(
            f"missing values in a feature matrix: {gaps}. Build it from a view "
            "that requires these columns"
        )
    return matrix


# ---------------------------------------------------------------------------
# Synthetic mode
# ---------------------------------------------------------------------------
def synthetic_inputs(cohort: pd.DataFrame) -> dict:
    """Split a synthetic (simulator) cohort into the named store inputs.

    The simulator emits every field in one table, the L2 "truth" columns
    included, because the L2 stubs read them back. Here the clinical, label and
    molecular inputs are cut from it by ownership; the L2 values come only from
    the stub outputs passed to :func:`build_feature_store`.
    """
    ids = _ids(cohort, "synthetic cohort")
    clinical = pd.DataFrame({"patient_id": ids.to_numpy()})
    for col in CLINICAL_COLUMNS:
        clinical[col] = cohort[col].to_numpy() if col in cohort.columns else pd.NA
    labels = pd.DataFrame({"patient_id": ids.to_numpy()})
    for col in LABEL_COLUMNS:
        labels[col] = cohort[col].to_numpy()
    molecular = pd.DataFrame({"patient_id": ids.to_numpy()})
    for kind, col in (("mgmt", "mgmt_methylated"), ("idh", "idh_mutated")):
        values = cohort[col] if col in cohort.columns else pd.Series(pd.NA, index=cohort.index)
        molecular[col] = values.to_numpy()
        molecular[f"{kind}_call_status"] = [
            "not_tested" if pd.isna(v) else "called" for v in values]
    molecular["mgmt_score"] = pd.NA
    molecular["mgmt_threshold"] = pd.NA
    return {"clinical": clinical, "labels": labels, "molecular": molecular,
            "expression_samples": ids.tolist()}


def build_synthetic_store(cohort: pd.DataFrame, subtypes: pd.DataFrame,
                          states: pd.DataFrame, tme: pd.DataFrame,
                          *, sources: dict[str, list[dict]] | None = None) -> FeatureStore:
    """The synthetic-mode store: simulator inputs plus the three L2 stub outputs."""
    return build_feature_store(
        **synthetic_inputs(cohort), subtypes=subtypes, states=states, tme=tme,
        mode="synthetic", sources=sources,
        scoring={"subtypes": "synthetic stub", "states": "synthetic stub",
                 "tme": "synthetic stub"},
    )


# ---------------------------------------------------------------------------
# The layout the L5 to L7 stubs consume
# ---------------------------------------------------------------------------
def stub_layer_frame(store: FeatureStore | pd.DataFrame) -> pd.DataFrame:
    """One-hot layout for the L5 to L7 stubs, derived from a store.

    Refuses (``FeatureStoreError``) where it would have to invent a value: a
    missing resection would become "biopsy" and a missing call 0 in the one-hot
    encoding, so a store with such gaps (the real TCGA store) cannot feed the
    stub layers. A frame already in this layout is returned unchanged.
    """
    frame = store.features if isinstance(store, FeatureStore) else store
    if all(c in frame.columns for c in FEATURE_COLUMNS) and "sex_male" in frame.columns:
        return frame
    needed = ["age", "kps", "sex", "extent_of_resection", "mgmt_methylated",
              "idh_mutated", *[f"subtype_prob_{s}" for s in
                               ("Classical", "Mesenchymal", "Proneural", "Neural")],
              *STATE_COLUMNS, *TME_COLUMNS]
    absent = [c for c in needed if c not in frame.columns]
    if absent:
        raise FeatureStoreError(
            f"the L5 to L7 stubs need {absent}, which this store does not carry "
            "(the real store carries single-sample subtype scores, not "
            "cohort-relative probabilities)"
        )
    gaps = {c: int(n) for c, n in frame[needed].isna().sum().items() if n}
    if gaps:
        raise FeatureStoreError(
            f"the L5 to L7 stubs would encode absent values as values {gaps}; "
            "refusing. Use a named view (gliostrata.features.views)"
        )
    out = pd.DataFrame({"patient_id": frame["patient_id"].astype(str).to_numpy()})
    out["age"] = frame["age"].astype(float).to_numpy()
    out["kps"] = frame["kps"].astype(float).to_numpy()
    out["sex_male"] = (frame["sex"].astype(str) == "M").astype(int).to_numpy()
    for level in ("STR", "GTR"):
        out[f"extent_of_resection_{level}"] = (
            frame["extent_of_resection"].astype(str) == level).astype(int).to_numpy()
    for col in ("mgmt_methylated", "idh_mutated"):
        out[col] = frame[col].astype(bool).astype(int).to_numpy()
    for col in FEATURE_COLUMNS:
        if col not in out.columns:
            out[col] = frame[col].astype(float).to_numpy()
    out = out[["patient_id", *FEATURE_COLUMNS]]
    out["treatment_arm"] = frame["treatment_arm"].astype(object).to_numpy()
    out["subtype"] = frame["subtype"].astype(object).to_numpy()
    out.attrs["feature_columns"] = FEATURE_COLUMNS
    out.attrs["version"] = "stub-layer-1"
    return out


def synthetic_feature_frame(cohort: pd.DataFrame, subtypes: pd.DataFrame,
                            states: pd.DataFrame, tme: pd.DataFrame) -> pd.DataFrame:
    """Synthetic store, then the stub-layer layout: what L5 to L7 consume offline."""
    return stub_layer_frame(build_synthetic_store(cohort, subtypes, states, tme))
