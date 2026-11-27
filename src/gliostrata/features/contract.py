"""L3 feature-store contract: vocabularies, errors, the store object, its files.

THE CONTRACT (decisions D1 to D5, docs/project_record/05_future_plan.md 5.4)
    * The store holds every cohort patient, one row each, and never filters.
      Model-ready populations are named views built from it
      (:mod:`gliostrata.features.views`).
    * Every column carries a role, a source layer, the input it came from, the
      SHA-256 of the file(s) behind that input, and a class from
      :data:`COLUMN_CLASSES`.
    * Every missing value carries exactly one reason code from
      :data:`MISSING_REASONS`. Absent never becomes a value: no imputation, no
      0 for a missing call, no "biopsy" for a missing resection.
    * Outcomes are labels, not features. They live in a separate label table
      keyed by ``patient_id`` and are refused in any feature matrix.

THE FILES (written by :func:`write_store`, all deterministic)
    ``<store>.parquet``          features table: ``patient_id`` and every
                                 feature and status column, one row per
                                 patient, nullable dtypes, nulls where absent.
    ``<store>_labels.parquet``   label table: ``patient_id`` and the outcome
                                 columns (role ``label``).
    ``<store>_missing.parquet``  long table ``patient_id, column, reason``: one
                                 row per null cell of the two tables above,
                                 sorted by column then patient.
    ``<store>_columns.json``     store-level metadata (mode, version, the
                                 vocabularies) and one record per column:
                                 role, source layer, source input, source
                                 files with SHA-256, class, dtype, missing
                                 count, scoring note, description.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

STORE_VERSION = "1.0.0"

#: The mode a store was built in. ``synthetic`` stores come from the offline
#: simulator and the L2 stubs and are labelled as such in every file.
MODES: tuple[str, ...] = ("synthetic", "real")

#: Column roles. ``label`` columns live only in the label table.
ROLES: tuple[str, ...] = ("id", "feature", "status", "label")

#: Column classes (decision D1).
COLUMN_CLASSES: dict[str, str] = {
    "measured": "a value read from an assay (expression availability, a probe "
                "beta value, a sequenced variant)",
    "derived_from_measurement": "computed from a measurement by a stated rule "
                                "(a thresholded MGMT call, an expression-derived "
                                "subtype, state or TME score)",
    "clinical_record": "copied from the clinical record (age, sex, KPS, "
                       "treatment arm, outcome)",
    "not_available": "the source provides no value for any patient of this cohort",
}

#: Missing-value reason codes (decision D1). Exactly one per null cell.
MISSING_REASONS: dict[str, str] = {
    "not_recorded_at_source": "the raw source record has no entry for this field",
    "recorded_unknown_at_source": "the raw source record holds an explicit "
                                  "'unknown' / 'not reported' entry",
    "lost_in_harmonisation": "the raw source has a value but the harmonised "
                             "value is missing",
    "missing_in_source_table": "missing in the harmonised source table; the raw "
                               "record was not supplied, so it is not traced further",
    "not_available_in_cohort": "the field does not exist for any patient of "
                               "this cohort",
    "no_expression_sample": "the patient has no expression sample, so no "
                            "expression-derived value exists",
    "scoring_failed": "the patient has an expression sample but the scorer "
                      "returned no value",
    "not_tested": "no assay data exist for this patient (no methylation array, "
                  "no sequenced sample of any sample type)",
    "tested_not_callable": "assay data exist but no admissible call could be made "
                           "(a required probe absent or NA, or the only assayed "
                           "sample is not a primary tumour)",
}

#: Call status: three states, never collapsed (decision D3).
CALL_STATUSES: tuple[str, ...] = ("called", "not_tested", "tested_not_callable")

#: Outcome columns. Labels; forbidden in any feature matrix.
LABEL_COLUMNS: tuple[str, ...] = ("overall_survival_months", "vital_status")


class FeatureStoreError(ValueError):
    """Base class: the store contract was violated."""


class DuplicateSourceError(FeatureStoreError):
    """A column was supplied by more than one input."""


class ShadowColumnError(FeatureStoreError):
    """A null-only column under the name of a computed feature was supplied by
    an input that does not own it."""


class ColumnOwnershipError(FeatureStoreError):
    """An input supplied a column it does not own."""


class StubInRealModeError(FeatureStoreError):
    """A stub or cohort-relative output was offered to a real-mode build."""


class LabelInFeatureMatrixError(FeatureStoreError):
    """An outcome label was requested as a feature."""


class UnexplainedMissingError(FeatureStoreError):
    """A null cell has no reason code."""


def sha256_file(path: Path | str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def refuse_external(path: Path | str) -> None:
    """The store is built from the development cohort only (no CGGA)."""
    if "cgga" in str(path).lower():
        raise FeatureStoreError(
            f"refusing {path}: the held-out cohort is read only through "
            "gliostrata.evaluation.read_guard, never by the feature-store build"
        )


@dataclass(frozen=True)
class ColumnSpec:
    name: str
    role: str
    source_layer: str
    source_input: str
    klass: str
    description: str
    scoring: str | None = None
    source_files: tuple[dict, ...] = ()

    def record(self, frame: pd.DataFrame) -> dict:
        col = frame[self.name]
        return {
            "name": self.name,
            "role": self.role,
            "source_layer": self.source_layer,
            "source_input": self.source_input,
            "source_files": [dict(f) for f in self.source_files],
            "class": self.klass,
            "dtype": str(col.dtype),
            "n_missing": int(col.isna().sum()),
            "scoring": self.scoring,
            "description": self.description,
        }


@dataclass
class FeatureStore:
    """The assembled store. Build it with ``build_feature_store``."""

    mode: str
    features: pd.DataFrame
    labels: pd.DataFrame
    missing: pd.DataFrame
    columns: dict[str, ColumnSpec]
    meta: dict = field(default_factory=dict)

    @property
    def patient_ids(self) -> list[str]:
        return self.features["patient_id"].astype(str).tolist()

    def column_records(self) -> list[dict]:
        records = []
        for name, spec in self.columns.items():
            frame = self.labels if spec.role == "label" else self.features
            records.append(spec.record(frame))
        return records

    def reason_counts(self) -> dict[str, dict[str, int]]:
        out: dict[str, dict[str, int]] = {}
        if self.missing.empty:
            return out
        grouped = self.missing.groupby(["column", "reason"], sort=True).size()
        for (col, reason), n in grouped.items():
            out.setdefault(str(col), {})[str(reason)] = int(n)
        return out

    def metadata_document(self) -> dict:
        return {
            "store_version": STORE_VERSION,
            "mode": self.mode,
            "n_patients": len(self.features),
            "vocabularies": {
                "roles": list(ROLES),
                "classes": COLUMN_CLASSES,
                "missing_reasons": MISSING_REASONS,
                "call_statuses": list(CALL_STATUSES),
            },
            "meta": self.meta,
            "columns": self.column_records(),
        }


def sibling_paths(store_path: Path | str) -> dict[str, Path]:
    """Default paths of the four store files, derived from the features path."""
    path = Path(store_path)
    stem = path.with_suffix("")
    return {
        "features": path,
        "labels": Path(f"{stem}_labels.parquet"),
        "missing": Path(f"{stem}_missing.parquet"),
        "columns": Path(f"{stem}_columns.json"),
    }


def _write_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = frame.copy()
    frame.attrs = {}
    frame.to_parquet(path, index=False, engine="pyarrow")


def write_store(store: FeatureStore, paths: dict[str, Path | str]) -> dict[str, dict]:
    """Write the four files; return ``{name: {path, bytes, sha256}}``."""
    paths = {k: Path(v) for k, v in paths.items()}
    _write_parquet(store.features, paths["features"])
    _write_parquet(store.labels, paths["labels"])
    _write_parquet(store.missing, paths["missing"])
    paths["columns"].parent.mkdir(parents=True, exist_ok=True)
    paths["columns"].write_text(
        json.dumps(store.metadata_document(), indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return {
        name: {"path": str(p).replace("\\", "/"), "bytes": p.stat().st_size,
               "sha256": sha256_file(p)}
        for name, p in paths.items()
    }


def read_store_metadata(columns_path: Path | str) -> dict:
    return json.loads(Path(columns_path).read_text(encoding="utf-8"))
