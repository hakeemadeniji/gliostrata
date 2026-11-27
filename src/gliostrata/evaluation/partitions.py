"""Deterministic, freezable cohort partitioning (defensibility point 2).

Why hashing and not shuffling
-----------------------------
A shuffle-based split (``train_test_split``, ``df.sample(frac=...)``) is a
function of the *whole table*: the random stream is consumed row by row, so the
assignment of patient *P* depends on how many rows precede it and how many rows
exist at all. Change an upstream QC filter, add a cohort, sort differently, or
upgrade pandas/NumPy, and every patient can move between train and test. For a
partition that is supposed to stay frozen for the life of the project that is
fatal: once a model has trained on a patient who was previously in the test
set, the external test is contaminated and cannot be repaired.

Here the split is a pure function of ``(seed, patient_id)``:

    u = sha256(f"{seed}|{patient_id}") mapped into [0, 1)
    test  if u < test_fraction  else  train

Consequences, all of them properties the test suite asserts:

* **Row order independent** -- shuffling the input cannot change the split.
* **Cohort size independent** -- adding or removing *other* patients cannot
  move a patient already assigned.
* **Library-version independent** -- SHA-256 and integer arithmetic only; no
  NumPy/pandas RNG is involved.
* **Reconstructable** -- the split can be recomputed from the seed alone,
  without the original table.

Stratification
--------------
Stratification salts the hash with the stratum label
(``sha256(f"{seed}|{stratum}|{patient_id}")``) instead of ranking patients
within each stratum. Ranking would give exact per-stratum quotas but would
re-introduce set dependence (a patient's rank changes when other patients are
added), which is precisely the property this module exists to avoid. Salted
hashing keeps strict per-patient determinism and makes each stratum's expected
test fraction equal to ``test_fraction`` independently of the others; the
achieved fraction is close to, but not exactly, the target. That trade-off is
deliberate and is recorded in ``docs/preregistration.md``.

All-test partitions
-------------------
A test fraction of exactly ``1.0`` holds out every patient: the training list is
empty and the test list is every unique patient id. That is the correct seal for
an external cohort, which has no training data. It is produced by the same hash
rule (every ``u`` lies in ``[0, 1)`` and is therefore below 1.0), so it keeps the
determinism properties above and needs no second code path. The CLI spells it
``--all-test``.

Sealing
-------
:func:`write_manifest` with ``seal=True`` marks a manifest as frozen. Re-running
against a sealed manifest refuses to overwrite it unless ``force=True`` is
passed explicitly, and the refusal explains the scientific consequence rather
than just reporting that a file exists.
"""
from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

__all__ = [
    "MANIFEST_SCHEMA_VERSION",
    "PartitionSealedError",
    "PartitionResult",
    "PartitionManifest",
    "DriftReport",
    "hash_unit_interval",
    "assign_split",
    "partition_cohort",
    "stratum_balance",
    "sha256_over_ids",
    "file_content_sha256",
    "build_manifest",
    "freeze_partition",
    "write_manifest",
    "read_manifest",
    "verify_partition",
]

MANIFEST_SCHEMA_VERSION = "1.0"
_SPLITS = ("train", "test")


class PartitionSealedError(RuntimeError):
    """Raised when a write would overwrite a sealed (frozen) partition manifest."""


# ---------------------------------------------------------------------------
# Deterministic hashing primitives
# ---------------------------------------------------------------------------


def hash_unit_interval(key: str) -> float:
    """Map an arbitrary string deterministically and uniformly into ``[0, 1)``.

    Uses the top 64 bits of the SHA-256 digest. Stable across processes,
    platforms, Python versions and library versions (unlike :func:`hash`,
    which is salted per process).
    """
    digest = hashlib.sha256(key.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") / 2**64


def _split_key(patient_id: str, seed: int, stratum: str | None) -> str:
    if stratum is None:
        return f"{seed}|{patient_id}"
    return f"{seed}|{stratum}|{patient_id}"


def assign_split(
    patient_id: str,
    *,
    seed: int,
    test_fraction: float = 0.30,
    stratum: str | None = None,
) -> str:
    """Return ``"train"`` or ``"test"`` for one patient.

    Depends only on its arguments -- never on any other patient, on row order,
    or on cohort size. A ``test_fraction`` of exactly ``1.0`` assigns ``"test"``
    always: an all-test partition, for an external cohort with no training data.
    """
    if not 0.0 < test_fraction <= 1.0:
        raise ValueError(f"test_fraction must be in (0, 1], got {test_fraction!r}.")
    u = hash_unit_interval(_split_key(str(patient_id), seed, stratum))
    return "test" if u < test_fraction else "train"


def sha256_over_ids(ids: Iterable[str]) -> str:
    """SHA-256 over the newline-joined, sorted, de-duplicated id list."""
    payload = "\n".join(sorted({str(i) for i in ids}))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def file_content_sha256(path: str | Path, *, chunk_size: int = 1 << 20) -> str:
    """SHA-256 over a file's bytes (streamed, so large Parquet files are fine)."""
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Partitioning
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class PartitionResult:
    """A computed partition: sorted id lists plus the parameters that made it."""

    seed: int
    test_fraction: float
    patient_id_col: str
    stratify_by: str | None
    train_ids: tuple[str, ...]
    test_ids: tuple[str, ...]
    strata: Mapping[str, str] = field(default_factory=dict)

    @property
    def n_train(self) -> int:
        return len(self.train_ids)

    @property
    def n_test(self) -> int:
        return len(self.test_ids)

    @property
    def n_total(self) -> int:
        return self.n_train + self.n_test

    def assignments(self) -> dict[str, str]:
        """Return ``{patient_id: "train"|"test"}``."""
        out = {pid: "train" for pid in self.train_ids}
        out.update({pid: "test" for pid in self.test_ids})
        return out

    def ids(self, split: str) -> tuple[str, ...]:
        if split not in _SPLITS:
            raise ValueError(f"split must be one of {_SPLITS}, got {split!r}.")
        return self.train_ids if split == "train" else self.test_ids

    def label_frame(self) -> pd.DataFrame:
        """Return a two-column DataFrame of ``(patient_id, split)``, sorted by id."""
        items = sorted(self.assignments().items())
        return pd.DataFrame(items, columns=[self.patient_id_col, "split"])


def partition_cohort(
    cohort: pd.DataFrame,
    *,
    seed: int,
    test_fraction: float = 0.30,
    patient_id_col: str = "patient_id",
    stratify_by: str | None = None,
    na_stratum: str = "__NA__",
) -> PartitionResult:
    """Partition a cohort into train/test deterministically by hashed patient id.

    Parameters
    ----------
    cohort:
        Table with one or more rows per patient. Duplicate patient ids are
        collapsed -- a patient is assigned once, so multi-sample patients can
        never straddle the split (which would leak the test set).
    seed:
        Integer seed baked into the hash. Part of the frozen manifest; changing
        it produces a different partition and is a pre-registration amendment.
    test_fraction:
        Target fraction of patients held out. Achieved fraction is close to but
        not exactly this value, by construction (see module docstring). The one
        exception is exactly ``1.0``, which holds out every patient: an
        all-test partition with an empty training list, for an external cohort
        that has no training data.
    stratify_by:
        Optional column whose value stratifies the hash. Must be constant per
        patient; a patient with conflicting values raises.
    na_stratum:
        Label used for missing stratification values, so that NaN does not
        silently create per-row strata.

    Raises
    ------
    ValueError
        For a missing id/stratum column, an empty cohort, null patient ids, a
        ``test_fraction`` outside ``(0, 1]``, or a patient with inconsistent
        stratum values.
    """
    if patient_id_col not in cohort.columns:
        raise ValueError(
            f"Cohort has no '{patient_id_col}' column; columns are {list(cohort.columns)}."
        )
    if not 0.0 < test_fraction <= 1.0:
        raise ValueError(f"test_fraction must be in (0, 1], got {test_fraction!r}.")
    if len(cohort) == 0:
        raise ValueError("Cannot partition an empty cohort.")

    ids = cohort[patient_id_col]
    if ids.isna().any():
        raise ValueError(
            f"{int(ids.isna().sum())} row(s) have a null {patient_id_col}; "
            "a frozen partition cannot be defined over unidentified patients."
        )

    strata: dict[str, str] = {}
    if stratify_by is not None:
        if stratify_by not in cohort.columns:
            raise ValueError(
                f"stratify_by column '{stratify_by}' not in cohort; "
                f"columns are {list(cohort.columns)}."
            )
        for pid, group in cohort.groupby(ids.astype(str), sort=True)[stratify_by]:
            labels = {na_stratum if pd.isna(v) else str(v) for v in group}
            if len(labels) > 1:
                raise ValueError(
                    f"Patient '{pid}' has conflicting '{stratify_by}' values "
                    f"{sorted(labels)}; stratification must be constant per patient."
                )
            strata[str(pid)] = labels.pop()

    unique_ids = sorted({str(i) for i in ids})
    train: list[str] = []
    test: list[str] = []
    for pid in unique_ids:
        stratum = strata.get(pid) if stratify_by is not None else None
        target = test if assign_split(
            pid, seed=seed, test_fraction=test_fraction, stratum=stratum
        ) == "test" else train
        target.append(pid)

    return PartitionResult(
        seed=int(seed),
        test_fraction=float(test_fraction),
        patient_id_col=patient_id_col,
        stratify_by=stratify_by,
        train_ids=tuple(train),
        test_ids=tuple(test),
        strata=dict(strata),
    )


def stratum_balance(result: PartitionResult) -> dict[str, dict[str, float | int]]:
    """Per-stratum achieved test fraction, for reporting and balance checks.

    Returns ``{stratum: {"n": int, "n_test": int, "test_fraction": float}}``.
    Empty when the partition was not stratified.
    """
    if result.stratify_by is None:
        return {}
    assignments = result.assignments()
    out: dict[str, dict[str, float | int]] = {}
    for pid, stratum in result.strata.items():
        row = out.setdefault(stratum, {"n": 0, "n_test": 0, "test_fraction": 0.0})
        row["n"] = int(row["n"]) + 1
        if assignments.get(pid) == "test":
            row["n_test"] = int(row["n_test"]) + 1
    for row in out.values():
        row["test_fraction"] = int(row["n_test"]) / int(row["n"]) if row["n"] else 0.0
    return out


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class PartitionManifest:
    """The on-disk record that freezes a partition.

    Everything needed to (a) recompute the partition from scratch and (b) prove
    that the cohort it was computed from has not changed.
    """

    schema_version: str
    partition_name: str
    seed: int
    test_fraction: float
    patient_id_col: str
    stratify_by: str | None
    cohort_source: str
    cohort_sha256: str
    cohort_rows: int
    created_utc: str
    sealed: bool
    n_train: int
    n_test: int
    train_ids: tuple[str, ...]
    test_ids: tuple[str, ...]
    train_ids_sha256: str
    test_ids_sha256: str
    all_ids_sha256: str
    stratum_balance: Mapping[str, Mapping[str, float | int]] = field(default_factory=dict)
    notes: str = ""

    def to_dict(self) -> dict[str, object]:
        row = asdict(self)
        row["train_ids"] = list(self.train_ids)
        row["test_ids"] = list(self.test_ids)
        return row

    @classmethod
    def from_dict(cls, data: Mapping[str, object]) -> PartitionManifest:
        known = {f for f in cls.__dataclass_fields__}  # type: ignore[attr-defined]
        missing = known - set(data) - {"stratum_balance", "notes"}
        if missing:
            raise ValueError(f"Manifest is missing required field(s): {sorted(missing)}.")
        kwargs = {k: v for k, v in data.items() if k in known}
        kwargs["train_ids"] = tuple(str(i) for i in kwargs.get("train_ids", ()))  # type: ignore[arg-type]
        kwargs["test_ids"] = tuple(str(i) for i in kwargs.get("test_ids", ()))  # type: ignore[arg-type]
        return cls(**kwargs)  # type: ignore[arg-type]

    def assignments(self) -> dict[str, str]:
        out = {pid: "train" for pid in self.train_ids}
        out.update({pid: "test" for pid in self.test_ids})
        return out

    def as_partition_params(self) -> dict[str, object]:
        """The arguments needed to recompute this partition."""
        return {
            "seed": self.seed,
            "test_fraction": self.test_fraction,
            "patient_id_col": self.patient_id_col,
            "stratify_by": self.stratify_by,
        }


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_manifest(
    result: PartitionResult,
    *,
    cohort_source: str | Path,
    partition_name: str,
    cohort_rows: int,
    sealed: bool = True,
    notes: str = "",
    cohort_sha256: str | None = None,
) -> PartitionManifest:
    """Assemble a :class:`PartitionManifest` (does not write anything)."""
    source = Path(cohort_source)
    return PartitionManifest(
        schema_version=MANIFEST_SCHEMA_VERSION,
        partition_name=partition_name,
        seed=result.seed,
        test_fraction=result.test_fraction,
        patient_id_col=result.patient_id_col,
        stratify_by=result.stratify_by,
        cohort_source=source.as_posix(),
        cohort_sha256=cohort_sha256 if cohort_sha256 is not None else file_content_sha256(source),
        cohort_rows=int(cohort_rows),
        created_utc=_utc_now_iso(),
        sealed=bool(sealed),
        n_train=result.n_train,
        n_test=result.n_test,
        train_ids=result.train_ids,
        test_ids=result.test_ids,
        train_ids_sha256=sha256_over_ids(result.train_ids),
        test_ids_sha256=sha256_over_ids(result.test_ids),
        all_ids_sha256=sha256_over_ids((*result.train_ids, *result.test_ids)),
        stratum_balance=stratum_balance(result),
        notes=notes,
    )


_SEAL_REFUSAL = (
    "Refusing to overwrite the sealed partition manifest at {path}.\n"
    "\n"
    "This is not a file-permission problem -- it is the point of the seal. "
    "'{name}' was frozen at {created} and the models trained since then have "
    "been fitted on its training split. Rewriting it now would silently move "
    "patients between train and test. Any patient moved from train to test "
    "would be scored by a model that has already seen them, which makes the "
    "'external' test result meaningless and is not recoverable by retraining "
    "-- the information has already leaked into every hyper-parameter and "
    "model-selection decision made on this partition.\n"
    "\n"
    "If the partition genuinely must change (for example the source cohort was "
    "re-released), that is a pre-registration amendment, not a re-run: record "
    "the reason and date in docs/preregistration.md, archive the old manifest "
    "alongside the new one, re-run every model from scratch, and only then "
    "pass force=True (CLI: --force)."
)


def write_manifest(
    manifest: PartitionManifest,
    path: str | Path,
    *,
    force: bool = False,
) -> Path:
    """Write a manifest to JSON, refusing to clobber a sealed one.

    Raises
    ------
    PartitionSealedError
        If ``path`` already holds a manifest with ``sealed: true`` and
        ``force`` is not ``True``. The message explains the scientific
        consequence of overwriting, not merely that the file exists.
    """
    p = Path(path)
    if p.exists() and not force:
        try:
            existing = read_manifest(p)
        except (ValueError, json.JSONDecodeError):
            existing = None
        if existing is not None and existing.sealed:
            raise PartitionSealedError(
                _SEAL_REFUSAL.format(
                    path=p, name=existing.partition_name, created=existing.created_utc
                )
            )
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(
        json.dumps(manifest.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return p


def read_manifest(path: str | Path) -> PartitionManifest:
    """Read a manifest back from JSON."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Partition manifest not found: {p}")
    return PartitionManifest.from_dict(json.loads(p.read_text(encoding="utf-8")))


def freeze_partition(
    cohort: pd.DataFrame,
    *,
    cohort_source: str | Path,
    manifest_path: str | Path,
    partition_name: str,
    seed: int,
    test_fraction: float = 0.30,
    patient_id_col: str = "patient_id",
    stratify_by: str | None = None,
    seal: bool = True,
    force: bool = False,
    notes: str = "",
) -> tuple[PartitionResult, PartitionManifest, Path]:
    """Partition a cohort, build its manifest and write it, honouring the seal.

    Convenience wrapper over :func:`partition_cohort`, :func:`build_manifest`
    and :func:`write_manifest`. The seal check happens *before* anything is
    written.
    """
    result = partition_cohort(
        cohort,
        seed=seed,
        test_fraction=test_fraction,
        patient_id_col=patient_id_col,
        stratify_by=stratify_by,
    )
    manifest = build_manifest(
        result,
        cohort_source=cohort_source,
        partition_name=partition_name,
        cohort_rows=len(cohort),
        sealed=seal,
        notes=notes,
    )
    written = write_manifest(manifest, manifest_path, force=force)
    return result, manifest, written


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class DriftReport:
    """Result of checking a live cohort against a frozen manifest."""

    ok: bool
    partition_name: str
    manifest_created_utc: str
    cohort_source: str
    cohort_hash_matches: bool
    expected_cohort_sha256: str
    observed_cohort_sha256: str | None
    expected_rows: int
    observed_rows: int
    added_ids: tuple[str, ...]
    removed_ids: tuple[str, ...]
    moved_ids: tuple[tuple[str, str, str], ...]
    expected_n_train: int
    expected_n_test: int
    observed_n_train: int
    observed_n_test: int
    messages: tuple[str, ...]

    @property
    def n_drifted(self) -> int:
        return len(self.added_ids) + len(self.removed_ids) + len(self.moved_ids)

    def describe(self) -> str:
        """Multi-line human-readable report."""
        head = (
            f"Partition '{self.partition_name}' (frozen {self.manifest_created_utc}): "
            f"{'VERIFIED' if self.ok else 'DRIFTED'}"
        )
        return "\n".join((head, *(f"  - {m}" for m in self.messages)))

    def to_dict(self) -> dict[str, object]:
        row = asdict(self)
        row["added_ids"] = list(self.added_ids)
        row["removed_ids"] = list(self.removed_ids)
        row["moved_ids"] = [list(m) for m in self.moved_ids]
        row["messages"] = list(self.messages)
        return row


def verify_partition(
    manifest: PartitionManifest | str | Path,
    cohort: pd.DataFrame,
    *,
    cohort_source: str | Path | None = None,
    max_listed: int = 20,
) -> DriftReport:
    """Check that ``cohort`` still reproduces the frozen partition.

    Recomputes the partition from the manifest's own parameters and compares
    it, patient by patient, against the frozen id lists. Reports precisely what
    drifted: ids added, ids removed, and ids that moved between splits (the
    dangerous case -- a train->test move means the model has already seen a
    test patient).

    ``cohort_source`` is optional; when given, the file's content hash is also
    compared against the manifest's, which catches a cohort that was rewritten
    in a way the id list alone cannot see (e.g. an edited survival column).
    """
    if not isinstance(manifest, PartitionManifest):
        manifest = read_manifest(manifest)

    messages: list[str] = []
    result = partition_cohort(
        cohort,
        seed=manifest.seed,
        test_fraction=manifest.test_fraction,
        patient_id_col=manifest.patient_id_col,
        stratify_by=manifest.stratify_by,
    )
    frozen = manifest.assignments()
    current = result.assignments()

    added = tuple(sorted(set(current) - set(frozen)))
    removed = tuple(sorted(set(frozen) - set(current)))
    moved = tuple(
        (pid, frozen[pid], current[pid])
        for pid in sorted(set(frozen) & set(current))
        if frozen[pid] != current[pid]
    )

    observed_hash: str | None = None
    hash_matches = True
    if cohort_source is not None:
        src = Path(cohort_source)
        if src.exists():
            observed_hash = file_content_sha256(src)
            hash_matches = observed_hash == manifest.cohort_sha256
            if not hash_matches:
                messages.append(
                    f"Source file content hash changed: manifest recorded "
                    f"{manifest.cohort_sha256[:16]}..., file is now "
                    f"{observed_hash[:16]}.... The cohort file was rewritten; the id "
                    f"list may still match while column values have changed."
                )
        else:
            hash_matches = False
            messages.append(
                f"Recorded cohort source '{src}' does not exist, so the content "
                f"hash could not be checked."
            )

    def _preview(items: Sequence[object]) -> str:
        shown = [str(i) for i in items[:max_listed]]
        suffix = f" (+{len(items) - max_listed} more)" if len(items) > max_listed else ""
        return ", ".join(shown) + suffix

    if added:
        messages.append(
            f"{len(added)} patient id(s) present in the cohort but not in the frozen "
            f"partition: {_preview(added)}"
        )
    if removed:
        messages.append(
            f"{len(removed)} patient id(s) in the frozen partition but missing from the "
            f"cohort: {_preview(removed)}"
        )
    if moved:
        to_test = [m for m in moved if m[2] == "test"]
        messages.append(
            f"{len(moved)} patient id(s) changed split: "
            f"{_preview([f'{p}: {a}->{b}' for p, a, b in moved])}"
        )
        if to_test:
            messages.append(
                f"CONTAMINATION RISK: {len(to_test)} patient(s) moved train->test and "
                f"have therefore already been seen by models fitted on this partition."
            )
    if manifest.cohort_rows != len(cohort):
        messages.append(
            f"Row count changed: manifest recorded {manifest.cohort_rows}, cohort now "
            f"has {len(cohort)}. (Harmless on its own -- the hash split is row-count "
            f"independent -- but it means upstream filtering changed.)"
        )

    ok = not added and not removed and not moved and hash_matches
    if ok and not messages:
        messages.append(
            f"All {result.n_total} patient ids reproduce their frozen split "
            f"({result.n_train} train / {result.n_test} test); id-list and source "
            f"content hashes match."
        )
    elif ok:
        messages.append(
            f"All {result.n_total} patient ids reproduce their frozen split "
            f"({result.n_train} train / {result.n_test} test)."
        )

    return DriftReport(
        ok=ok,
        partition_name=manifest.partition_name,
        manifest_created_utc=manifest.created_utc,
        cohort_source=manifest.cohort_source,
        cohort_hash_matches=hash_matches,
        expected_cohort_sha256=manifest.cohort_sha256,
        observed_cohort_sha256=observed_hash,
        expected_rows=manifest.cohort_rows,
        observed_rows=len(cohort),
        added_ids=added,
        removed_ids=removed,
        moved_ids=moved,
        expected_n_train=manifest.n_train,
        expected_n_test=manifest.n_test,
        observed_n_train=result.n_train,
        observed_n_test=result.n_test,
        messages=tuple(messages),
    )
