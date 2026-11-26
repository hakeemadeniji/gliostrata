"""Pre-registration and frozen-evaluation machinery (defensibility points 1-3).

The GlioStrata protocol (Chapter 8) wraps every test in a ten-point
defensibility stack. The first three points are infrastructure, and they must
exist *before* any model is trained:

1. **Pre-registration** -- the analysis plan and thresholds are deposited on
   OSF before test data is touched. See ``docs/preregistration.md`` for the
   human-readable draft and :mod:`gliostrata.evaluation.thresholds` for the
   machine-readable, hash-protected threshold table.
2. **Frozen test sets** -- external cohort partitions are sealed and never
   recomputed. See :mod:`gliostrata.evaluation.partitions`.
3. **Pre-specified thresholds** -- every threshold is committed before
   evaluation, and evaluating an unregistered metric raises.

Once a model has seen a test partition, that partition's value as an external
test is permanently gone; nothing downstream can repair it. These modules are
deliberately strict for that reason.

Typical use::

    from gliostrata.evaluation import freeze_partition, verify_partition, evaluate

    result, manifest, path = freeze_partition(
        cohort,
        cohort_source="data/raw/tcga_gbm/clinical_cohort.parquet",
        manifest_path="data/partitions/tcga_gbm_clinical.json",
        partition_name="tcga_gbm_clinical_v1",
        seed=20251126,
    )
    report = verify_partition(manifest, cohort)   # later, before every eval run
    verdict = evaluate("T2_survival_cindex_cgga", 0.68, 0.661, 0.702)
"""
from __future__ import annotations

from gliostrata.evaluation.partitions import (
    MANIFEST_SCHEMA_VERSION,
    DriftReport,
    PartitionManifest,
    PartitionResult,
    PartitionSealedError,
    assign_split,
    build_manifest,
    file_content_sha256,
    freeze_partition,
    hash_unit_interval,
    partition_cohort,
    read_manifest,
    sha256_over_ids,
    stratum_balance,
    verify_partition,
    write_manifest,
)
from gliostrata.evaluation.thresholds import (
    THRESHOLD_TABLE_SHA256,
    THRESHOLDS,
    CIBound,
    Direction,
    Status,
    Threshold,
    UnregisteredMetricError,
    Verdict,
    canonical_table,
    evaluate,
    evaluate_all,
    format_threshold_table,
    get_threshold,
    registered_metrics,
    threshold_table_checksum,
)

__all__ = [
    # thresholds
    "THRESHOLDS",
    "THRESHOLD_TABLE_SHA256",
    "CIBound",
    "Direction",
    "Status",
    "Threshold",
    "UnregisteredMetricError",
    "Verdict",
    "canonical_table",
    "evaluate",
    "evaluate_all",
    "format_threshold_table",
    "get_threshold",
    "registered_metrics",
    "threshold_table_checksum",
    # partitions
    "MANIFEST_SCHEMA_VERSION",
    "DriftReport",
    "PartitionManifest",
    "PartitionResult",
    "PartitionSealedError",
    "assign_split",
    "build_manifest",
    "file_content_sha256",
    "freeze_partition",
    "hash_unit_interval",
    "partition_cohort",
    "read_manifest",
    "sha256_over_ids",
    "stratum_balance",
    "verify_partition",
    "write_manifest",
]
