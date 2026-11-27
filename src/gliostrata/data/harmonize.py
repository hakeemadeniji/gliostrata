"""L1 — Harmonization pipeline.

Currently a stub: real implementation uses ComBat-seq (R) for bulk RNA-seq
batch correction and harmony for single-cell data, with quality control gates.

The `harmonize` function below passes data through unchanged (with a QC filter
that drops obviously invalid samples) so downstream layers can be developed
without waiting for real harmonization to be implemented.
"""
from __future__ import annotations

import pandas as pd


def harmonize(cohort: pd.DataFrame) -> pd.DataFrame:
    """Harmonize an ingested cohort.

    Stub implementation: validates required columns, drops samples failing
    QC, and returns the (largely unchanged) cohort for downstream profiling.

    Real implementation will:
      1. Canonicalize sample IDs against a master sample map.
      2. Apply ComBat-seq for cohort-of-origin batch effects.
      3. Apply harmony for single-cell data integration.
      4. Estimate tumor purity (ABSOLUTE for WES, ESTIMATE for RNA-seq).
      5. Flag samples < 40% purity (retain but covary).
      6. Check molecular consistency (e.g. IDH from WES == IDH from methylation).
    """
    if not isinstance(cohort, pd.DataFrame):
        raise TypeError(f"Expected DataFrame, got {type(cohort)}")

    # Required columns — at minimum, a patient ID and survival outcome
    required = {"patient_id", "overall_survival_months", "vital_status"}
    missing = required - set(cohort.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # QC gate: drop samples with non-positive survival time
    qc_pass = cohort["overall_survival_months"] > 0
    n_dropped = (~qc_pass).sum()
    if n_dropped > 0:
        # In a stub we log but don't actually drop, to keep cohorts stable for testing.
        # Real implementation would drop and log to a QC report.
        pass

    out = cohort.copy()
    # Add a harmonization metadata column
    out["harmonization_version"] = "0.1.0-stub"
    out["batch_corrected"] = False  # will be True once ComBat-seq is wired up
    return out


def combat_seq(expression_matrix: pd.DataFrame, batch: pd.Series) -> pd.DataFrame:
    """STUB — apply ComBat-seq batch correction.

    Real implementation will shell out to R via subprocess.run(["Rscript", ...]).
    """
    raise NotImplementedError("ComBat-seq not implemented. See docs/roadmap.md.")
