"""L1 — Quality control gates.

Stub: real implementation checks sequencing depth, mapping rate, contamination,
tumor purity, and molecular consistency (e.g. IDH from WES == IDH from methylation).

The `qc_report` function below returns a per-sample QC summary that tests can use.
"""
from __future__ import annotations

import pandas as pd


def qc_report(cohort: pd.DataFrame) -> pd.DataFrame:
    """Generate a per-sample QC report.

    Returns a DataFrame with one row per sample and QC columns:
      - qc_pass: bool, whether the sample passes all gates
      - qc_flags: list of failure reasons (empty if pass)
      - estimated_purity: float, fraction (stub: returns 0.6 + noise)
    """
    n = len(cohort)
    report = pd.DataFrame({
        "patient_id": cohort["patient_id"].values,
        "qc_pass": [True] * n,
        "qc_flags": [[] for _ in range(n)],
        "estimated_purity": [0.60] * n,  # stub — real impl uses ABSOLUTE/ESTIMATE
    })

    # Gate 1: survival time must be positive
    bad_survival = cohort["overall_survival_months"] <= 0
    report.loc[bad_survival.values, "qc_pass"] = False
    for i in bad_survival[bad_survival].index:
        report.at[i, "qc_flags"].append("non_positive_survival")

    return report
