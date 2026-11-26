"""L0 — Public data ingestion.

Currently a stub: real TCGA/CGGA/CPTAC ingestion requires TCGAbiolinks (R) called
via subprocess, or direct HTTP download. The `generate_synthetic_cohort` function
below is fully working and powers all CLI dry-runs and tests.

See docs/roadmap.md for implementation order.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# Synthetic cohort schema — kept consistent across L0–L3 so downstream modules
# can rely on stable column names.
COHORT_COLUMNS = [
    "patient_id",
    "age",
    "sex",
    "kps",  # Karnofsky Performance Score
    "extent_of_resection",  # biopsy / STR / GTR
    "mgmt_methylated",
    "idh_mutated",
    "subtype",  # Verhaak: Classical / Mesenchymal / Proneural / Neural
    # Neftel single-cell state fractions (sum ~ 1 over malignant cells)
    "state_npc", "state_opc", "state_ac", "state_mes",
    # TME fractions (sum ~ 1 over all cells incl. non-malignant)
    "tme_immune", "tme_stromal", "tme_endothelial",
    # Longitudinal outcome
    "overall_survival_months",
    "vital_status",  # 1 = death, 0 = censored
    "treatment_arm",  # Stupp / experimental
]


def generate_synthetic_cohort(n: int = 50, seed: int = 42) -> pd.DataFrame:
    """Generate a synthetic GBM cohort with realistic feature distributions.

    Distributions are calibrated to roughly match published TCGA-GBM summary
    statistics (age, sex, MGMT frequency, subtype proportions, survival).
    Numbers are illustrative — not for inference.
    """
    rng = np.random.default_rng(seed)

    # Subtype proportions (Verhaak, approx)
    subtype_probs = {"Classical": 0.27, "Mesenchymal": 0.32, "Proneural": 0.30, "Neural": 0.11}
    subtypes = rng.choice(list(subtype_probs.keys()), size=n, p=list(subtype_probs.values()))

    # MGMT methylation (~ 45% in unselected GBM)
    mgmt = rng.binomial(1, 0.45, size=n).astype(bool)
    # IDH mutation (~ 8% in unselected GBM, ~ 80% in lower-grade)
    idh = rng.binomial(1, 0.08, size=n).astype(bool)

    # Demographics
    age = rng.normal(60, 12, size=n).clip(20, 90).astype(int)
    sex = rng.binomial(1, 0.56, size=n)  # ~ 56% male
    kps = rng.normal(80, 12, size=n).clip(40, 100).astype(int)
    eor = rng.choice(["biopsy", "STR", "GTR"], size=n, p=[0.20, 0.35, 0.45])

    # Neftel state fractions — depend on subtype (MES subtype has more MES-state cells)
    states = np.zeros((n, 4))
    for i, st in enumerate(subtypes):
        if st == "Mesenchymal":
            base = np.array([0.05, 0.10, 0.20, 0.65])
        elif st == "Classical":
            base = np.array([0.20, 0.15, 0.55, 0.10])
        elif st == "Proneural":
            base = np.array([0.30, 0.40, 0.20, 0.10])
        else:  # Neural
            base = np.array([0.25, 0.30, 0.30, 0.15])
        noisy = base + rng.normal(0, 0.05, size=4)
        noisy = np.clip(noisy, 0.01, None)
        states[i] = noisy / noisy.sum()
    state_npc, state_opc, state_ac, state_mes = states.T

    # TME fractions
    tme_immune = rng.uniform(0.05, 0.30, size=n)
    tme_stromal = rng.uniform(0.10, 0.25, size=n)
    tme_endothelial = rng.uniform(0.02, 0.10, size=n)

    # Overall survival — depends on subtype, MGMT, IDH, age, KPS
    # Base median 14 months, modifiers per published literature
    base_log_os = np.log(14)
    log_os = (
        base_log_os
        + 0.25 * mgmt.astype(float)         # MGMT methylated → +~30% survival
        + 1.0 * idh.astype(float)           # IDH mutated → ~3× survival
        + 0.30 * (subtypes == "Proneural").astype(float)
        - 0.20 * (subtypes == "Mesenchymal").astype(float)
        - 0.015 * (age - 60)
        + 0.015 * (kps - 80)
        - 0.20 * (eor == "biopsy").astype(float)
        + 0.10 * (eor == "GTR").astype(float)
        + rng.normal(0, 0.35, size=n)        # individual variation
    )
    os_months = np.exp(log_os)
    # Censoring ~ 30%
    vital = rng.binomial(1, 0.70, size=n)
    # Censored patients have longer observed survival
    observed_os = np.where(vital.astype(bool), os_months, os_months * rng.uniform(1.1, 1.5, size=n))

    treatment = rng.choice(["Stupp", "experimental"], size=n, p=[0.8, 0.2])

    df = pd.DataFrame({
        "patient_id": [f"GS-SYN-{i+1:04d}" for i in range(n)],
        "age": age,
        "sex": ["M" if s == 1 else "F" for s in sex],
        "kps": kps,
        "extent_of_resection": eor,
        "mgmt_methylated": mgmt,
        "idh_mutated": idh,
        "subtype": subtypes,
        "state_npc": state_npc,
        "state_opc": state_opc,
        "state_ac": state_ac,
        "state_mes": state_mes,
        "tme_immune": tme_immune,
        "tme_stromal": tme_stromal,
        "tme_endothelial": tme_endothelial,
        "overall_survival_months": observed_os,
        "vital_status": vital,
        "treatment_arm": treatment,
    })
    return df[COHORT_COLUMNS]


def download_tcga(out_dir: str, *, n_samples: int = 600) -> None:
    """STUB — download TCGA-GBM data.

    Real implementation:
      1. Use TCGAbiolinks (R) via subprocess.run(["Rscript", ...]) to query GDC.
      2. Save raw RNA-seq (FPKM), methylation, and clinical data to out_dir.
      3. Requires R + TCGAbiolinks installed. Controlled-access data needs
         dbGaP API key in .env.

    See docs/roadmap.md for the implementation plan.
    """
    raise NotImplementedError(
        "TCGA ingestion is a stub. Use generate_synthetic_cohort() for development. "
        "See docs/roadmap.md."
    )


def download_cgga(out_dir: str) -> None:
    """STUB — download CGGA data."""
    raise NotImplementedError("CGGA ingestion is a stub. See docs/roadmap.md.")


def download_cptac(out_dir: str) -> None:
    """STUB — download CPTAC-GBM proteomics data."""
    raise NotImplementedError("CPTAC ingestion is a stub. See docs/roadmap.md.")
