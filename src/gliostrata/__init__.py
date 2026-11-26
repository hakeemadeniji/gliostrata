"""GlioStrata — patient-specific computational framework for stem-cell-delivered GBM therapy.

Eight-layer architecture (L0–L7):
    L0 — Public data sources (TCGA, CGGA, CPTAC, scRNA atlases, GLASS, TCIA)
    L1 — Ingestion & harmonization
    L2 — Tumor profiling & subtyping (Verhaak + Neftel + TME + radiomics)
    L3 — Unified patient feature store
    L4 — Mechanistic simulation (PhysiCell tumor growth + NSC migration + payload kinetics)
    L5 — Outcome prediction (DeepSurv + response + pseudoprogression)
    L6 — Therapeutic configuration optimizer (BoTorch Bayesian optimization)
    L7 — Wet-lab handoff & recalibration
"""

__version__ = "0.1.0"

from gliostrata.config import load_config

__all__ = ["__version__", "load_config"]
