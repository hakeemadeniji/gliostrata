# GlioStrata — System Architecture

## Eight-Layer Pipeline

```
L0  Ingestion from public data sources
    TCGA-GBM · CGGA · CPTAC · Neftel scRNA-seq · GLASS · TCIA · BraTS · DepMap · GDSC
    GDC REST API · direct download
        ↓
L1  Harmonization
    ComBat · Sample-ID canonicalization · QC gate
        ↓
L2  Tumor Profiling & Subtyping
    Verhaak subtype · Neftel state deconvolution · TME deconvolution · MGMT / radiomics
        ↓
L3  Unified Patient Feature Store
    Genomic + transcriptomic + proteomic + radiomic + clinical + longitudinal
        ↓
        ┌────────────────────────┬────────────────────────┐
        ↓                        ↓                        ↓
L4  Simulation Engine           L5  Outcome Prediction
    (mechanistic, per-patient)      (statistical / ML)
    PhysiCell tumor growth          DeepSurv Cox-PH NN
    NSC chemotaxis (SDF-1/VEGF/uPA)  Response classifier
    Payload diffusion + kill        Pseudoprogression classifier
        └────────────────────────┬────────────────────────┘
                                 ↓
L6  Therapeutic Configuration Optimizer
    BoTorch Bayesian optimization over:
    cell type · payload · dose · injection sites · schedule
                                 ↓
L7  Wet-Lab Handoff & Recalibration
    Auto-generated protocol (pandoc) · Predicted-outcome table
    REDCap results entry · Bayesian recalibration feedback loop
```

## Design Principles

1. **Layer independence** — each layer is runnable and publishable on its own.
2. **Public-data-only by default** — IRB-exempt secondary analysis.
3. **BU infrastructure target** — SCC HPC, Singularity, Snakemake.
4. **Wet-lab output as structured artifact** — protocol document + predicted-outcome table.

This is the **target** architecture. Layer numbers match the code and `dvc.yaml`
(L0 ingest, L1 harmonize, L2 profile, L3 feature store, L4 simulation, L5
prediction, L6 optimization, L7 handoff). For what is implemented today see
`docs/roadmap.md`; for what is reproducible, `docs/reproducibility.md`.

See `docs/protocol.md` for the full research protocol.
