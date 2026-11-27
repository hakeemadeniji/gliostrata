# GlioStrata — Implementation Roadmap

Status as of 2026-11-27 (stage 3b). Four statuses are used:

- **Working / Implemented** — does what it says, with tests.
- **Implemented, not wired** — real, tested code that the CLI pipeline does not
  call yet. Wiring it in is a change to `src/gliostrata/cli.py` and `dvc.yaml`.
- **Implemented, wired (real mode)** — called by `gliostrata build-feature-store --mode real`
  on TCGA; the synthetic pipeline still uses the stub of the same layer.
- **Stub** — returns synthetic or schema-only output and says so in its docstring.

## Status by layer

| Layer | Module | Status | Notes |
|-------|--------|--------|-------|
| L0 | `scripts/download_tcga_{clinical,expression,mutations,methylation}.py` | Implemented | Real GDC REST ingestion; Snakemake rules and `dvc.yaml` stages |
| L0 | `scripts/download_cgga.py` | Implemented | Pre-registered external cohort, RNA-seq batches; the microarray batch is not ingested yet |
| L0 | `data/ingest.py` | Stub | `generate_synthetic_cohort` works; the `download_*` functions raise |
| L1 | `data/harmonize.py` | Stub | Schema normalisation + QC rows; ComBat / harmony calls TODO |
| L1 | `data/clinical.py` | Working (TCGA) | `gliostrata harmonize-clinical`: clinical table with owned columns, a separate label table, a reason code for every missing value; nothing imputed |
| L1 | `data/qc.py` | Stub | Purity / consistency checks TODO |
| L2 | `profiling/centroid.py` | Implemented, wired (real mode) | Single-sample nearest centroid over the published Verhaak 840-gene signature (831 genes resolved); two frozen centroid sets with checksums. The per-patient subtype. Measured, not validated |
| L2 | `profiling/verhaak.py` | Implemented | Cohort-relative signature scoring; descriptive cohort distributions only |
| L2 | `profiling/deconvolve.py` | Implemented, wired (real mode) | Marker-based state and TME scoring from bulk expression; cohort-relative by default |
| L2 | `profiling/{subtype,states,tme}.py` | Stub | Synthetic calls; this is what `gliostrata profile` runs |
| L3 | `features/{store,contract,views,build}.py` | Working | `gliostrata build-feature-store --mode real\|synthetic`. The store holds every cohort patient (594 on TCGA) and never filters; outcomes in a separate label table; named views profiling (283), transferable model (570) and molecular model (298) with per-split flow counts |
| L4 | `simulation/tumor.py`, `simulation/nsc.py` | Stub | PhysiCell wrapper and chemotaxis solver TODO; synthetic curves provided |
| L5 | `prediction/cox.py` | Implemented, not wired | Fitted Cox-PH. External C-index and Brier score on CGGA miss the T2 thresholds (`docs/cgga_validation.md`) |
| L5 | `prediction/survival.py` | Stub | Synthetic hazard model; this is what the CLI uses |
| L6 | `optimization/config_space.py` | Working | Configuration space fully defined |
| L6 | `optimization/bayesian.py` | Working | Ranked enumeration |
| L6 | `optimization/bayes_opt.py` | Implemented, not wired | Gaussian-process surrogate + expected improvement; BoTorch backend TODO |
| L7 | `handoff/protocol.py` | Working | JSON protocol documents; not bench-executable |
| — | `evaluation/{thresholds,partitions,read_guard}.py` | Working | Pre-registered threshold registry, frozen partitions, access ledger and read guard |

## Recommended implementation order

1. **L0 ingestion** — done for TCGA (direct GDC REST API in `scripts/`, not
   TCGAbiolinks) and CGGA RNA-seq. CPTAC is not started.
2. **L1 harmonization** — ComBat-seq for batch correction (open).
3. **L2 subtype classifier** — single-sample nearest centroid over the published
   Verhaak 840-gene signature (`profiling/centroid.py`), wired into
   `build-feature-store --mode real`. Not the protocol's random forest + 1D-CNN.
   Its T1 evaluation on CGGA is not run.
4. **L2 state deconvolution** — `CIBERSORTx` (web tool) → import results.
   `profiling/deconvolve.py` is the marker-based alternative; the CIBERSORTx
   signature-matrix builder is not implemented.
5. **L3 feature store** — working. `gliostrata build-feature-store --mode real` runs
   L2 and L3 on the real TCGA inputs (all 594 patients in the store; three named
   views with their flow counts). The DVC stage and the Snakemake rule (with
   `source=tcga`) run the same command; `scripts/build_tcga_feature_store.py` is a
   thin wrapper over it. The synthetic mode stays the default for tests.
6. **L5 survival model** — implemented as a hand-written Cox-PH
   (`prediction/cox.py`, Breslow partial likelihood). DeepSurv is open.
7. **L4 simulation** — install PhysiCell, wrap with `pyMCDS`.
8. **L6 Bayesian optimizer** — a GP / expected-improvement optimizer exists
   (`optimization/bayes_opt.py`); the `botorch` acquisition is open.
9. **L7 protocol generator** — convert JSON output to `.docx` via `python-docx`.

## What the tests cover

1113 tests in 34 modules (1112 passed, 1 skipped) run without the optional dependencies. They
pass without network access and without heavy dependencies (scanpy, lifelines,
botorch and PhysiCell are not required). Tests that read the real TCGA
data skip when it is absent; no test reads CGGA; one workflow test skips when Snakemake is not
installed.

- `test_cli.py` — CLI, including `freeze-partition` and `access-ledger`
- `test_data.py`, `test_profiling.py`, `test_features.py`, `test_simulation.py`,
  `test_optimization.py` — the synthetic L0→L7 path, including the end-to-end pipeline test
- `test_evaluation.py` — threshold registry (checksum), partition freeze / verify / drift
- `test_workflow.py` — the Snakemake, DVC and CI definitions agree with the code
- `test_verhaak.py`, `test_deconvolve.py` — expression-based L2 modules
- `test_cox.py`, `test_validate_cox_cgga.py`, `test_validate_brier_cgga.py` — Cox model
  and the two T2 validation scripts (C-index, Brier score)
- `test_tcga_download.py`, `test_tcga_expression.py`, `test_tcga_mutations.py`,
  `test_tcga_methylation.py` — the GDC ingest scripts
- `test_cgga.py`, `test_cgga_eligibility.py` — CGGA ingest and external-partition eligibility
- `test_join_molecular_calls.py`, `test_idh_sensitivity_cox.py`,
  `test_idh_sensitivity_subtype.py` — molecular joins and the sensitivity analysis
- `test_export_cibersortx.py`, `test_bayes_opt.py`, `test_stats.py`, `test_read_guard.py` —
  export script, GP optimizer, statistics helpers, and the access ledger with its read gates
