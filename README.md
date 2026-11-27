# GlioStrata

**A patient-specific computational framework for optimizing stem-cell-delivered therapy in glioblastoma.**

GlioStrata is an open-source research framework that integrates public multi-omic GBM
data (TCGA, CGGA, CPTAC, single-cell atlases, GLASS, TCIA) with mechanistic tumor
simulation and Bayesian configuration optimization to recommend, per patient, the
**cell type, payload, dose, injection sites, and schedule** most likely to extend
survival beyond the Stupp protocol baseline.

The project is designed as a Master's thesis protocol (Health Informatics, Boston
University) with a **computational-first** deliverable and an **optional wet-lab
extension** running in parallel. See `docs/protocol.md` for the full research protocol
and `docs/architecture.md` for the system architecture.

---

## Status

GlioStrata is run as a **pre-registered** study: the evaluation plan, its thresholds
and the frozen train/test partitions are fixed in this repository. The OSF
registration itself is drafted (`docs/preregistration.md`) but **not yet deposited**.

The pipeline runs end to end on synthetic data. The scientific layers are at
different stages, stated per layer:

| Layer | Module | Status |
|-------|--------|--------|
| L0 — Ingestion | `scripts/download_tcga_*.py`, `scripts/download_cgga.py` | **Working** — real GDC and CGGA downloads. `gliostrata ingest` itself only implements `--source synthetic`. |
| L1 — Harmonization | `gliostrata.data.harmonize`, `gliostrata.data.qc`, `gliostrata.data.clinical` | `harmonize` is a stub (schema normalisation and QC rows; no batch correction). `harmonize-clinical` is **working** on TCGA: a clinical table with explicitly owned columns and a separate label table, every missing value with a reason code, nothing imputed. |
| L2 — Tumor profiling | `gliostrata.profiling.centroid`, `gliostrata.profiling.deconvolve`, `gliostrata.profiling.verhaak` | **Implemented** from bulk expression; run on TCGA by `gliostrata build-feature-store --mode real`. Subtype: single-sample nearest centroid over the published Verhaak 840-gene signature (831 genes resolved), centroids frozen from the TCGA v2 training split. State and TME scores are still cohort-relative. On 100 simulated mixtures from the Neftel 2019 atlas the state scores recover the imposed proportions (T1's deconvolution row: mean per-state r 0.9521, 95% CI 0.9430 to 0.9589, against 0.75; run before the OSF deposit and optimistic by design, see `docs/t1_deconvolution_spec.md`). Measured, not validated on tumours. What `gliostrata profile` runs (`subtype`, `states`, `tme`) is still a stub. |
| L3 — Feature store | `gliostrata.features.store`, `contract`, `views`, `build` | **Working**. Built from real data by `gliostrata build-feature-store --mode real`: all 594 TCGA patients, one row each, never filtered; outcomes in a separate label table; every column with its source and class, every missing value with a reason code. Model-ready populations are named views: profiling 283, transferable model 570, molecular model 298 (`docs/project_record/06_evidence_reports.md`, stage 3b). |
| L4 — Simulation | `gliostrata.simulation.tumor`, `nsc` | Stub (synthetic curves) |
| L5 — Outcome prediction | `gliostrata.prediction.cox` | **Fitted Cox-PH model**, not wired into the pipeline. Neither its external C-index nor its Brier score on CGGA met the registered T2 thresholds (`docs/cgga_validation.md`). `gliostrata.prediction.survival`, which the CLI uses, is a stub. |
| L6 — Optimizer | `gliostrata.optimization.config_space`, `bayesian` | **Working** (ranked enumeration). `optimization.bayes_opt` adds a Gaussian-process optimizer; the BoTorch backend is TODO. |
| L7 — Wet-lab handoff | `gliostrata.handoff.protocol` | **Working** (JSON protocol generator) |
| — Evaluation | `gliostrata.evaluation` | **Working** — pre-registered threshold registry, frozen partitions, drift verification, access ledger and read guard |

The test suite has 1130 tests in 35 modules: 1129 pass and 1 is skipped (the skipped test needs Snakemake, an
optional dependency). It passes with all outbound sockets blocked. Tests that need the real
TCGA data skip when it is absent, and no test reads CGGA. Heavy optional dependencies
(scanpy, lifelines, botorch, PhysiCell) are NOT required.

---

## Quick start

### Prerequisites

- Python 3.10+
- `pip` (or `conda`/`mamba` — see `environment.yml`)

### Install

```bash
# Clone or unzip the repository, then:
cd gliostrata
pip install -e ".[dev]"
```

For a minimal install without dev tools:

```bash
pip install -e .
```

For the full stack (heavy optional deps for actual scientific implementation):

```bash
pip install -e ".[dev,heavy]"
```

### Verify the install

```bash
# Smoke test
gliostrata hello

# Show available commands
gliostrata --help

# Run the full L0→L7 pipeline on synthetic data
gliostrata pipeline-dry-run --n 20
```

### Run the test suite

```bash
pytest
```

---

## CLI commands

| Command | Description |
|---------|-------------|
| `gliostrata hello` | Smoke test — print greeting + version. |
| `gliostrata version` | Print package version. |
| `gliostrata synthetic-data --n 50 --out <path>` | Generate a synthetic GBM cohort (Parquet). |
| `gliostrata pipeline-dry-run [--n 20] [--config <path>]` | Run the full L0→L7 pipeline on synthetic data and print a summary. |
| `gliostrata ingest --source <tcga\|cgga\|cptac\|synthetic> --out <path>` | L0 — ingest from a source. Only `synthetic` is implemented here; for TCGA and CGGA use `scripts/download_tcga_*.py` and `scripts/download_cgga.py`. |
| `gliostrata harmonize --input <path> --out <path>` | L1 — harmonize ingested data (stub). |
| `gliostrata profile --input <path> --out <path>` | L2 — tumor profiling (stub; not the expression-based code, which `build-feature-store --mode real` runs). |
| `gliostrata harmonize-clinical [--cohort <path>] [--raw <path>]` | L1 — TCGA clinical table with owned columns (`data/interim/tcga_gbm_clinical_harmonized.parquet`) and a separate label table (`data/interim/tcga_gbm_labels.parquet`). Every missing value is traced to the raw GDC record for its reason code. Refuses any path containing "cgga". Also writes `data/interim/tcga_gbm_clinical_harmonized_provenance.json` (SHA-256 of inputs and outputs). |
| `gliostrata build-feature-store --mode synthetic --input <path> --out <path>` | L2 + L3 on synthetic data (the default mode; labelled synthetic in every file it writes). |
| `gliostrata build-feature-store --mode real` | L2 + L3 on TCGA: single-sample subtype against the frozen centroid set `verhaak840_tcga_v2train`, state and TME scores, molecular calls. Writes the store, its label table, its missing-value reason table, its column metadata and a report of the three views with per-split flow counts. Every input defaults to its repository path. A stub L2 output, a missing centroid set, a centroid set other than `verhaak840_tcga_v2train`, a missing or mismatched L1 provenance sidecar, or any path containing "cgga" stops it before anything is written. |
| `gliostrata optimize --input <path> --out <dir> [--n-top 5] [--seed 42]` | L6 — therapeutic configuration optimization (ranked enumeration). |
| `gliostrata generate-protocol --input <configs.parquet> --feature-store <fs.parquet> --out <dir>` | L7 — wet-lab handoff protocol generation. |
| `gliostrata freeze-partition --input <cohort.parquet> --manifest <path.json> [--verify]` | Seal a train/test partition, or verify a sealed one still reproduces. Run `--verify` before every evaluation. |
| `gliostrata access-ledger` | Print every recorded access to the held-out cohort, and whether the registration is deposited. |

### End-to-end example (synthetic data)

```bash
# 1. Generate a synthetic cohort
gliostrata synthetic-data --n 50 --out data/interim/synthetic_cohort.parquet

# 2. Profile tumors (L2)
gliostrata profile --input data/interim/synthetic_cohort.parquet \
                   --out data/processed/features.parquet

# 3. Optimize therapeutic configurations (L6)
gliostrata optimize --input data/processed/features.parquet \
                    --out data/results/configs/ \
                    --n-top 5 --seed 42

# 4. Generate wet-lab protocol documents (L7)
gliostrata generate-protocol \
    --input data/results/configs/ranked_configurations.parquet \
    --feature-store data/processed/features.parquet \
    --out data/results/protocols/

# Each patient now has a protocol_<patient_id>.json in data/results/protocols/
```

---

## Project structure

```
gliostrata/
├── README.md                      ← you are here
├── pyproject.toml                 ← package config + dependencies
├── requirements.txt / requirements-dev.txt / environment.yml
├── Makefile                       ← dev commands (make test, make ci, ...)
├── dvc.yaml                       ← DVC pipeline definition (not initialised; see docs/reproducibility.md)
├── .github/workflows/ci.yml       ← lint, tests, partition guardrail, workflow dry run
├── src/gliostrata/                ← the Python package
│   ├── cli.py                     ← Click-based CLI
│   ├── config.py                  ← YAML config loading + typed Config dataclass
│   ├── data/                      ← L0-L1
│   │   ├── ingest.py              ← synthetic cohort generator (download_* are stubs)
│   │   ├── harmonize.py           ← L1: schema harmonization (stub: no batch correction)
│   │   ├── qc.py                  ← L1: quality control gates (stub)
│   │   └── clinical.py            ← L1: TCGA clinical table with owned columns, label table, reason codes
│   ├── profiling/                 ← L2
│   │   ├── centroid.py            ← single-sample Verhaak 840-gene nearest centroid (the per-patient subtype)
│   │   ├── verhaak.py             ← cohort-relative Verhaak signature scoring (descriptive distributions only)
│   │   ├── deconvolve.py          ← marker-based cell-state / TME scoring
│   │   ├── subtype.py             ← stubs: what `gliostrata profile` runs
│   │   ├── states.py              ←   (synthetic calls)
│   │   └── tme.py                 ←   (synthetic calls)
│   ├── features/                  ← L3: store.py (assembly), contract.py (files, classes, reason codes), views.py (named views), build.py (real and synthetic runs)
│   ├── simulation/                ← L4 (stubs)
│   │   ├── tumor.py               ← PhysiCell growth wrapper
│   │   └── nsc.py                 ← NSC migration + payload kinetics
│   ├── prediction/                ← L5
│   │   ├── cox.py                 ← fitted Cox-PH model
│   │   └── survival.py            ← stub used by `gliostrata optimize`
│   ├── optimization/              ← L6
│   │   ├── config_space.py        ← configuration space + constraints
│   │   ├── bayesian.py            ← ranked enumeration
│   │   └── bayes_opt.py           ← Gaussian-process Bayesian optimization
│   ├── evaluation/                ← pre-registration machinery
│   │   ├── thresholds.py          ← pre-registered thresholds T1-T5 (checksum-frozen)
│   │   ├── read_guard.py          ← access ledger and read guard for the held-out cohort
│   │   └── partitions.py          ← deterministic, sealable partitions
│   ├── handoff/protocol.py        ← L7: wet-lab protocol document generator
│   └── utils/                     ← io.py, stats.py
├── scripts/                       ← standalone scripts
│   ├── download_tcga_{clinical,expression,mutations,methylation}.py   ← L0, network
│   ├── download_cgga.py           ← L0: external validation cohort, network
│   ├── join_molecular_calls.py    ← join IDH / MGMT calls onto the cohort
│   ├── build_tcga_feature_store.py  ← thin wrapper: `gliostrata build-feature-store --mode real`
│   ├── build_subtype_centroids.py, prepare_subtype_reference.py  ← the frozen 840-gene centroid sets
│   ├── subtype_scoring_set_dependence.py, compare_subtype_labels.py  ← descriptive subtype checks
│   ├── export_cibersortx.py       ← CIBERSORTx mixture file
│   ├── validate_cox_cgga.py       ← T2 external validation (C-index)
│   ├── validate_brier_cgga.py     ← T2 external validation (Brier score)
│   ├── idh_sensitivity_cox.py, idh_sensitivity_subtype.py
│   ├── cgga_external_eligibility.py  ← patients eligible for the sealed external partition
│   ├── generate_synthetic_data.py
│   └── run_full_pipeline.py
├── tests/                         ← 35 test modules, 1130 tests
├── configs/                       ← YAML configs (default.yaml holds the T1-T5 thresholds; registration.yaml the registration state; subtype_centroids/ the frozen centroid sets with their checksums)
├── workflow/                      ← Snakemake workflow (Snakefile + rules/)
├── data/partitions/               ← sealed partition manifests (in git)
├── data/access_ledger.jsonl        ← every access to the held-out cohort (in git, append-only)
├── docs/
│   ├── preregistration.md         ← the OSF pre-registration (not yet deposited)
│   ├── osf_final_fields.md        ← paste-ready OSF registration text
│   ├── project_record/           ← full project record: findings, work log, commits, plan
│   ├── protocol.md                ← full research protocol
│   ├── architecture.md, roadmap.md, reproducibility.md
│   └── cgga_ingestion.md, cgga_validation.md, mgmt_idh_calling.md
├── notebooks/
│   └── 01_data_exploration.ipynb
└── data/                          ← local data (gitignored; .gitkeep markers)
    ├── raw/                       ← L0 downloads
    ├── interim/                   ← L1 harmonized
    ├── processed/                 ← L3 feature store
    ├── external/                  ← reference data
    └── results/                   ← L6-L7 outputs
```

---

## Configuration

Layer configs are in `configs/`. The master config is `configs/default.yaml` and
includes the **pre-registered success thresholds** (T1–T5). Do not modify these
thresholds without filing a new version on the Open Science Framework (the deposit
is pending; see "Pre-registration" below).

`src/gliostrata/evaluation/thresholds.py` is the authoritative, checksum-frozen
registry of these thresholds (a test fails on any change to the table), and
`docs/preregistration.md` is the registered text.

Key thresholds (from `configs/default.yaml`):

```yaml
thresholds:
  t1_subtype_accuracy_cgga: 0.88       # 95% CI lower bound
  t2_survival_cindex_cgga: 0.65
  t3_tumor_growth_rmse_pct: 15.0
  t4_optimizer_beats_random_pct: 25.0
  t5_rmst_gain_months: 3.0
  t5_success_fraction: 0.50             # ≥ 50% of held-out with ≥ 3-month RMST gain
  t5_harm_fraction_max: 0.10            # ≤ 10% predicted harm
```

Environment variables (`.env`, copy from `.env.example`):

- `DBGAP_API_KEY` — only needed for TCGA controlled-access data
- `SCC_USERNAME`, `SCC_PROJECT_DIR` — only for BU Shared Computing Cluster
- `MLFLOW_TRACKING_URI` — only for experiment tracking
- `OSF_TOKEN` — only for protocol pre-registration

The scaffold runs end-to-end without ANY of these set.

---

## Development

```bash
# Install with dev deps
pip install -e ".[dev]"

# Run tests
pytest                                # all tests
pytest tests/test_optimization.py     # one file
pytest -k "personalization"           # by name pattern
pytest -m "not slow"                  # skip slow tests

# Lint + format
ruff check src tests
black src tests
ruff check --fix src tests

# Type check (optional)
mypy src/gliostrata
```

`Makefile` shortcuts:

```bash
make install       # pip install -e .
make dev-install   # pip install -e ".[dev]"
make test          # pytest
make lint             # ruff, full project config (advisory backlog)
make lint-strict      # ruff, errors only -- the subset CI blocks on
make verify-partition # check the frozen train/test split still reproduces
make ci               # everything CI runs, locally
make synthetic        # generate synthetic cohort
make pipeline         # run pipeline dry-run
make clean            # remove build artifacts
```

---

## Data sources (all public, IRB-exempt secondary analysis)

| Source | Modality | Cohort Size | Access |
|--------|----------|-------------|--------|
| TCGA-GBM | Bulk WES, RNA-seq, methylation, clinical | ~600 | Open (GDC) |
| CGGA | Bulk RNA-seq, methylation, clinical | >1,000 | Open (no registration) |
| CPTAC-GBM | Proteomics | ~100 | Open (GDC) |
| Neftel et al. 2019 | scRNA-seq | ~7,000 cells | Open (GEO) |
| Darmanis et al. 2017 | scRNA-seq | ~3,600 cells | Open (GEO) |
| GLASS | Longitudinal bulk | >800 samples | Open (application) |
| TCIA | Radiomics | ~260 patients | Open |
| BraTS | Segmented MRI | >1,200 scans | Open |
| DepMap / GDSC / CTRP | Drug response | >1,000 lines | Open |

See `docs/protocol.md` §4 for the full data catalog and harmonization strategy.

---

## Pre-registration

The evaluation plan is fixed before evaluation, and the repository enforces it:

- `docs/preregistration.md` is the pre-registration. `docs/osf_final_fields.md`
  holds the paste-ready OSF registration text.
- `configs/default.yaml` and `src/gliostrata/evaluation/thresholds.py` hold the
  thresholds T1–T5; a test fails if the threshold table changes.
- `data/partitions/*.json` are sealed train/test splits. `tcga_gbm_clinical_v2.json`
  (594 patients) is the current development partition; `tcga_gbm_clinical_v1.json`
  is a superseded historical record of the 583-patient cohort.

Rules that keep the registration meaningful:

- Run `make verify-partition` before every evaluation. Never re-freeze a sealed
  manifest to silence a drift warning; if the cohort changes, freeze a new,
  separately named partition.
- CGGA is the external validation cohort. Do not read it for model or feature
  selection (`scripts/download_cgga.py` says the same).
- Every access to CGGA is listed in `data/access_ledger.jsonl` (print it with
  `gliostrata access-ledger`). A script that would inspect values calls
  `gliostrata.evaluation.read_guard.record_access` first. The guard refuses the read
  until `configs/registration.yaml` carries the OSF DOI, the partition is sealed (or
  the scoring set is fixed by registered rule), and the analysis has not been read.
- Results so far: T2's C-index and Brier criteria are **not met** on CGGA
  (`docs/cgga_validation.md`). The IDH sensitivity analysis is in
  `docs/mgmt_idh_calling.md` §6.3.

---

## Roadmap

See `docs/roadmap.md` for the per-layer status and the recommended order. In short:
L0 ingestion is implemented for TCGA (clinical, RNA-seq, IDH, MGMT) and CGGA RNA-seq;
the expression-based L2 modules, the Cox model (L5) and the GP optimizer (L6) exist
but are not wired into the pipeline; L1 batch correction, L4 simulation, the BoTorch
backend and a bench-executable protocol are open.

---

## Citation

If you use GlioStrata in your research, cite:

> GlioStrata: A Patient-Specific Computational Framework for Optimizing
> Stem-Cell-Delivered Therapy in Glioblastoma. Master's Thesis Research Protocol,
> Boston University Health Informatics Program, v0.1.0, November 2026.

And see `docs/protocol.md` Appendix C for the primary-literature reference list.

---

## License

MIT — see [LICENSE](LICENSE).
