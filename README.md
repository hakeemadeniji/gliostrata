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

## Status — v0.1.0 (scaffold)

This is a **first-draft scaffold**. The repository structure, CLI, tests, and synthetic-
data pipeline all work end-to-end. The actual scientific implementations (TCGA ingestion,
Verhaak classifier training, PhysiCell simulation, BoTorch optimization, DeepSurv survival
model) are stubs that return synthetic data — clearly marked in each module's docstring
and in `docs/roadmap.md`.

| Layer | Module | Status |
|-------|--------|--------|
| L0 — Public data sources | `gliostrata.data.ingest` | Stub (synthetic cohort generator works) |
| L1 — Harmonization | `gliostrata.data.harmonize`, `qc` | Stub |
| L2 — Tumor profiling | `gliostrata.profiling.*` | Stub (subtype, states, TME) |
| L3 — Feature store | `gliostrata.features.store` | **Working** |
| L4 — Simulation | `gliostrata.simulation.tumor`, `nsc` | Stub (synthetic curves) |
| L5 — Outcome prediction | `gliostrata.prediction.survival` | Stub (Cox-PH linear) |
| L6 — Optimizer | `gliostrata.optimization.config_space`, `bayesian` | **Working** (random search; BoTorch TODO) |
| L7 — Wet-lab handoff | `gliostrata.handoff.protocol` | **Working** (JSON protocol generator) |

**62 tests pass** on a clean install. Heavy optional dependencies (scanpy, lifelines,
botorch, PhysiCell) are NOT required.

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
| `gliostrata ingest --source <tcga\|cgga\|cptac\|synthetic> --out <path>` | L0 — ingest from a source (stub for real sources). |
| `gliostrata harmonize --input <path> --out <path>` | L1 — harmonize ingested data (stub). |
| `gliostrata profile --input <path> --out <path>` | L2 — tumor profiling (stub). |
| `gliostrata build-feature-store --input <path> --out <path>` | L3 — assemble unified feature store. |
| `gliostrata optimize --input <path> --out <dir> [--n-top 5] [--seed 42]` | L6 — Bayesian configuration optimization. |
| `gliostrata generate-protocol --input <configs.parquet> --feature-store <fs.parquet> --out <dir>` | L7 — wet-lab handoff protocol generation. |

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
├── requirements.txt               ← core deps (pip)
├── requirements-dev.txt           ← dev deps (pytest, black, ruff, mypy)
├── environment.yml                ← conda environment
├── Makefile                       ← common dev commands (make test, make install, ...)
├── dvc.yaml                       ← DVC pipeline definition
├── LICENSE                        ← MIT
├── .gitignore
├── .env.example                   ← env var template (dbGaP key, etc.)
├── src/gliostrata/                ← the Python package
│   ├── __init__.py
│   ├── __main__.py                ← allows `python -m gliostrata`
│   ├── cli.py                     ← Click-based CLI
│   ├── config.py                  ← YAML config loading + typed Config dataclass
│   ├── data/                      ← L0-L1
│   │   ├── ingest.py              ← L0: download from public sources (stub)
│   │   ├── harmonize.py           ← L1: harmonize + batch correct (stub)
│   │   └── qc.py                  ← L1: quality control gates (stub)
│   ├── profiling/                 ← L2
│   │   ├── subtype.py             ← Verhaak subtype classifier (stub)
│   │   ├── states.py              ← Neftel single-cell state deconvolution (stub)
│   │   └── tme.py                 ← Tumor microenvironment deconvolution (stub)
│   ├── features/                  ← L3
│   │   └── store.py               ← Unified patient feature store (working)
│   ├── simulation/                ← L4
│   │   ├── tumor.py               ← PhysiCell tumor growth wrapper (stub)
│   │   └── nsc.py                 ← NSC migration + payload kinetics (stub)
│   ├── prediction/                ← L5
│   │   └── survival.py            ← DeepSurv survival model (stub: Cox-PH linear)
│   ├── optimization/              ← L6
│   │   ├── config_space.py        ← Configuration space + constraints (working)
│   │   └── bayesian.py            ← Bayesian optimizer (working: random search; BoTorch TODO)
│   ├── handoff/                   ← L7
│   │   └── protocol.py            ← Wet-lab protocol document generator (working)
│   └── utils/
│       ├── io.py                  ← Parquet/JSON I/O helpers
│       └── stats.py               ← Bootstrap CI, CCC, RMST, BH-FDR
├── tests/                         ← pytest suite (62 tests, all pass)
│   ├── conftest.py                ← shared fixtures
│   ├── test_cli.py                ← CLI smoke tests
│   ├── test_data.py               ← L0-L1 tests
│   ├── test_profiling.py          ← L2 tests
│   ├── test_features.py           ← L3 tests
│   ├── test_simulation.py         ← L4-L5 tests
│   └── test_optimization.py       ← L6-L7 tests + end-to-end pipeline test
├── configs/                       ← YAML configs for each layer
│   ├── default.yaml               ← master config + pre-registered thresholds
│   ├── data/                      ← per-source configs (tcga, cgga)
│   ├── profiling/                 ← subtype, states
│   ├── simulation/                ← tumor, nsc
│   ├── prediction/                ← survival
│   └── optimization/              ← bayesian
├── workflow/                      ← Snakemake workflow (scaffold)
│   ├── Snakefile
│   └── rules/
├── scripts/                       ← standalone helper scripts
│   ├── generate_synthetic_data.py
│   └── run_full_pipeline.py
├── docs/
│   ├── protocol.md                ← full research protocol (markdown version)
│   ├── architecture.md            ← 8-layer architecture overview
│   └── roadmap.md                 ← implementation status + recommended order
├── notebooks/                     ← exploration notebooks
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
includes the **pre-registered success thresholds** (T1–T5) that are committed to the
Open Science Framework before any model touches the frozen test partition. Do not
modify these thresholds without filing a new OSF version.

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
make lint          # ruff check
make format        # black + ruff --fix
make synthetic     # generate synthetic cohort
make pipeline      # run pipeline dry-run
make clean         # remove build artifacts
```

---

## Data sources (all public, IRB-exempt secondary analysis)

| Source | Modality | Cohort Size | Access |
|--------|----------|-------------|--------|
| TCGA-GBM | Bulk WES, RNA-seq, methylation, clinical | ~600 | Open (GDC) |
| CGGA | Bulk RNA-seq, methylation, clinical | >1,000 | Open (registration) |
| CPTAC-GBM | Proteomics | ~100 | Open (GDC) |
| Neftel et al. 2019 | scRNA-seq | ~7,000 cells | Open (GEO) |
| Darmanis et al. 2017 | scRNA-seq | ~3,600 cells | Open (GEO) |
| GLASS | Longitudinal bulk | >800 samples | Open (application) |
| TCIA | Radiomics | ~260 patients | Open |
| BraTS | Segmented MRI | >1,200 scans | Open |
| DepMap / GDSC / CTRP | Drug response | >1,000 lines | Open |

See `docs/protocol.md` §4 for the full data catalog and harmonization strategy.

---

## Roadmap

See `docs/roadmap.md` for the recommended implementation order. Summary:

1. L0 ingestion — TCGA via `TCGAbiolinks` (R, called via `subprocess`)
2. L1 harmonization — `ComBat-seq` for batch correction
3. L2 subtype classifier — random forest ensemble on Verhaak signature genes
4. L2 state deconvolution — `CIBERSORTx` (web tool) → import results
5. L3 feature store — already working; just wire up real L2 outputs
6. L5 survival model — `lifelines.CoxPHFitter` first, then DeepSurv
7. L4 simulation — install PhysiCell, wrap with `pyMCDS`
8. L6 Bayesian optimizer — swap random search for `botorch` acquisition
9. L7 protocol generator — convert JSON output to `.docx` via `python-docx`

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
