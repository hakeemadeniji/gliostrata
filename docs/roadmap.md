# GlioStrata — Implementation Roadmap

This scaffold is **v0.1.0**: structure, CLI, tests, and stub modules.
The actual scientific implementations are TODO. Use this document to track
what Claude Code (or you) should implement next.

## Status by layer

| Layer | Module | Status | Notes |
|-------|--------|--------|-------|
| L0 | `data/ingest.py` | Stub | Function signatures defined; download not implemented |
| L1 | `data/harmonize.py` | Stub | ComBat / harmony calls TODO |
| L1 | `data/qc.py` | Stub | Purity / consistency checks TODO |
| L2 | `profiling/subtype.py` | Stub | Returns synthetic subtype calls |
| L2 | `profiling/states.py` | Stub | Returns synthetic state fractions |
| L2 | `profiling/tme.py` | Stub | Returns synthetic TME fractions |
| L3 | `features/store.py` | Working | Assembles feature store from L2 outputs (synthetic OK) |
| L4 | `simulation/tumor.py` | Stub | PhysiCell wrapper TODO; synthetic growth curve provided |
| L4 | `simulation/nsc.py` | Stub | Chemotaxis solver TODO; synthetic migration provided |
| L5 | `prediction/survival.py` | Stub | DeepSurv TODO; synthetic hazard model provided |
| L6 | `optimization/config_space.py` | Working | Configuration space fully defined |
| L6 | `optimization/bayesian.py` | Working | Random-search optimizer (BoTorch wrapper is TODO) |
| L7 | `handoff/protocol.py` | Working | Generates JSON protocol document |

## Recommended implementation order

1. **L0 ingestion** — start with TCGA via `TCGAbiolinks` (R, called via `subprocess`)
2. **L1 harmonization** — ComBat-seq for batch correction
3. **L2 subtype classifier** — random forest ensemble on Verhaak signature genes
4. **L2 state deconvolution** — `CIBERSORTx` (web tool) → import results
5. **L3 feature store** — already working; just wire up real L2 outputs
6. **L5 survival model** — `lifelines.CoxPHFitter` first, then DeepSurv
7. **L4 simulation** — install PhysiCell, wrap with `pyMCDS`
8. **L6 Bayesian optimizer** — swap random search for `botorch` acquisition
9. **L7 protocol generator** — convert JSON output to `.docx` via `python-docx`

## What the tests cover

- `test_cli.py` — CLI smoke tests (gliostrata --help, hello, version, synthetic-data)
- `test_data.py` — synthetic data generation + harmonize stub
- `test_profiling.py` — subtype/states/tme stubs return valid outputs
- `test_features.py` — feature store assembly from L2 outputs
- `test_optimization.py` — config space + optimizer end-to-end on synthetic data
- `test_cli.py::test_pipeline_dry_run` — full pipeline on synthetic data

All tests pass without heavy dependencies (scanpy, lifelines, botorch not required).
