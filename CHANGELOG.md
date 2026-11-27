# Changelog

All notable changes to GlioStrata will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added — pre-registration and the real data layers
- Pre-registration machinery: checksum-frozen threshold registry
  (`evaluation/thresholds.py`), frozen partitions with drift detection
  (`evaluation/partitions.py`), and the `freeze-partition` CLI command,
  including the `--all-test` seal for external cohorts.
- `docs/preregistration.md` (OSF pre-registration, not yet deposited) and
  `docs/osf_final_fields.md` (paste-ready OSF text).
- L0 ingestion: TCGA-GBM bulk RNA-seq, IDH mutation calls and MGMT methylation
  calls (`scripts/download_tcga_*.py`); CGGA RNA-seq ingestion
  (`scripts/download_cgga.py`) and the external-partition eligibility listing
  (`scripts/cgga_external_eligibility.py`).
- L2 from bulk expression, not wired into the pipeline: signature-based Verhaak
  subtype scoring (`profiling/verhaak.py`) and marker-based state/TME scoring
  (`profiling/deconvolve.py`). `scripts/join_molecular_calls.py` and
  `scripts/export_cibersortx.py`.
- L5: fitted Cox proportional-hazards model (`prediction/cox.py`), its T2
  external validation on CGGA (`scripts/validate_cox_cgga.py`) and the IDH
  sensitivity analyses (`scripts/idh_sensitivity_*.py`).
- T2's Brier row (score at 18 months): `scripts/validate_brier_cgga.py`, evaluated once on
  a second read of CGGA and not met (0.2617, 95% CI 0.2391 to 0.2829, against 0.25).
- `docs/project_record/`: a full project record (overview, findings, work log, commit
  history and future plan). Every number in it was cross-checked against its source.
- Read discipline for the held-out cohort: an append-only access ledger
  (`data/access_ledger.jsonl`, with every access to CGGA to date back-filled and
  sourced), a read guard (`evaluation/read_guard.py`) that refuses a value read
  until the registration is deposited and, where one applies, the partition is
  sealed and the analysis not already read, and `gliostrata access-ledger` to
  print the ledger.
- `configs/registration.yaml`: the machine-readable registration state. The OSF
  DOI is `null` until the deposit; a test keeps this file and the
  pre-registration header in agreement.

- L6: Gaussian-process / expected-improvement optimizer
  (`optimization/bayes_opt.py`).
- Snakemake workflow, DVC pipeline definition and GitHub Actions CI.

### Changed
- The TCGA-GBM clinical cohort is 594 patients (was 583): only contradictory
  diagnoses are dropped. Partition `tcga_gbm_clinical_v2` is the current sealed
  split; `v1` is kept unmodified as a superseded historical record.
- Pre-registration corrected before deposit: the T1 harmonisation rule
  (per-sample z-scoring), two-sided Clopper-Pearson bounds, the T4 tie rule and
  the cohort facts.
- `requests` is declared as a dependency; the download scripts import it.
- Ruff settings moved to `[tool.ruff.lint]`, so ruff no longer warns.
- `gliostrata ingest` and `data/ingest.py` point to the real download scripts
  instead of the roadmap.
- README, `docs/roadmap.md`, `docs/reproducibility.md` and
  `docs/mgmt_idh_calling.md` brought back in line with the code.
- CI partition guardrail: a manifest is active because `ACTIVE_MANIFESTS` lists
  it, not because it is the one manifest that is not v1. A new sealed manifest
  is added to that list in the commit that adds it, so the external partition
  will not be mislabelled "superseded".

### Fixed
- `test_argmax_matches_input_subtype` computed an argmax and asserted nothing
  about it; it now asserts the match rate.
- Lint backlog: unused imports, unsorted import blocks, an `== False`
  comparison and an unused variable.

### Added — repository scaffold (v0.1.0)
- Project structure with `src/gliostrata` package layout (L0–L7 modules).
- Click-based CLI (`gliostrata` command) with smoke-test, synthetic-data,
  pipeline-dry-run, and per-layer stub commands.
- Pytest test suite with synthetic-data fixtures; all tests pass on a clean install.
- Configuration files (`configs/*.yaml`) for every layer.
- Snakemake workflow scaffold (`workflow/Snakefile` + `rules/`).
- Documentation: `README.md`, `docs/protocol.md`, `docs/architecture.md`.
- DVC pipeline definition (`dvc.yaml`).
- MIT license.

### Notes
- L0–L7 modules are **stubs** with realistic function signatures and synthetic-data
  implementations where feasible. The actual scientific implementations
  (TCGA ingestion, Verhaak classifier, PhysiCell simulation, BoTorch optimizer,
  DeepSurv survival model, protocol generator) are TODO — see `docs/roadmap.md`.
- Heavy optional dependencies (scanpy, lifelines, botorch, physiCell) are NOT
  required to run the scaffold or the tests.
