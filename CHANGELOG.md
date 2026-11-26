# Changelog

All notable changes to GlioStrata will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
