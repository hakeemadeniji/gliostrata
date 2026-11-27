.PHONY: help install dev-install workflow-install test test-fast lint lint-strict format \
        clean clean-fixture synthetic pipeline pipeline-fixture dag dag-svg \
        ingest harmonize-real features-real cibersortx verify-partition \
        dvc-init dvc-dag dvc-repro ci

# Everything below is config-driven so a reviewer can run the whole pipeline
# against a tiny offline fixture without reading the CLI docs. See
# docs/reproducibility.md for what is actually reproducible today.

SNAKEFILE := workflow/Snakefile
CORES     ?= 2

# The offline fixture: synthetic L0, no network, everything under .test/.
# Keep the roots repo-relative — Snakemake's metadata store on Windows cannot
# handle long absolute out-of-tree paths.
FIXTURE_CONFIG := source=synthetic cohort=fixture n_synthetic=12 \
                  with_expression=false verify_partition=false \
                  raw_root=.test/raw interim_root=.test/interim \
                  processed_root=.test/processed results_root=.test/results

COHORT    ?= data/raw/tcga_gbm/clinical_cohort.parquet
# v2 is the CURRENT partition (594 patients). v1 is a superseded historical
# record of the 583-patient cohort and must never be re-frozen.
MANIFEST  ?= data/partitions/tcga_gbm_clinical_v2.json

help:
	@echo "GlioStrata development commands:"
	@echo ""
	@echo " Setup"
	@echo "  make install           - install core package"
	@echo "  make dev-install       - install with dev dependencies"
	@echo "  make workflow-install  - install with snakemake + dvc (the 'workflow' extra)"
	@echo ""
	@echo " Quality"
	@echo "  make test              - run the full pytest suite"
	@echo "  make test-fast         - run the suite without slow/heavy tests"
	@echo "  make lint              - ruff, project config (currently ~148 known findings)"
	@echo "  make lint-strict       - ruff, errors-only set that CI blocks on"
	@echo "  make format            - black + ruff --fix"
	@echo "  make ci                - everything CI runs, locally"
	@echo ""
	@echo " Pipeline (offline)"
	@echo "  make synthetic         - generate a synthetic GBM cohort (50 patients)"
	@echo "  make pipeline          - in-memory L0->L7 dry run on synthetic data"
	@echo "  make pipeline-fixture  - full Snakemake run on a 12-patient fixture -> .test/"
	@echo "  make dag               - print the real-cohort Snakemake DAG (executes nothing)"
	@echo "  make dag-svg           - render the DAG to dag.svg (needs graphviz)"
	@echo ""
	@echo " Pipeline (real cohort -- NETWORK: GDC, ~1.2 GB for expression)"
	@echo "  make ingest            - L0: download TCGA-GBM clinical + expression"
	@echo "  make harmonize-real    - L1 on the real clinical cohort"
	@echo "  make features-real     - L0->L3 on the real cohort"
	@echo "  make cibersortx        - export the CIBERSORTx mixture file (manual upload)"
	@echo ""
	@echo " Guardrails"
	@echo "  make verify-partition  - check the frozen train/test split still reproduces"
	@echo ""
	@echo " DVC"
	@echo "  make dvc-init          - initialise .dvc/ (read docs/reproducibility.md first)"
	@echo "  make dvc-dag           - print the DVC stage graph"
	@echo "  make dvc-repro         - reproduce the DVC pipeline (NETWORK + large)"
	@echo ""
	@echo " Cleanup"
	@echo "  make clean             - remove build artifacts and caches"
	@echo "  make clean-fixture     - remove .test/ and .snakemake/"

# ── Setup ─────────────────────────────────────────────────────────────────
install:
	pip install -e .

dev-install:
	pip install -e ".[dev]"

workflow-install:
	pip install -e ".[workflow]"

# ── Quality ───────────────────────────────────────────────────────────────
test:
	pytest

test-fast:
	pytest -m "not slow and not heavy"

lint:
	ruff check src tests scripts

# The subset CI blocks on: syntax errors, invalid comparisons, misplaced
# statements, undefined names. Passes clean today.
lint-strict:
	ruff check --select E9,F63,F7,F82 src tests scripts

format:
	black src tests
	ruff check --fix src tests

ci: lint-strict test dag pipeline-fixture verify-partition

# ── Pipeline (offline) ────────────────────────────────────────────────────
synthetic:
	gliostrata synthetic-data --n 50 --out data/interim/synthetic_cohort.parquet

pipeline:
	gliostrata pipeline-dry-run --config configs/default.yaml

pipeline-fixture:
	snakemake --cores $(CORES) --snakefile $(SNAKEFILE) --config $(FIXTURE_CONFIG)

dag:
	snakemake -n --cores 1 --snakefile $(SNAKEFILE)

dag-svg:
	snakemake --dag --snakefile $(SNAKEFILE) | dot -Tsvg > dag.svg
	@echo "Wrote dag.svg"

# ── Pipeline (real cohort — hits the GDC REST API) ────────────────────────
ingest:
	snakemake --cores $(CORES) --snakefile $(SNAKEFILE) ingest

harmonize-real:
	snakemake --cores $(CORES) --snakefile $(SNAKEFILE) \
		data/interim/tcga_gbm_harmonized.parquet

features-real:
	snakemake --cores $(CORES) --snakefile $(SNAKEFILE) features

# Manual round trip to the CIBERSORTx web tool: this only builds the upload
# file. Nothing in the pipeline consumes the fractions it returns yet.
cibersortx:
	snakemake --cores $(CORES) --snakefile $(SNAKEFILE) cibersortx

# ── Guardrails ────────────────────────────────────────────────────────────
# A sealed partition is a pre-registration commitment. Run this before every
# evaluation, and never "fix" a failure by re-freezing the manifest.
verify-partition:
	gliostrata freeze-partition --input $(COHORT) --manifest $(MANIFEST) --verify

# ── DVC ───────────────────────────────────────────────────────────────────
dvc-init:
	@echo "This writes .dvc/ and installs git hooks. Read the 'DVC' section of"
	@echo "docs/reproducibility.md first, including the storage implications."
	dvc init

dvc-dag:
	dvc dag

dvc-repro:
	dvc repro

# ── Cleanup ───────────────────────────────────────────────────────────────
clean:
	rm -rf build/ dist/ *.egg-info src/*.egg-info
	rm -rf .pytest_cache .mypy_cache .ruff_cache
	find . -type d -name __pycache__ -exec rm -rf {} +

clean-fixture:
	rm -rf .test/ .snakemake/ dag.svg
