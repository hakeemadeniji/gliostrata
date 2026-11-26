.PHONY: help install dev-install test lint format clean synthetic pipeline

help:
	@echo "GlioStrata development commands:"
	@echo "  make install      - install core package"
	@echo "  make dev-install  - install with dev dependencies"
	@echo "  make test         - run pytest test suite"
	@echo "  make lint         - run ruff linter"
	@echo "  make format       - run black + ruff --fix"
	@echo "  make synthetic    - generate a synthetic GBM cohort (50 patients)"
	@echo "  make pipeline     - run end-to-end pipeline on synthetic data"
	@echo "  make clean        - remove build artifacts and caches"

install:
	pip install -e .

dev-install:
	pip install -e ".[dev]"

test:
	pytest

lint:
	ruff check src tests

format:
	black src tests
	ruff check --fix src tests

synthetic:
	gliostrata synthetic-data --n 50 --out data/interim/synthetic_cohort.parquet

pipeline:
	gliostrata pipeline-dry-run --config configs/default.yaml

clean:
	rm -rf build/ dist/ *.egg-info src/*.egg-info
	rm -rf .pytest_cache .mypy_cache .ruff_cache
	find . -type d -name __pycache__ -exec rm -rf {} +
