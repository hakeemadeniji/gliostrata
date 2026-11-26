"""CLI smoke tests — verify the gliostrata command works end-to-end."""
from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from gliostrata import __version__
from gliostrata.cli import cli


def test_cli_help():
    """`gliostrata --help` returns 0 and lists commands."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "GlioStrata" in result.output
    assert "synthetic-data" in result.output
    assert "pipeline-dry-run" in result.output


def test_cli_version():
    """`gliostrata version` prints the package version."""
    runner = CliRunner()
    result = runner.invoke(cli, ["version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_cli_hello():
    """`gliostrata hello` smoke test."""
    runner = CliRunner()
    result = runner.invoke(cli, ["hello"])
    assert result.exit_code == 0
    assert "ready" in result.output.lower()


def test_cli_synthetic_data(tmp_path):
    """`gliostrata synthetic-data` writes a Parquet file."""
    runner = CliRunner()
    out = tmp_path / "cohort.parquet"
    result = runner.invoke(cli, ["synthetic-data", "--n", "10", "--out", str(out)])
    assert result.exit_code == 0, result.output
    assert out.exists()
    assert "Generated 10 synthetic patients" in result.output


def test_cli_pipeline_dry_run():
    """`gliostrata pipeline-dry-run` executes the full L0→L7 pipeline on synthetic data."""
    runner = CliRunner()
    result = runner.invoke(cli, ["pipeline-dry-run", "--n", "10"])
    assert result.exit_code == 0, result.output
    # Check that every layer was touched
    for layer_marker in ["L0", "L1", "L2", "L3", "L5", "L6", "L7"]:
        assert layer_marker in result.output
    assert "Pipeline complete" in result.output


def test_cli_pipeline_dry_run_without_out_writes_nothing(tmp_path):
    """Without --out the dry run is purely in-memory: exit 0, no files written."""
    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=tmp_path) as cwd:
        result = runner.invoke(cli, ["pipeline-dry-run", "--n", "10"])
        assert result.exit_code == 0, result.output
        assert "no files were written" in result.output
        assert "--out" in result.output
        leftovers = [p for p in Path(cwd).rglob("*") if p.is_file()]
        assert leftovers == [], leftovers


def test_cli_pipeline_dry_run_with_out_writes_artifacts(tmp_path):
    """With --out the dry run writes the pipeline artifacts to that directory."""
    runner = CliRunner()
    n = 10
    out_dir = tmp_path / "results"
    result = runner.invoke(cli, ["pipeline-dry-run", "--n", str(n), "--out", str(out_dir)])
    assert result.exit_code == 0, result.output
    for name in ["cohort.parquet", "harmonized.parquet", "feature_store.parquet",
                 "survival_predictions.parquet", "ranked_configurations.parquet"]:
        assert (out_dir / name).exists(), f"missing {name}"
    protocols = sorted((out_dir / "protocols").glob("protocol_*.json"))
    assert len(protocols) == n
    assert str(out_dir) in result.output
    assert "Pipeline complete" in result.output


def test_cli_ingest_synthetic(tmp_path):
    """`gliostrata ingest --source synthetic` works; real sources are stubs."""
    runner = CliRunner()
    out_dir = tmp_path / "raw"
    result = runner.invoke(cli, ["ingest", "--source", "synthetic",
                                  "--out", str(out_dir), "--n", "5"])
    assert result.exit_code == 0
    assert (out_dir / "synthetic_cohort.parquet").exists()


def test_cli_ingest_tcga_stub(tmp_path):
    """`gliostrata ingest --source tcga` reports it's a stub (does NOT crash)."""
    runner = CliRunner()
    result = runner.invoke(cli, ["ingest", "--source", "tcga", "--out", str(tmp_path)])
    # The stub message should be shown, not raise an exception
    assert "STUB" in result.output


def test_cli_manual_layer_by_layer(tmp_path):
    """Regression: running L0→L7 layer-by-layer via the CLI (parquet round-trips)
    must preserve upstream columns (e.g. `sex`) so build-feature-store succeeds."""
    import pandas as pd

    runner = CliRunner()
    n = 30
    cohort = tmp_path / "cohort.parquet"
    harmonized = tmp_path / "harmonized.parquet"
    profiled = tmp_path / "profiled.parquet"
    features = tmp_path / "features.parquet"
    configs_dir = tmp_path / "configs"
    protocols_dir = tmp_path / "protocols"

    steps = [
        ["synthetic-data", "--n", str(n), "--seed", "42", "--out", str(cohort)],
        ["harmonize", "--input", str(cohort), "--out", str(harmonized)],
        ["profile", "--input", str(harmonized), "--out", str(profiled)],
        ["build-feature-store", "--input", str(profiled), "--out", str(features)],
        ["optimize", "--input", str(features), "--out", str(configs_dir),
         "--n-top", "5", "--seed", "42"],
    ]
    for args in steps:
        result = runner.invoke(cli, args)
        assert result.exit_code == 0, f"{args[0]} failed:\n{result.output}\n{result.exception!r}"

    # Each intermediate layer's schema must be a superset of its input's.
    cohort_cols = set(pd.read_parquet(cohort).columns)
    harmonized_cols = set(pd.read_parquet(harmonized).columns)
    profiled_cols = set(pd.read_parquet(profiled).columns)
    assert cohort_cols <= harmonized_cols
    assert harmonized_cols <= profiled_cols

    assert features.exists()
    config_files = list(configs_dir.glob("*.parquet"))
    assert len(config_files) == 1

    result = runner.invoke(cli, [
        "generate-protocol", "--input", str(config_files[0]),
        "--feature-store", str(features), "--out", str(protocols_dir),
    ])
    assert result.exit_code == 0, f"generate-protocol failed:\n{result.output}\n{result.exception!r}"
    assert len(list(protocols_dir.glob("*.json"))) == n
