"""Click-based CLI for GlioStrata.

Commands:
    gliostrata hello                  — smoke test
    gliostrata version                — show version
    gliostrata synthetic-data         — generate a synthetic GBM cohort
    gliostrata pipeline-dry-run       — run end-to-end pipeline on synthetic data
    gliostrata ingest                 — L0: ingest from a public source (stub)
    gliostrata harmonize              — L1: harmonize + QC (stub)
    gliostrata profile                — L2: tumor profiling (stub)
    gliostrata build-feature-store    — L3: assemble unified feature store
    gliostrata optimize               — L6: Bayesian configuration optimization
    gliostrata generate-protocol      — L7: wet-lab protocol document generation
"""
from __future__ import annotations

from pathlib import Path

import click

from gliostrata import __version__
from gliostrata.config import default_config_path, load_config


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(version=__version__, prog_name="gliostrata")
def cli() -> None:
    """GlioStrata — patient-specific stem-cell-delivered therapy optimizer for GBM."""


@cli.command()
def hello() -> None:
    """Smoke test — print a greeting."""
    click.echo(f"GlioStrata v{__version__} — ready. Try `gliostrata --help`.")


@cli.command()
def version() -> None:
    """Print the package version."""
    click.echo(__version__)


@cli.command("synthetic-data")
@click.option("--n", default=50, show_default=True, type=int, help="Number of synthetic patients.")
@click.option("--seed", default=42, show_default=True, type=int, help="Random seed.")
@click.option("--out", "out_path", required=True, type=click.Path(),
              help="Output Parquet path (e.g. data/interim/synthetic_cohort.parquet).")
def synthetic_data(n: int, seed: int, out_path: str) -> None:
    """Generate a synthetic GBM cohort for development and testing."""
    from gliostrata.data.ingest import generate_synthetic_cohort
    from gliostrata.utils.io import write_parquet

    df = generate_synthetic_cohort(n=n, seed=seed)
    written = write_parquet(df, out_path)
    click.echo(f"Generated {len(df)} synthetic patients -> {written}")
    click.echo(f"Columns: {list(df.columns)}")


@cli.command("pipeline-dry-run")
@click.option("--config", "config_path", default=None, type=click.Path(),
              help="Path to YAML config (defaults to configs/default.yaml).")
@click.option("--n", default=20, show_default=True, type=int,
              help="Number of synthetic patients for the dry run.")
@click.option("--out", "out_path", default=None, type=click.Path(),
              help="Optional output directory for the pipeline artifacts. "
                   "If omitted, nothing is written (pure in-memory run).")
def pipeline_dry_run(config_path: str | None, n: int, out_path: str | None) -> None:
    """Run the full L0->L7 pipeline on synthetic data and report results.

    This is the end-to-end smoke test of the scaffold. It:
      1. Generates a synthetic cohort (L0).
      2. Harmonizes (L1 stub).
      3. Profiles tumors (L2 stubs).
      4. Builds the feature store (L3).
      5. Predicts survival (L5 stub).
      6. Optimizes therapeutic configurations (L6).
      7. Emits a wet-lab protocol (L7).

    By default this is a pure in-memory run: no files are written. Pass
    --out DIR to persist the pipeline artifacts to DIR.
    """
    cfg = load_config(config_path) if config_path else load_config(default_config_path())
    click.echo(f"Loaded config: {cfg.project_name} v{cfg.version}, seed={cfg.seed}")
    click.echo("Running end-to-end pipeline on synthetic data...\n")

    from gliostrata.data.ingest import generate_synthetic_cohort
    from gliostrata.data.harmonize import harmonize
    from gliostrata.profiling.subtype import classify_subtypes
    from gliostrata.profiling.states import deconvolve_states
    from gliostrata.profiling.tme import deconvolve_tme
    from gliostrata.features.store import build_feature_store
    from gliostrata.prediction.survival import predict_survival
    from gliostrata.optimization.bayesian import optimize_configurations
    from gliostrata.handoff.protocol import generate_protocols

    # L0
    cohort = generate_synthetic_cohort(n=n, seed=cfg.seed)
    click.echo(f"  L0  ingestion       : {len(cohort)} patients (synthetic)")
    # L1
    harmonized = harmonize(cohort)
    click.echo(f"  L1  harmonization   : {len(harmonized)} samples after QC")
    # L2
    subtypes = classify_subtypes(harmonized, seed=cfg.seed)
    states = deconvolve_states(harmonized, seed=cfg.seed)
    tme = deconvolve_tme(harmonized, seed=cfg.seed)
    click.echo(f"  L2  profiling       : {subtypes['subtype'].nunique()} subtypes, "
               f"states+TME fractions computed")
    # L3
    feature_store = build_feature_store(harmonized, subtypes, states, tme)
    click.echo(f"  L3  feature store   : {feature_store.shape[0]} patients × "
               f"{feature_store.shape[1]} features")
    # L5
    survival_pred = predict_survival(feature_store, seed=cfg.seed)
    click.echo(f"  L5  survival pred   : median predicted OS = "
               f"{survival_pred['predicted_survival_months'].median():.1f} months")
    # L6
    configs = optimize_configurations(feature_store, survival_pred, n_top=3, seed=cfg.seed)
    click.echo(f"  L6  optimizer       : {len(configs)} ranked configs per patient "
               f"({configs['rank'].nunique()} ranks)")
    # L7
    protocols = generate_protocols(configs, feature_store, survival_pred)
    click.echo(f"  L7  handoff         : {len(protocols)} protocol documents generated")

    if out_path is None:
        click.echo("\nPipeline complete. In-memory run only -- no files were written. "
                   "Run with --out DIR to write the artifacts.")
        return

    from gliostrata.utils.io import ensure_dir, write_json, write_parquet
    from gliostrata.handoff.protocol import protocols_to_json

    out_dir = ensure_dir(out_path)
    write_parquet(cohort, out_dir / "cohort.parquet")
    write_parquet(harmonized, out_dir / "harmonized.parquet")
    write_parquet(feature_store, out_dir / "feature_store.parquet")
    write_parquet(survival_pred, out_dir / "survival_predictions.parquet")
    write_parquet(configs, out_dir / "ranked_configurations.parquet")
    protocols_dir = ensure_dir(out_dir / "protocols")
    protocols_data = protocols_to_json(protocols)
    for pid, proto in protocols_data.items():
        write_json(proto, protocols_dir / f"protocol_{pid}.json")
    click.echo(f"\nPipeline complete. Wrote 5 Parquet files and "
               f"{len(protocols_data)} protocol JSON documents -> {out_dir}")


@cli.command()
@click.option("--source", type=click.Choice(["tcga", "cgga", "cptac", "synthetic"]),
              default="synthetic", show_default=True,
              help="Data source to ingest.")
@click.option("--out", "out_path", required=True, type=click.Path(),
              help="Output directory or file path.")
@click.option("--n", default=50, show_default=True, type=int,
              help="Number of samples (only used for --source synthetic).")
def ingest(source: str, out_path: str, n: int) -> None:
    """L0 — ingest data from a public source.

    NOTE: TCGA/CGGA/CPTAC ingestion is a STUB. Only --source synthetic works.
    Real ingestion requires TCGAbiolinks (R) for TCGA, direct download for CGGA.
    See docs/roadmap.md.
    """
    from gliostrata.utils.io import write_parquet

    if source == "synthetic":
        from gliostrata.data.ingest import generate_synthetic_cohort
        df = generate_synthetic_cohort(n=n)
        out = Path(out_path)
        if out.is_dir() or not out.suffix:
            out = out / "synthetic_cohort.parquet"
        write_parquet(df, out)
        click.echo(f"Ingested synthetic cohort ({len(df)} patients) -> {out}")
    else:
        click.echo(f"STUB: ingestion of {source} is not yet implemented. "
                   f"Use --source synthetic for now. See docs/roadmap.md.")


@cli.command("harmonize")
@click.option("--input", "in_path", required=True, type=click.Path(exists=True))
@click.option("--out", "out_path", required=True, type=click.Path())
def harmonize_cmd(in_path: str, out_path: str) -> None:
    """L1 — harmonize ingested data (stub)."""
    from gliostrata.utils.io import read_parquet, write_parquet
    from gliostrata.data.harmonize import harmonize

    df = read_parquet(in_path)
    out = harmonize(df)
    write_parquet(out, out_path)
    click.echo(f"Harmonized {len(out)} samples -> {out_path}")


@cli.command()
@click.option("--input", "in_path", required=True, type=click.Path(exists=True))
@click.option("--out", "out_path", required=True, type=click.Path())
def profile(in_path: str, out_path: str) -> None:
    """L2 — tumor profiling (subtype + state + TME) (stub).

    Writes the input cohort enriched with L2 columns (subtype probabilities,
    state/TME fractions, profiling metadata). The output schema is always a
    superset of the input schema — L3 assembly happens in `build-feature-store`.
    """
    from gliostrata.utils.io import read_parquet, write_parquet
    from gliostrata.profiling.subtype import classify_subtypes
    from gliostrata.profiling.states import deconvolve_states
    from gliostrata.profiling.tme import deconvolve_tme

    df = read_parquet(in_path)
    profiled = df.copy()
    for layer_out in (classify_subtypes(df), deconvolve_states(df), deconvolve_tme(df)):
        new_cols = [c for c in layer_out.columns if c not in profiled.columns]
        if new_cols:
            profiled = profiled.merge(
                layer_out[["patient_id", *new_cols]], on="patient_id", how="left"
            )
    write_parquet(profiled, out_path)
    click.echo(f"Profiled {len(profiled)} patients -> {out_path}")


@cli.command("build-feature-store")
@click.option("--input", "in_path", required=True, type=click.Path(exists=True))
@click.option("--out", "out_path", required=True, type=click.Path())
def build_feature_store_cmd(in_path: str, out_path: str) -> None:
    """L3 — assemble the unified patient feature store."""
    from gliostrata.utils.io import read_parquet, write_parquet
    from gliostrata.profiling.subtype import classify_subtypes
    from gliostrata.profiling.states import deconvolve_states
    from gliostrata.profiling.tme import deconvolve_tme
    from gliostrata.features.store import build_feature_store

    df = read_parquet(in_path)
    feature_store = build_feature_store(
        df,
        classify_subtypes(df),
        deconvolve_states(df),
        deconvolve_tme(df),
    )
    write_parquet(feature_store, out_path)
    click.echo(f"Feature store: {feature_store.shape} -> {out_path}")


@cli.command()
@click.option("--input", "in_path", required=True, type=click.Path(exists=True),
              help="Feature store Parquet file.")
@click.option("--out", "out_path", required=True, type=click.Path(),
              help="Output directory for ranked configurations.")
@click.option("--n-top", default=5, show_default=True, type=int,
              help="Number of top configurations to return per patient.")
@click.option("--seed", default=42, show_default=True, type=int)
def optimize(in_path: str, out_path: str, n_top: int, seed: int) -> None:
    """L6 — Bayesian configuration optimization.

    Generates a ranked set of therapeutic configurations per patient.
    """
    from gliostrata.utils.io import read_parquet, ensure_dir, write_parquet
    from gliostrata.prediction.survival import predict_survival
    from gliostrata.optimization.bayesian import optimize_configurations

    feature_store = read_parquet(in_path)
    survival_pred = predict_survival(feature_store, seed=seed)
    configs = optimize_configurations(feature_store, survival_pred, n_top=n_top, seed=seed)
    out_dir = ensure_dir(out_path)
    write_parquet(configs, out_dir / "ranked_configurations.parquet")
    click.echo(f"Optimized configs for {configs['patient_id'].nunique()} patients "
               f"-> {out_dir / 'ranked_configurations.parquet'}")


@cli.command("generate-protocol")
@click.option("--input", "in_path", required=True, type=click.Path(exists=True),
              help="Ranked configurations Parquet (output of `gliostrata optimize`).")
@click.option("--feature-store", "fs_path", required=True, type=click.Path(exists=True))
@click.option("--out", "out_path", required=True, type=click.Path(),
              help="Output directory for protocol JSON documents.")
def generate_protocol_cmd(in_path: str, fs_path: str, out_path: str) -> None:
    """L7 — generate wet-lab handoff protocol documents."""
    from gliostrata.utils.io import read_parquet, ensure_dir, write_json
    from gliostrata.prediction.survival import predict_survival
    from gliostrata.handoff.protocol import generate_protocols, protocols_to_json

    configs = read_parquet(in_path)
    feature_store = read_parquet(fs_path)
    survival_pred = predict_survival(feature_store)
    protocols = generate_protocols(configs, feature_store, survival_pred)
    out_dir = ensure_dir(out_path)
    protocols_data = protocols_to_json(protocols)
    for pid, proto in protocols_data.items():
        write_json(proto, out_dir / f"protocol_{pid}.json")
    click.echo(f"Generated {len(protocols_data)} protocol documents -> {out_dir}")


if __name__ == "__main__":
    cli()
