"""Click-based CLI for GlioStrata.

Commands:
    gliostrata hello                  — smoke test
    gliostrata version                — show version
    gliostrata synthetic-data         — generate a synthetic GBM cohort
    gliostrata pipeline-dry-run       — run end-to-end pipeline on synthetic data
    gliostrata ingest                 — L0: ingest from a public source (stub)
    gliostrata harmonize              — L1: harmonize + QC (stub)
    gliostrata profile                — L2: tumor profiling (stub)
    gliostrata harmonize-clinical     — L1: clinical table with owned columns + label table
    gliostrata build-feature-store    — L2 + L3: feature store (--mode synthetic | real)
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

    from gliostrata.data.harmonize import harmonize
    from gliostrata.data.ingest import generate_synthetic_cohort
    from gliostrata.features.store import build_synthetic_store, stub_layer_frame
    from gliostrata.handoff.protocol import generate_protocols
    from gliostrata.optimization.bayesian import optimize_configurations
    from gliostrata.prediction.survival import predict_survival
    from gliostrata.profiling.states import deconvolve_states
    from gliostrata.profiling.subtype import classify_subtypes
    from gliostrata.profiling.tme import deconvolve_tme

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
    # L3 (synthetic mode), then the one-hot layout the L5-L7 stubs consume
    store = build_synthetic_store(harmonized, subtypes, states, tme)
    feature_store = stub_layer_frame(store)
    click.echo(f"  L3  feature store   : {feature_store.shape[0]} patients × "
               f"{feature_store.shape[1]} features (synthetic)")
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

    from gliostrata.handoff.protocol import protocols_to_json
    from gliostrata.utils.io import ensure_dir, write_json, write_parquet

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

    NOTE: only --source synthetic is implemented here. TCGA and CGGA ingestion
    are real, but they live in standalone scripts that talk to the GDC and CGGA
    directly: scripts/download_tcga_*.py and scripts/download_cgga.py.
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
        click.echo(f"STUB: ingestion of {source} is not implemented in the CLI. "
                   f"Use scripts/download_tcga_*.py (TCGA) or scripts/download_cgga.py "
                   f"(CGGA) directly, or --source synthetic.")


@cli.command("harmonize")
@click.option("--input", "in_path", required=True, type=click.Path(exists=True))
@click.option("--out", "out_path", required=True, type=click.Path())
def harmonize_cmd(in_path: str, out_path: str) -> None:
    """L1 — harmonize ingested data (stub)."""
    from gliostrata.data.harmonize import harmonize
    from gliostrata.utils.io import read_parquet, write_parquet

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
    from gliostrata.profiling.states import deconvolve_states
    from gliostrata.profiling.subtype import classify_subtypes
    from gliostrata.profiling.tme import deconvolve_tme
    from gliostrata.utils.io import read_parquet, write_parquet

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


@cli.command("harmonize-clinical")
@click.option("--cohort", "cohort_path", default="data/raw/tcga_gbm/clinical_cohort.parquet",
              show_default=True, type=click.Path(), help="L0 clinical cohort.")
@click.option("--raw", "raw_path", default="data/raw/tcga_gbm/clinical_raw.parquet",
              show_default=True, type=click.Path(),
              help="Raw GDC clinical records, used only to trace why a value is missing.")
@click.option("--clinical-out", default="data/interim/tcga_gbm_clinical_harmonized.parquet",
              show_default=True, type=click.Path())
@click.option("--labels-out", default="data/interim/tcga_gbm_labels.parquet",
              show_default=True, type=click.Path())
def harmonize_clinical_cmd(cohort_path: str, raw_path: str, clinical_out: str,
                           labels_out: str) -> None:
    """L1 — clinical table with explicitly owned columns, and a separate label table.

    The clinical table carries no L2 or molecular name; a column under such a
    name that holds a value stops the command. Outcomes go to the label table.
    Every missing value gets a reason code. Nothing is imputed. Refuses any
    path containing "cgga".

    \b
    Also writes a provenance sidecar next to the clinical table,
    <clinical stem>_provenance.json (by default
    data/interim/tcga_gbm_clinical_harmonized_provenance.json): the SHA-256 of
    the L0 inputs and of both outputs. build-feature-store --mode real reads it
    and stops if it is missing or no longer matches the clinical table.
    """
    from gliostrata.features.build import run_harmonize_clinical
    from gliostrata.features.contract import FeatureStoreError

    try:
        r = run_harmonize_clinical(cohort_path, raw_path, clinical_out, labels_out)
    except (FeatureStoreError, ValueError) as exc:
        raise click.ClickException(str(exc)) from None
    click.echo(f"Harmonized clinical table: {r['n_patients']} patients -> {clinical_out}")
    click.echo(f"Label table -> {labels_out}")
    click.echo(f"Provenance (input and output SHA-256) -> {r['outputs']['provenance']}")
    for table in ("clinical_reasons", "label_reasons"):
        for col, counts in r[table].items():
            click.echo(f"  missing {col:<22} {counts}")


@cli.command("build-feature-store")
@click.option("--mode", type=click.Choice(["synthetic", "real"]), default="synthetic",
              show_default=True,
              help="synthetic: offline simulator + L2 stubs (tests, fixture runs). "
                   "real: TCGA expression through the real L2 (single-sample subtype, "
                   "state and TME scores) and L3.")
@click.option("--input", "in_path", default=None, type=click.Path(),
              help="synthetic mode: harmonized synthetic cohort Parquet.")
@click.option("--clinical", default=None, type=click.Path(),
              help="real mode: L1 clinical table (harmonize-clinical).")
@click.option("--labels", default=None, type=click.Path(),
              help="real mode: L1 label table (harmonize-clinical).")
@click.option("--expression", default=None, type=click.Path(),
              help="real mode: gene x sample TPM matrix.")
@click.option("--annotation", default=None, type=click.Path(),
              help="real mode: gene annotation of the TPM matrix.")
@click.option("--mgmt-calls", default=None, type=click.Path(),
              help="real mode: MGMT call table.")
@click.option("--idh-calls", default=None, type=click.Path(),
              help="real mode: IDH call table.")
@click.option("--mgmt-manifest", default=None, type=click.Path(),
              help="real mode: methylation GDC file manifest (which patients were assayed).")
@click.option("--idh-manifest", default=None, type=click.Path(),
              help="real mode: mutation GDC file manifest (which patients were sequenced).")
@click.option("--centroids", default=None, type=click.Path(),
              help="real mode: frozen subtype centroid set (JSON); only "
                   "verhaak840_tcga_v2train is accepted.")
@click.option("--partition", default=None, type=click.Path(),
              help="real mode: frozen development partition, for per-split counts.")
@click.option("--out", "out_path", default=None, type=click.Path(),
              help="Features table (Parquet). Real-mode default "
                   "data/processed/tcga_gbm_feature_store.parquet.")
@click.option("--labels-out", default=None, type=click.Path(),
              help="Label table (default: <out stem>_labels.parquet).")
@click.option("--missing-out", default=None, type=click.Path(),
              help="Missing-value reason table (default: <out stem>_missing.parquet).")
@click.option("--columns-out", default=None, type=click.Path(),
              help="Column metadata JSON (default: <out stem>_columns.json).")
@click.option("--report", "report_path", default=None, type=click.Path(),
              help="Build report JSON: inputs with SHA-256, views and their flow counts.")
def build_feature_store_cmd(mode: str, in_path: str | None, clinical: str | None,
                            labels: str | None, expression: str | None,
                            annotation: str | None, mgmt_calls: str | None,
                            idh_calls: str | None, mgmt_manifest: str | None,
                            idh_manifest: str | None, centroids: str | None,
                            partition: str | None, out_path: str | None,
                            labels_out: str | None, missing_out: str | None,
                            columns_out: str | None, report_path: str | None) -> None:
    """L2 + L3 — run profiling and assemble the feature store.

    \b
    The store holds every cohort patient, one row each, and never filters.
    It is written as four files: the features table, a separate label table
    (outcomes are labels, never features), a missing-value reason table and
    a column-metadata JSON (role, source layer, source file SHA-256, class).
    Model-ready populations are named views; the report gives their flow
    counts, split by the frozen development partition in real mode.

    \b
    --mode synthetic (default): --input HARMONIZED --out STORE. Offline.
    --mode real: TCGA only. Every input defaults to its repository path.
    A stub L2 output, a missing centroid set, or any path containing "cgga"
    stops the command before anything is written.
    """
    from gliostrata.features.build import REAL_DEFAULTS, run_real, run_synthetic
    from gliostrata.features.contract import FeatureStoreError

    try:
        if mode == "synthetic":
            if not in_path or not out_path:
                raise click.UsageError("--mode synthetic needs --input and --out")
            store = run_synthetic(in_path, out_path, labels_out=labels_out,
                                  missing_out=missing_out, columns_out=columns_out,
                                  report=report_path)
            click.echo(f"Feature store (SYNTHETIC): {len(store.features)} patients x "
                       f"{store.features.shape[1]} columns -> {out_path}")
            return
        if in_path:
            raise click.UsageError("--input is for --mode synthetic; real mode takes "
                                   "--clinical, --labels, --expression, ...")
        d = REAL_DEFAULTS
        r = run_real(
            clinical=clinical or d["clinical"], labels=labels or d["labels"],
            expression=expression or d["expression"],
            annotation=annotation or d["annotation"],
            mgmt_calls=mgmt_calls or d["mgmt_calls"], idh_calls=idh_calls or d["idh_calls"],
            mgmt_manifest=mgmt_manifest or d["mgmt_manifest"],
            idh_manifest=idh_manifest or d["idh_manifest"],
            centroids=centroids or d["centroids"], partition=partition or d["partition"],
            out=out_path or d["out"],
            labels_out=labels_out or (None if out_path else d["labels_out"]),
            missing_out=missing_out or (None if out_path else d["missing_out"]),
            columns_out=columns_out or (None if out_path else d["columns_out"]),
            report=report_path or d["report"],
        )
    except FeatureStoreError as exc:
        raise click.ClickException(str(exc)) from None

    click.echo("Feature store (REAL, development cohort TCGA). Measured, not validated.")
    click.echo(f"  patients in the store              {r['n_patients']}")
    c = r["counts"]
    click.echo(f"  expression samples                 {c['n_expression_samples']} "
               f"({c['n_expression_samples_outside_cohort']} outside the cohort, not scored)")
    click.echo(f"  centroid set                       {r['centroid_set']['set_id']}")
    for name, view in r["views"].items():
        click.echo(f"\n  view {name}: {view['n']} patients")
        for step in view["flow"]:
            split = ""
            if "n_remaining_train" in step:
                split = (f"  (train {step['n_remaining_train']}, "
                         f"test {step['n_remaining_test']})")
            click.echo(f"    {step['step']:<26} -{step['n_excluded']:<4} "
                       f"-> {step['n_remaining']}{split}")
    click.echo("\nWrote:")
    for rec in r["outputs"].values():
        click.echo(f"  {rec['path']}  sha256 {rec['sha256']}")
    click.echo(f"  {report_path or REAL_DEFAULTS['report']}")


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
    from gliostrata.features.contract import FeatureStoreError
    from gliostrata.features.store import stub_layer_frame
    from gliostrata.optimization.bayesian import optimize_configurations
    from gliostrata.prediction.survival import predict_survival
    from gliostrata.utils.io import ensure_dir, read_parquet, write_parquet

    try:
        feature_store = stub_layer_frame(read_parquet(in_path))
    except FeatureStoreError as exc:
        raise click.ClickException(f"L6 is a stub that cannot use this store: {exc}") from None
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
    from gliostrata.features.contract import FeatureStoreError
    from gliostrata.features.store import stub_layer_frame
    from gliostrata.handoff.protocol import generate_protocols, protocols_to_json
    from gliostrata.prediction.survival import predict_survival
    from gliostrata.utils.io import ensure_dir, read_parquet, write_json

    configs = read_parquet(in_path)
    try:
        feature_store = stub_layer_frame(read_parquet(fs_path))
    except FeatureStoreError as exc:
        raise click.ClickException(f"L7 is a stub that cannot use this store: {exc}") from None
    survival_pred = predict_survival(feature_store)
    protocols = generate_protocols(configs, feature_store, survival_pred)
    out_dir = ensure_dir(out_path)
    protocols_data = protocols_to_json(protocols)
    for pid, proto in protocols_data.items():
        write_json(proto, out_dir / f"protocol_{pid}.json")
    click.echo(f"Generated {len(protocols_data)} protocol documents -> {out_dir}")


if __name__ == "__main__":
    cli()


# ---------------------------------------------------------------------------
# Appended below the `__main__` guard on purpose: this file is edited
# concurrently by several agents, so this command is added as a pure append
# with no change to anything above. Click registers it at import time
# regardless of position, so the placement is cosmetic only.
# ---------------------------------------------------------------------------


@cli.command("freeze-partition")
@click.option("--input", "in_path", required=True, type=click.Path(exists=True),
              help="Cohort Parquet file to partition (must have a patient_id column).")
@click.option("--manifest", "manifest_path", required=True, type=click.Path(),
              help="Path to the partition manifest JSON (e.g. data/partitions/cgga_v1.json).")
@click.option("--name", "partition_name", default=None,
              help="Partition name recorded in the manifest. Defaults to the manifest stem.")
@click.option("--seed", default=20251126, show_default=True, type=int,
              help="Seed baked into the patient-id hash. Changing it is a "
                   "pre-registration amendment, not a re-run.")
@click.option("--test-fraction", default=0.30, show_default=True, type=float,
              help="Target fraction of patients held out.")
@click.option("--all-test", "all_test", is_flag=True,
              help="Hold out EVERY patient (test fraction 1.0): an empty training "
                   "split and a test split holding every patient id. This is the "
                   "seal for an external cohort, which has no training data. "
                   "Exclusive with --test-fraction and --stratify-by.")
@click.option("--stratify-by", default=None,
              help="Optional column to stratify on (e.g. subtype). Must be constant per patient.")
@click.option("--verify", "verify_only", is_flag=True,
              help="Do not write anything: check that the cohort still reproduces the "
                   "frozen manifest and report precisely what drifted.")
@click.option("--no-seal", is_flag=True,
              help="Write an UNSEALED draft manifest that may be freely overwritten. "
                   "Never use this for a partition any model will be trained against.")
@click.option("--force", is_flag=True,
              help="Overwrite a sealed manifest. Only legitimate as a documented "
                   "pre-registration amendment with every model retrained from scratch.")
@click.option("--notes", default="", help="Free-text provenance note stored in the manifest.")
@click.pass_context
def freeze_partition_cmd(ctx: click.Context, in_path: str, manifest_path: str,
                         partition_name: str | None,
                         seed: int, test_fraction: float, all_test: bool,
                         stratify_by: str | None,
                         verify_only: bool, no_seal: bool, force: bool, notes: str) -> None:
    """Freeze (or verify) a train/test partition by seeded hash of patient_id.

    The split is a pure function of (seed, patient_id), so it is stable under
    row reordering, cohort growth and library upgrades -- unlike a shuffle,
    which silently re-assigns patients whenever upstream filtering changes the
    row count. Once a model has trained on a partition, moving a patient from
    train to test destroys that patient's value as an external test forever,
    which is why a sealed manifest refuses to be overwritten.

    Verify before every evaluation run:

        gliostrata freeze-partition --input COHORT --manifest M.json --verify

    Seal an external cohort, which has no training split, with every patient
    held out (test fraction 1.0, empty training list):

        gliostrata freeze-partition --input COHORT --manifest M.json --all-test
    """
    from gliostrata.evaluation import (
        PartitionSealedError,
        freeze_partition,
        read_manifest,
        verify_partition,
    )
    from gliostrata.utils.io import read_parquet

    if all_test:
        # Refuse contradictory instructions before anything is read or written.
        # click.core.ParameterSource exists from click 8.0; the top-level alias
        # does not exist in click 8.1, the declared floor.
        from click.core import ParameterSource

        for name, flag in (("test_fraction", "--test-fraction"),
                           ("stratify_by", "--stratify-by")):
            if ctx.get_parameter_source(name) == ParameterSource.COMMANDLINE:
                raise click.UsageError(
                    "--all-test holds out every patient, so there is nothing to "
                    f"split; do not also pass {flag}."
                )
        test_fraction = 1.0

    cohort = read_parquet(in_path)

    if verify_only:
        report = verify_partition(read_manifest(manifest_path), cohort, cohort_source=in_path)
        click.echo(report.describe())
        if not report.ok:
            raise click.ClickException(
                f"Partition drift detected ({report.n_drifted} id(s) affected). "
                f"Do not evaluate against this partition until the drift is explained."
            )
        return

    name = partition_name or Path(manifest_path).stem
    try:
        result, manifest, written = freeze_partition(
            cohort,
            cohort_source=in_path,
            manifest_path=manifest_path,
            partition_name=name,
            seed=seed,
            test_fraction=test_fraction,
            stratify_by=stratify_by,
            seal=not no_seal,
            force=force,
            notes=notes,
        )
    except PartitionSealedError as exc:
        raise click.ClickException(str(exc)) from None

    click.echo(f"Partition '{name}' {'sealed' if not no_seal else 'written (UNSEALED draft)'} "
               f"-> {written}")
    click.echo(f"  patients      : {result.n_total} "
               f"({result.n_train} train / {result.n_test} test, "
               f"achieved test fraction {result.n_test / result.n_total:.4f})")
    click.echo(f"  seed          : {manifest.seed}")
    click.echo(f"  cohort sha256 : {manifest.cohort_sha256}")
    click.echo(f"  test ids sha256: {manifest.test_ids_sha256}")
    if manifest.stratum_balance:
        click.echo("  stratum balance:")
        for stratum, row in sorted(manifest.stratum_balance.items()):
            click.echo(f"    {stratum:<20} n={row['n']:<6} test={row['n_test']:<6} "
                       f"({row['test_fraction']:.3f})")


# ---------------------------------------------------------------------------
# Access ledger. Appended, like freeze-partition above, as a pure append.
# ---------------------------------------------------------------------------


@cli.command("access-ledger")
@click.option("--ledger", "ledger_path", default=None, type=click.Path(exists=True),
              help="Ledger file to print. Defaults to data/access_ledger.jsonl.")
def access_ledger_cmd(ledger_path: str | None) -> None:
    """Print every recorded access to the held-out cohort.

    The ledger (data/access_ledger.jsonl) is append-only and version-controlled.
    Every access is in it: schema looks, dry runs, checksum-only accesses and
    value reads. Run this before any evaluation on the external cohort.

    A read that inspects values is refused by gliostrata.evaluation.read_guard
    until the registration is deposited, the partition is sealed or the
    analysis's registered rule names a frozen scoring-set manifest, and every
    input's checksum is the one on disk. Entries that are observations or notes
    are marked "not an access". Verify a partition with:

        gliostrata freeze-partition --input COHORT --manifest M.json --verify
    """
    from gliostrata.evaluation.read_guard import (
        LEDGER_PATH,
        format_ledger,
        is_deposited,
        load_ledger,
        registration_state,
    )

    click.echo(format_ledger(load_ledger(ledger_path or LEDGER_PATH)))
    state = registration_state()
    if is_deposited(state):
        click.echo(f"Registration: deposited, DOI {state['osf_registration_doi']}.")
    else:
        click.echo("Registration: NOT deposited. A read that inspects values is refused.")
