# ---------------------------------------------------------------------------
# workflow/rules/common.smk — configuration defaults and path helpers.
#
# Every path in the workflow is derived from the config dict, so the whole
# pipeline can be pointed at a small test fixture instead of the full cohort:
#
#     snakemake -n --cores 1 --snakefile workflow/Snakefile \
#         --config source=synthetic cohort=fixture n_synthetic=12 \
#                  raw_root=.test/raw interim_root=.test/interim \
#                  processed_root=.test/processed results_root=.test/results \
#                  partitions_root=.test/partitions with_expression=false
#
# Nothing here is science; it is pure plumbing.
# ---------------------------------------------------------------------------
from pathlib import PurePosixPath


def _as_bool(value) -> bool:
    """Coerce a --config value to bool. `--config x=false` arrives as the str 'false'."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


# ── Configuration defaults ─────────────────────────────────────────────────
# `source` selects the L0 entry point:
#   tcga      — real GDC ingestion via scripts/download_tcga_*.py (network)
#   synthetic — `gliostrata synthetic-data` (offline; the CI / fixture path)
config.setdefault("source", "tcga")
config.setdefault("cohort", "tcga_gbm")
config.setdefault("gdc_project", "TCGA-GBM")

config.setdefault("raw_root", "data/raw")
config.setdefault("interim_root", "data/interim")
config.setdefault("processed_root", "data/processed")
config.setdefault("results_root", "data/results")
config.setdefault("partitions_root", "data/partitions")

config.setdefault("seed", 42)
config.setdefault("n_synthetic", 50)
config.setdefault("n_top", 5)

# Expression ingest is a large (~1.2 GB cached) network download. `expression_limit`
# caps the number of GDC files pulled; 0 means the full cohort.
config.setdefault("with_expression", config["source"] == "tcga")
config.setdefault("expression_limit", 0)

# Partition manifest to verify before any evaluation. Verification is a hard
# scientific gate, not a build step: see workflow/rules/partition.smk.
#
# v2 is the CURRENT partition (594 patients). v1 (583 patients) is retained
# unmodified as a superseded historical record of the pre-filter-revision
# cohort and must never be re-frozen — see docs/reproducibility.md §4.
config.setdefault("partition_manifest", "tcga_gbm_clinical_v2.json")
config.setdefault("verify_partition", config["source"] == "tcga")

# CIBERSORTx mixture export — a manual round trip to an external web tool,
# deliberately outside `rule all`. See workflow/rules/export.smk.
config.setdefault("cibersortx_gene_set", "protein_coding")

SOURCE = str(config["source"]).strip().lower()
if SOURCE not in {"tcga", "synthetic"}:
    raise WorkflowError(
        f"config['source'] must be 'tcga' or 'synthetic', got {config['source']!r}"
    )

COHORT = str(config["cohort"])
WITH_EXPRESSION = _as_bool(config["with_expression"])
VERIFY_PARTITION = _as_bool(config["verify_partition"])

RAW = f"{config['raw_root']}/{COHORT}"
EXPR = f"{RAW}/expression"
INTERIM = str(config["interim_root"])
PROCESSED = str(config["processed_root"])
RESULTS = f"{config['results_root']}/{COHORT}"
PARTITIONS = str(config["partitions_root"])

# ── Layer artifacts ───────────────────────────────────────────────────────
# L0 — ingestion
CLINICAL_RAW = f"{RAW}/clinical_raw.parquet"
CLINICAL_COHORT = f"{RAW}/clinical_cohort.parquet"
SYNTHETIC_COHORT = f"{INTERIM}/{COHORT}_synthetic_cohort.parquet"

EXPR_FILE_MANIFEST = f"{EXPR}/file_manifest.parquet"
EXPR_GENE_ANNOTATION = f"{EXPR}/gene_annotation.parquet"
EXPR_COUNTS = f"{EXPR}/counts.parquet"
EXPR_TPM = f"{EXPR}/tpm.parquet"
EXPR_SAMPLE_METADATA = f"{EXPR}/sample_metadata.parquet"

# L0 -- MGMT/IDH call tables (scripts/download_tcga_mutations.py and
# scripts/download_tcga_methylation.py). Produced but not consumed by L2.
MUTATIONS = f"{RAW}/mutations"
IDH_FILE_MANIFEST = f"{MUTATIONS}/file_manifest.parquet"
IDH_CALLS = f"{MUTATIONS}/idh_calls.parquet"

METHYLATION = f"{RAW}/methylation"
METH_FILE_MANIFEST = f"{METHYLATION}/file_manifest.parquet"
MGMT_CALLS = f"{METHYLATION}/mgmt_calls.parquet"

# L1 — harmonization
HARMONIZED = f"{INTERIM}/{COHORT}_harmonized.parquet"

# L2 — profiling
PROFILED = f"{PROCESSED}/{COHORT}_profiled.parquet"

# L3 — feature store
FEATURE_STORE = f"{PROCESSED}/{COHORT}_feature_store.parquet"

# L6 — optimization  (`gliostrata optimize --out DIR` writes DIR/ranked_configurations.parquet)
CONFIG_DIR = f"{RESULTS}/configs"
RANKED_CONFIGS = f"{CONFIG_DIR}/ranked_configurations.parquet"

# L7 — handoff
PROTOCOL_DIR = f"{RESULTS}/protocols"

# External-tool export (manual round trip; not in `rule all`)
CIBERSORTX_DIR = f"{PROCESSED}/cibersortx"
CIBERSORTX_MIXTURE = (
    f"{CIBERSORTX_DIR}/tcga_gbm_mixture_{config['cibersortx_gene_set']}.txt"
)

# Guardrail / completion stamps
PARTITION_MANIFEST = f"{PARTITIONS}/{config['partition_manifest']}"
PARTITION_STAMP = f"{RESULTS}/partition_verified.txt"
PIPELINE_STAMP = f"{RESULTS}/pipeline_complete.txt"


def cohort_parquet(_wildcards=None) -> str:
    """The L0 cohort table the rest of the pipeline consumes, per `source`."""
    return CLINICAL_COHORT if SOURCE == "tcga" else SYNTHETIC_COHORT


def expression_targets() -> list:
    """L0 expression matrices, or nothing when expression ingest is disabled."""
    if not WITH_EXPRESSION:
        return []
    return [
        EXPR_FILE_MANIFEST,
        EXPR_GENE_ANNOTATION,
        EXPR_COUNTS,
        EXPR_TPM,
        EXPR_SAMPLE_METADATA,
    ]


def partition_targets() -> list:
    """Partition-drift stamp, or nothing when no frozen manifest applies."""
    if not VERIFY_PARTITION:
        return []
    if not PurePosixPath(PARTITION_MANIFEST).name:
        return []
    return [PARTITION_STAMP]
