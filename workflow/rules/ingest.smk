# ---------------------------------------------------------------------------
# workflow/rules/ingest.smk — L0 ingestion.
#
# Two mutually exclusive entry points, selected by config['source']:
#
#   source=tcga       REAL ingestion from the GDC REST API. Fully implemented.
#                     These rules TOUCH THE NETWORK and are the only rules in
#                     the workflow that do. They must never run in CI.
#   source=synthetic   `gliostrata synthetic-data` — offline, reproducible from a
#                     seed. This is NOT data; it is a schema-conformant fixture
#                     used to exercise the plumbing.
#
# NOTE on `gliostrata ingest --source tcga`: that CLI subcommand is a STUB and
# prints a "not yet implemented" message. Real TCGA ingestion lives in
# scripts/download_tcga_*.py, which is what these rules call. The old dvc.yaml
# called the stub, which is why it never produced its declared outputs.
# ---------------------------------------------------------------------------


# ── L0a: TCGA-GBM clinical ingest (REAL — network) ────────────────────────
# One GDC /cases query, mapped onto the GlioStrata cohort schema. No file
# inputs: the upstream dependency is the GDC API, which Snakemake cannot hash.
# Re-run deliberately with `--forcerun ingest_tcga_clinical`.
rule ingest_tcga_clinical:
    output:
        raw=CLINICAL_RAW,
        cohort=CLINICAL_COHORT,
    params:
        out_dir=RAW,
        project=config["gdc_project"],
    log:
        f"{RESULTS}/logs/L0_ingest_tcga_clinical.log",
    shell:
        "python scripts/download_tcga_clinical.py "
        "--out {params.out_dir} --project {params.project} > {log} 2>&1"


# ── L0b: TCGA-GBM bulk RNA-seq ingest (REAL — network, large) ─────────────
# Downloads STAR gene counts per sample into {EXPR}/_cache/ (~1.2 GB for the
# full cohort) and assembles gene x sample matrices. The cache is deliberately
# excluded from both git and DVC; see .dvcignore and docs/reproducibility.md.
#
# The matrices are produced but NOT YET CONSUMED by L2 — see profile.smk.
rule ingest_tcga_expression:
    input:
        cohort=CLINICAL_COHORT,
    output:
        file_manifest=EXPR_FILE_MANIFEST,
        gene_annotation=EXPR_GENE_ANNOTATION,
        counts=EXPR_COUNTS,
        tpm=EXPR_TPM,
        sample_metadata=EXPR_SAMPLE_METADATA,
    params:
        out_dir=EXPR,
        project=config["gdc_project"],
        limit=(
            f"--limit {config['expression_limit']}"
            if int(config["expression_limit"]) > 0
            else ""
        ),
    log:
        f"{RESULTS}/logs/L0_ingest_tcga_expression.log",
    shell:
        "python scripts/download_tcga_expression.py "
        "--out {params.out_dir} --project {params.project} "
        "--clinical {input.cohort} {params.limit} > {log} 2>&1"


# -- L0c: TCGA-GBM IDH mutation calling (REAL -- network) ------------------
# Mirrors the dvc.yaml stage of the same name. The script downloads and
# calls only; it joins nothing to the cohort (the cohort is read solely to
# list patients without a MAF as not-callable). The *_report.json beside the
# outputs and the _cache/ download scratch are not declared outputs.
#
# NOT PART OF `rule all` or `rule ingest`, on purpose: nothing downstream
# consumes the call tables yet, and they are a network download. Build
# explicitly:
#     snakemake --cores 1 --snakefile workflow/Snakefile ingest_tcga_mutations
rule ingest_tcga_mutations:
    input:
        code="scripts/download_tcga_mutations.py",
        cohort=CLINICAL_COHORT,
    output:
        file_manifest=IDH_FILE_MANIFEST,
        calls=IDH_CALLS,
    params:
        out_dir=MUTATIONS,
    log:
        f"{RESULTS}/logs/L0_ingest_tcga_mutations.log",
    shell:
        "python scripts/download_tcga_mutations.py "
        "--out {params.out_dir} > {log} 2>&1"


# -- L0d: TCGA-GBM MGMT methylation calling (REAL -- network, ~2.3 GB) -----
# Same pattern as L0c. The calls are a two-probe beta score, NOT MGMT-STP27;
# see docs/mgmt_idh_calling.md before using the binary call. The script also
# reads expression/tpm.parquet when it exists, only for the coverage section
# of its JSON report, so tpm is deliberately not an input here.
#
# NOT PART OF `rule all` or `rule ingest`, for the same reasons as L0c.
rule ingest_tcga_methylation:
    input:
        code="scripts/download_tcga_methylation.py",
        cohort=CLINICAL_COHORT,
    output:
        file_manifest=METH_FILE_MANIFEST,
        calls=MGMT_CALLS,
    params:
        out_dir=METHYLATION,
    log:
        f"{RESULTS}/logs/L0_ingest_tcga_methylation.log",
    shell:
        "python scripts/download_tcga_methylation.py "
        "--out {params.out_dir} > {log} 2>&1"


# ── synthetic fixture cohort (offline) ───────────────────────────────
# Deterministic given --config seed=... . Used by CI and by `make pipeline`.
rule ingest_synthetic:
    output:
        SYNTHETIC_COHORT,
    params:
        n=config["n_synthetic"],
        seed=config["seed"],
    log:
        f"{RESULTS}/logs/L0_ingest_synthetic.log",
    shell:
        "gliostrata synthetic-data --n {params.n} --seed {params.seed} "
        "--out {output} > {log} 2>&1"
