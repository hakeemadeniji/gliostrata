# ---------------------------------------------------------------------------
# workflow/rules/harmonize.smk — L1 harmonization + QC.
#
# STUB LAYER. `gliostrata harmonize` calls gliostrata.data.harmonize.harmonize,
# which currently normalises the cohort schema and applies basic QC filters.
# The architecture's L1 scope — cross-cohort batch correction (ComBat/limma),
# platform harmonization, imputation of missing molecular covariates — is NOT
# implemented. Treat the output as a schema-conformant passthrough, not as a
# harmonized multi-cohort table.
# ---------------------------------------------------------------------------


rule harmonize:
    input:
        cohort=cohort_parquet,
    output:
        HARMONIZED,
    log:
        f"{RESULTS}/logs/L1_harmonize.log",
    shell:
        "gliostrata harmonize --input {input.cohort} --out {output} > {log} 2>&1"
