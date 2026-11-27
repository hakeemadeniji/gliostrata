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


# The real L1 clinical step (not a stub). `gliostrata harmonize-clinical` emits
# the clinical table with explicitly owned columns -- no L2 or molecular names;
# a column under such a name that carries a value stops it -- and a separate
# label table. Every missing value is traced to the raw GDC record and gets a
# reason code. Nothing is imputed. Used by build_feature_store when source=tcga.
rule harmonize_clinical:
    input:
        cohort=CLINICAL_COHORT,
        raw=CLINICAL_RAW,
    output:
        clinical=CLINICAL_HARMONIZED,
        labels=CLINICAL_LABELS,
        provenance=CLINICAL_PROVENANCE,
    log:
        f"{RESULTS}/logs/L1_harmonize_clinical.log",
    shell:
        "gliostrata harmonize-clinical --cohort {input.cohort} --raw {input.raw} "
        "--clinical-out {output.clinical} --labels-out {output.labels} > {log} 2>&1"
