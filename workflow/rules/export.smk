# ---------------------------------------------------------------------------
# workflow/rules/export.smk — exports to external tools.
#
# NOT PART OF `rule all`, on purpose.
#
# CIBERSORTx (https://cibersortx.stanford.edu) is a web tool. It cannot be
# driven from this pipeline, so L2 TME deconvolution is a MANUAL round trip:
# this rule produces the upload file, a human uploads it, and the returned
# fractions are re-ingested separately. There is therefore a hole in the DAG
# that no amount of Snakemake can close, and pretending otherwise by chaining
# a downstream rule onto this output would be a lie.
#
# Build it explicitly:
#     snakemake --cores 1 --snakefile workflow/Snakefile cibersortx_export
#
# The output is ~39 MB of tab-delimited text, regenerable from tpm.parquet in
# seconds. It lives under data/processed/, which is gitignored, and is a poor
# DVC candidate for the same reason (see docs/reproducibility.md).
# ---------------------------------------------------------------------------


rule cibersortx_export:
    """Export the bulk TPM matrix as a CIBERSORTx mixture file (manual upload)."""
    input:
        tpm=EXPR_TPM,
        annotation=EXPR_GENE_ANNOTATION,
        cohort=CLINICAL_COHORT,
    output:
        CIBERSORTX_MIXTURE,
    params:
        out_dir=CIBERSORTX_DIR,
        gene_set=config["cibersortx_gene_set"],
    log:
        f"{RESULTS}/logs/export_cibersortx.log",
    shell:
        "python scripts/export_cibersortx.py --tpm {input.tpm} "
        "--annotation {input.annotation} --cohort {input.cohort} "
        "--out {params.out_dir} --gene-set {params.gene_set} > {log} 2>&1"
