# ---------------------------------------------------------------------------
# workflow/rules/profile.smk — L2 tumor profiling.
#
# STUB LAYER — and the most consequential stub in the pipeline.
#
# `gliostrata profile` calls three functions:
#   gliostrata.profiling.subtype.classify_subtypes  — Verhaak subtype
#   gliostrata.profiling.states.deconvolve_states   — Neftel cell states
#   gliostrata.profiling.tme.deconvolve_tme         — TME composition
#
# None of them read an expression matrix. They derive their outputs from the
# clinical/synthetic cohort table alone, so the subtype calls and state/TME
# fractions produced here are NOT measurements. The real implementations need
# the L0 expression matrices (counts.parquet / tpm.parquet), which the ingest
# rules already produce.
#
# That gap is why the expression matrices are declared as a top-level target in
# the Snakefile rather than as an input to this rule: declaring them as an input
# would make the DAG claim a data dependency the code does not have.
#
# The REAL L2 is not this rule: with source=tcga, build_feature_store
# (`gliostrata build-feature-store --mode real`, features.smk) reads the
# expression input and runs the single-sample subtype classifier and the
# state / TME scoring. This stub rule's output is consumed by nothing.
# ---------------------------------------------------------------------------


rule profile:
    input:
        harmonized=HARMONIZED,
        # Intentionally NOT an input: EXPR_TPM / EXPR_COUNTS. See header.
    output:
        PROFILED,
    log:
        f"{RESULTS}/logs/L2_profile.log",
    shell:
        "gliostrata profile --input {input.harmonized} --out {output} > {log} 2>&1"
