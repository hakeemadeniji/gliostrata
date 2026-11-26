# ---------------------------------------------------------------------------
# workflow/rules/features.smk — L3 unified patient feature store.
#
# `gliostrata build-feature-store --input X --out Y` takes the HARMONIZED
# cohort and re-runs the three L2 profiling functions internally before
# assembling the store. It does NOT read the L2 `profiled.parquet` table.
#
# So the only honest input here is HARMONIZED. The previous Snakefile declared
# L2_subtypes/L2_states/L2_tme as inputs to this rule and then passed only the
# cohort on the command line — the DAG asserted a dependency the command never
# used, and the extra inputs merely forced needless re-runs.
#
# The duplicated L2 computation is real technical debt (two code paths that can
# drift apart). It belongs in src/, not in the workflow, so it is recorded here
# rather than papered over with a fake dependency edge.
# ---------------------------------------------------------------------------


rule build_feature_store:
    input:
        harmonized=HARMONIZED,
    output:
        FEATURE_STORE,
    log:
        f"{RESULTS}/logs/L3_feature_store.log",
    shell:
        "gliostrata build-feature-store "
        "--input {input.harmonized} --out {output} > {log} 2>&1"
