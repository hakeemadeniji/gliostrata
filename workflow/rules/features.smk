# ---------------------------------------------------------------------------
# workflow/rules/features.smk — L2 + L3: profiling and the unified feature store.
#
# `gliostrata build-feature-store` runs in one of two modes, chosen by
# config `source` (see feature_store_inputs / feature_store_args in common.smk):
#
#   source=synthetic  `--mode synthetic`: the harmonized synthetic cohort, the
#                     three L2 STUBS, and the store. Offline; what CI and the
#                     fixture run use. Every file it writes says "synthetic".
#   source=tcga       `--mode real`: the TCGA expression input (tpm.parquet and
#                     its gene annotation) through the real L2 -- single-sample
#                     nearest-centroid subtype calls against the frozen centroid
#                     set, and the state and TME enrichment scores -- and the L3
#                     store. A stub output consumed here, a missing centroid set
#                     or a path containing "cgga" stops the command. The SHA-256
#                     of every input, the expression matrix included, is written
#                     to the build report.
#
# The store keeps every cohort patient (one row each) and never filters. It is
# four files: features, labels (outcomes are labels, never features), the
# missing-value reasons, and the column metadata. Model-ready populations are
# named views whose flow counts are in the build report.
# ---------------------------------------------------------------------------


rule build_feature_store:
    input:
        unpack(feature_store_inputs),
    output:
        store=FEATURE_STORE,
        labels=FEATURE_LABELS,
        missing=FEATURE_MISSING,
        columns=FEATURE_COLUMNS_JSON,
        report=FEATURE_REPORT,
    params:
        args=feature_store_args,
    log:
        f"{RESULTS}/logs/L3_feature_store.log",
    shell:
        "gliostrata build-feature-store {params.args} "
        "--out {output.store} --labels-out {output.labels} "
        "--missing-out {output.missing} --columns-out {output.columns} "
        "--report {output.report} > {log} 2>&1"
