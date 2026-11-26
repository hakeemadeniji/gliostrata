# ---------------------------------------------------------------------------
# workflow/rules/partition.smk — frozen train/test partition verification.
#
# This is the one rule in the workflow that is a scientific guardrail rather
# than a build step. The split is a pure function of (seed, patient_id), so a
# sealed manifest must keep reproducing as the cohort is re-ingested. If it
# stops reproducing, a patient has moved between train and test and every
# downstream evaluation is contaminated.
#
# The rule therefore FAILS the DAG rather than regenerating anything:
# `gliostrata freeze-partition --verify` writes nothing and exits non-zero on
# drift. There is deliberately NO rule in this workflow that creates or
# overwrites a manifest — freezing a partition is a human, pre-registered act.
# ---------------------------------------------------------------------------


rule verify_partition:
    input:
        cohort=cohort_parquet,
        manifest=PARTITION_MANIFEST,
    output:
        touch(PARTITION_STAMP),
    log:
        f"{RESULTS}/logs/partition_verify.log",
    shell:
        "gliostrata freeze-partition --input {input.cohort} "
        "--manifest {input.manifest} --verify > {log} 2>&1"
