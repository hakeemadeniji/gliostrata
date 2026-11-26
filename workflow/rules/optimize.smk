# ---------------------------------------------------------------------------
# workflow/rules/optimize.smk — L5 survival prediction + L6 configuration search.
#
# There is no separate L5 rule: `gliostrata optimize` calls
# gliostrata.prediction.survival.predict_survival internally and does not write
# the predictions out. L5 therefore has no materialised artifact today. The old
# Snakefile invented a `data/results/L5_survival.parquet` rule driven by an
# inline `python -c` snippet — that artifact is not part of any real interface
# and is not reproduced here.
#
# STUB LAYERS:
#   L5 (gliostrata.prediction.survival) — not a fitted survival model over the
#       real cohort; see src/ and docs/reproducibility.md.
#   L6 (gliostrata.optimization.bayesian) — a ranked enumeration, not a
#       Bayesian-optimization loop against a validated response surface.
#
# Note the --out contract: `gliostrata optimize --out DIR` writes
# DIR/ranked_configurations.parquet, so this rule's output is the file inside
# the directory, with the directory passed as a param.
# ---------------------------------------------------------------------------


rule optimize:
    input:
        features=FEATURE_STORE,
        # Partition verification is a hard gate: with config verify_partition=true
        # the optimizer cannot run until the frozen split still reproduces.
        partition=partition_targets(),
    output:
        RANKED_CONFIGS,
    params:
        out_dir=CONFIG_DIR,
        n_top=config["n_top"],
        seed=config["seed"],
    log:
        f"{RESULTS}/logs/L6_optimize.log",
    shell:
        "gliostrata optimize --input {input.features} --out {params.out_dir} "
        "--n-top {params.n_top} --seed {params.seed} > {log} 2>&1"
