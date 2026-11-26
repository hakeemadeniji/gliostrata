# ---------------------------------------------------------------------------
# workflow/rules/handoff.smk — L7 wet-lab protocol generation.
#
# STUB LAYER. `gliostrata generate-protocol` emits one JSON document per
# patient from the ranked configurations and the feature store. The documents
# are structurally complete but their content is templated from stubbed L2/L5/L6
# outputs, so no protocol produced today is executable at a bench.
#
# The output is a directory of protocol_<patient_id>.json files whose number is
# only known after L6 runs, so it is declared with directory().
# ---------------------------------------------------------------------------


rule generate_protocol:
    input:
        configs=RANKED_CONFIGS,
        features=FEATURE_STORE,
    output:
        directory(PROTOCOL_DIR),
    log:
        f"{RESULTS}/logs/L7_generate_protocol.log",
    shell:
        "gliostrata generate-protocol --input {input.configs} "
        "--feature-store {input.features} --out {output} > {log} 2>&1"
