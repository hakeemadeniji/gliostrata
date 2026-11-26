"""Run the full GlioStrata L0→L7 pipeline on synthetic data.

Standalone script — equivalent to `gliostrata pipeline-dry-run` but writes
intermediate outputs to data/ so you can inspect them.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from gliostrata.data.harmonize import harmonize
from gliostrata.data.ingest import generate_synthetic_cohort
from gliostrata.features.store import build_feature_store
from gliostrata.handoff.protocol import generate_protocols
from gliostrata.optimization.bayesian import optimize_configurations
from gliostrata.prediction.survival import predict_survival
from gliostrata.profiling.states import deconvolve_states
from gliostrata.profiling.subtype import classify_subtypes
from gliostrata.profiling.tme import deconvolve_tme
from gliostrata.utils.io import ensure_dir, write_json, write_parquet


def main() -> None:
    parser = argparse.ArgumentParser(description="Run full GlioStrata pipeline on synthetic data.")
    parser.add_argument("--n", type=int, default=30, help="Number of synthetic patients.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-dir", type=str, default="data/results",
                        help="Output directory for pipeline artifacts.")
    args = parser.parse_args()

    out_dir = ensure_dir(args.out_dir)

    print(f"GlioStrata full pipeline run — n={args.n} synthetic patients\n")

    # L0
    cohort = generate_synthetic_cohort(n=args.n, seed=args.seed)
    write_parquet(cohort, out_dir / "L0_cohort.parquet")
    print(f"  L0  ingestion       : {len(cohort)} patients → {out_dir}/L0_cohort.parquet")

    # L1
    harmonized = harmonize(cohort)
    write_parquet(harmonized, out_dir / "L1_harmonized.parquet")
    print(f"  L1  harmonization   : {len(harmonized)} samples → {out_dir}/L1_harmonized.parquet")

    # L2
    subtypes = classify_subtypes(harmonized, seed=args.seed)
    states = deconvolve_states(harmonized, seed=args.seed)
    tme = deconvolve_tme(harmonized, seed=args.seed)
    write_parquet(subtypes, out_dir / "L2_subtypes.parquet")
    write_parquet(states, out_dir / "L2_states.parquet")
    write_parquet(tme, out_dir / "L2_tme.parquet")
    print(f"  L2  profiling       : subtypes/states/TME → {out_dir}/L2_*.parquet")

    # L3
    feature_store = build_feature_store(harmonized, subtypes, states, tme)
    write_parquet(feature_store, out_dir / "L3_feature_store.parquet")
    print(f"  L3  feature store   : {feature_store.shape} → {out_dir}/L3_feature_store.parquet")

    # L5
    survival_pred = predict_survival(feature_store, seed=args.seed)
    write_parquet(survival_pred, out_dir / "L5_survival.parquet")
    print(f"  L5  survival pred   : median predicted OS = "
          f"{survival_pred['predicted_survival_months'].median():.1f} months "
          f"→ {out_dir}/L5_survival.parquet")

    # L6
    configs = optimize_configurations(feature_store, survival_pred,
                                       n_top=5, n_candidates=50, seed=args.seed)
    write_parquet(configs, out_dir / "L6_ranked_configs.parquet")
    print(f"  L6  optimizer       : {len(configs)} configs "
          f"({configs['patient_id'].nunique()} patients × {configs['rank'].nunique()} ranks) "
          f"→ {out_dir}/L6_ranked_configs.parquet")

    # Personalization check
    from gliostrata.optimization.bayesian import evaluate_personalization
    pers = evaluate_personalization(configs)
    print(f"      personalization : {pers['n_distinct_top5_cell_payload_combos']} "
          f"distinct top-5 cell×payload combos "
          f"(threshold ≥ {pers['threshold_min_distinct']} → {'PASS' if pers['passes_threshold'] else 'FAIL'})")

    # L7
    protocols = generate_protocols(configs, feature_store, survival_pred, top_n=1)
    proto_dir = ensure_dir(out_dir / "L7_protocols")
    for pid, proto in protocols.items():
        write_json(proto, proto_dir / f"protocol_{pid}.json")
    print(f"  L7  handoff         : {len(protocols)} protocol docs → {proto_dir}/")

    # End-to-end success metric (T5 stub on synthetic data)
    gains = configs[configs["rank"] == 1]["survival_gain_vs_stupp_months"]
    success_fraction = (gains >= 3.0).mean()
    harm_fraction = (gains < 0).mean()
    print(f"\n  T5 stub metric      : {success_fraction:.0%} of patients with ≥ 3-month "
          f"gain (threshold ≥ 50%); {harm_fraction:.0%} predicted harm (threshold ≤ 10%)")
    print(f"\nPipeline complete. All artifacts in {out_dir}/")


if __name__ == "__main__":
    main()
