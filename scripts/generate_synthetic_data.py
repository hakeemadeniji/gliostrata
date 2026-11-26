"""Generate a synthetic GBM cohort and save it to Parquet.

Standalone script — equivalent to `gliostrata synthetic-data`.
Useful when you want to generate data without invoking the CLI.
"""
import argparse
from pathlib import Path

from gliostrata.data.ingest import generate_synthetic_cohort
from gliostrata.utils.io import write_parquet


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a synthetic GBM cohort.")
    parser.add_argument("--n", type=int, default=50, help="Number of patients (default: 50).")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42).")
    parser.add_argument("--out", type=str, default="data/interim/synthetic_cohort.parquet",
                        help="Output Parquet path.")
    args = parser.parse_args()

    df = generate_synthetic_cohort(n=args.n, seed=args.seed)
    written = write_parquet(df, args.out)
    print(f"Generated {len(df)} synthetic patients → {written}")
    print(f"Columns: {list(df.columns)}")
    print(f"\nSubtype distribution:")
    print(df["subtype"].value_counts())
    print(f"\nOverall survival summary (months):")
    print(df["overall_survival_months"].describe())


if __name__ == "__main__":
    main()
