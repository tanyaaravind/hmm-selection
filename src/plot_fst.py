"""
Quick FST plot for the example ABO allele-frequency table.

Outputs a small line/point plot of per-site FST for each population pair.

Usage:
  python src/plot_fst.py --af data/example_abo_af.csv --out results/fst_abo_example.png
"""

import argparse
import os
from typing import List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from stats import fst_from_af_table


def compute_pairs(af_df: pd.DataFrame, pairs: List[Tuple[str, str]]) -> pd.DataFrame:
    """Compute FST for each population pair and return a tidy dataframe."""
    frames = []
    for a, b in pairs:
        fst_df = fst_from_af_table(af_df, a, b)
        fst_df["pair"] = f"{a}-{b}"
        frames.append(fst_df)
    return pd.concat(frames, ignore_index=True)


def plot_fst(fst_df: pd.DataFrame, out_path: str):
    """Plot per-site FST for each pair and save to out_path."""
    sns.set_style("whitegrid")
    plt.figure(figsize=(10, 4))

    # Convert positions to kb relative to start for readability
    start = fst_df["pos"].min()
    fst_df = fst_df.copy()
    fst_df["pos_kb"] = (fst_df["pos"] - start) / 1000

    palette = {"YRI-CEU": "#1f77b4", "YRI-CHB": "#2ca02c", "CEU-CHB": "#d62728"}
    for pair, group in fst_df.groupby("pair"):
        plt.plot(group["pos_kb"], group["fst"], marker="o", label=pair, color=palette.get(pair))

    plt.xlabel(f"Position (kb) relative to {start:,}")
    plt.ylabel("FST")
    plt.title("Per-site FST across ABO example SNPs")
    plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved FST plot to: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Plot per-site FST for ABO example table.")
    parser.add_argument("--af", default="data/example_abo_af.csv", help="Path to allele frequency CSV.")
    parser.add_argument(
        "--out", default="results/fst_abo_example.png", help="Output PNG path (directories auto-created)."
    )
    args = parser.parse_args()

    af_df = pd.read_csv(args.af)
    pairs = [("YRI", "CEU"), ("YRI", "CHB"), ("CEU", "CHB")]
    fst_df = compute_pairs(af_df, pairs)
    plot_fst(fst_df, args.out)


if __name__ == "__main__":
    main()
