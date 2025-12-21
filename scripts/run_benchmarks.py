#!/usr/bin/env python3

from pathlib import Path
import sys
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import stats
AF_FILE = ROOT / "results" / "abo_expanded" / "allele_frequencies.csv"
OUT_DIR = ROOT / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

def main():
    print("Loading allele-frequency table:", AF_FILE)
    af = pd.read_csv(AF_FILE)

    print("Computing per-site FST (YRI vs CEU)...")
    fst_df = stats.fst_from_af_table(af, "YRI", "CEU")
    fst_out = OUT_DIR / "benchmarks_fst_yri_ceu.csv"
    fst_df.to_csv(fst_out, index=False)
    print(f"Saved FST table to {fst_out}")

    print("FST summary (YRI vs CEU):")
    print(fst_df["fst"].describe())
    top = fst_df.sort_values("fst", ascending=False).head(10)
    print("Top 10 FST sites:")
    print(top.to_string(index=False))

    print("Computing Tajima's D for YRI in non-overlapping windows (50 sites)...")
    yri = af[af["pop"] == "YRI"].sort_values(["chrom", "pos"]).reset_index(drop=True)
    positions = yri["pos"].values
    alt_counts = yri["alt_count"].values
    total = yri["total_alleles"].values

    window_size = 50
    rows = []
    for i in range(0, len(positions), window_size):
        win_pos = positions[i : i + window_size]
        win_alt = alt_counts[i : i + window_size]
        win_total = total[i : i + window_size]
        if len(win_pos) == 0:
            continue
        try:
            d = stats.tajimas_d_from_counts(win_alt, win_total)
        except Exception as e:
            d = float("nan")
        rows.append({
            "chrom": "9",
            "start_pos": int(win_pos[0]),
            "end_pos": int(win_pos[-1]),
            "n_sites": len(win_pos),
            "tajimas_d": d,
        })

    td_df = pd.DataFrame(rows)
    td_out = OUT_DIR / "benchmarks_tajimas_yri.csv"
    td_df.to_csv(td_out, index=False)
    print(f"Saved Tajima's D windows to {td_out}")
    print("Tajima's D summary:")
    print(td_df["tajimas_d"].describe())

if __name__ == "__main__":
    main()
