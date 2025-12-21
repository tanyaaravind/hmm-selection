#!/usr/bin/env python3

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import stats

OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
AF_FILE = OUT / "abo_expanded" / "allele_frequencies.csv"
DELTA_FILE = OUT / "abo_expanded" / "delta_af_nonzero.csv"

sns.set(style="whitegrid")


def compute_and_plot_fst(pop_a: str, pop_b: str):
    print(f"Computing FST for {pop_a} vs {pop_b}...")
    af = pd.read_csv(AF_FILE)
    fst_df = stats.fst_from_af_table(af, pop_a, pop_b)
    out_csv = OUT / f"benchmarks_fst_{pop_a.lower()}_{pop_b.lower()}.csv"
    fst_df.to_csv(out_csv, index=False)

    plt.figure(figsize=(6, 4))
    sns.histplot(fst_df["fst"].clip(lower=0), bins=60, kde=False)
    plt.title(f"FST distribution: {pop_a} vs {pop_b}")
    plt.xlabel("FST")
    plt.ylabel("Count")
    plt.tight_layout()
    hist_png = OUT / f"fst_hist_{pop_a.lower()}_{pop_b.lower()}.png"
    plt.savefig(hist_png)
    plt.close()

    plt.figure(figsize=(8, 3))
    plt.plot(fst_df["pos"], fst_df["fst"], marker=".", linestyle="", markersize=3)
    top10 = fst_df.sort_values("fst", ascending=False).head(10)
    plt.scatter(top10["pos"], top10["fst"], color="red", s=30, label="Top 10")
    plt.title(f"FST per-site: {pop_a} vs {pop_b}")
    plt.xlabel("Position (bp)")
    plt.ylabel("FST")
    plt.legend()
    plt.tight_layout()
    pos_png = OUT / f"fst_vs_pos_{pop_a.lower()}_{pop_b.lower()}.png"
    plt.savefig(pos_png)
    plt.close()

    print(f"Saved {out_csv}, {hist_png}, {pos_png}")
    return fst_df


def tajimas_sliding_and_annotate(window_size=50, step=10, top_delta_pct=95):
    print("Computing sliding Tajima's D for YRI with annotations...")
    af = pd.read_csv(AF_FILE)
    delta = pd.read_csv(DELTA_FILE)

    yri = af[af["pop"] == "YRI"].sort_values("pos").reset_index(drop=True)
    positions = yri["pos"].values
    alt_counts = yri["alt_count"].values
    total = yri["total_alleles"].values

    thresh = np.percentile(delta["delta_af"].values, top_delta_pct)
    high_delta_pos = set(delta[delta["delta_af"] >= thresh]["pos"].values)

    rows = []
    for start in range(0, len(positions) - window_size + 1, step):
        win_pos = positions[start : start + window_size]
        win_alt = alt_counts[start : start + window_size]
        win_total = total[start : start + window_size]
        try:
            d = stats.tajimas_d_from_counts(win_alt, win_total)
        except Exception:
            d = float("nan")
        overlap = any((p in high_delta_pos) for p in win_pos)
        rows.append({
            "start_pos": int(win_pos[0]),
            "end_pos": int(win_pos[-1]),
            "mid_pos": int((win_pos[0] + win_pos[-1]) / 2),
            "n_sites": len(win_pos),
            "tajimas_d": d,
            "overlap_high_delta": bool(overlap),
        })

    td_df = pd.DataFrame(rows)
    td_out = OUT / f"benchmarks_tajimas_yri_sliding_w{window_size}_s{step}.csv"
    td_df.to_csv(td_out, index=False)

    plt.figure(figsize=(10, 3))
    plt.plot(td_df["mid_pos"], td_df["tajimas_d"], marker="o", linestyle="-")
    for _, row in td_df[td_df["overlap_high_delta"]].iterrows():
        plt.axvspan(row["start_pos"], row["end_pos"], color="red", alpha=0.12)
    plt.title(f"Tajima's D (YRI) — window={window_size}, step={step}")
    plt.xlabel("Position (bp)")
    plt.ylabel("Tajima's D")
    plt.tight_layout()
    png = OUT / f"tajimas_yri_sliding_w{window_size}_s{step}.png"
    plt.savefig(png)
    plt.close()

    print(f"Saved Tajima's D table to {td_out} and plot to {png}")
    return td_df, thresh


def write_analysis(fst_results, td_df, thresh, top_k=10):
    print("Writing short analysis summary...")
    lines = []
    lines.append("# Benchmark summary\n")
    lines.append("## FST results (summary)\n")
    for key, df in fst_results.items():
        stats_series = df["fst"].describe()
        lines.append(f"**{key}** — count={int(stats_series['count'])}, mean={stats_series['mean']:.3f}, median={df['fst'].median():.3f}, max={stats_series['max']:.3f}\n")
        top = df.sort_values("fst", ascending=False).head(top_k)
        lines.append("Top sites (pos, fst):\n")
        for _, r in top.iterrows():
            lines.append(f"- {int(r['pos'])}: {r['fst']:.4f}\n")
        lines.append("\n")

    lines.append("## Tajima's D (YRI) — sliding windows\n")
    lines.append(f"Window size/step: see file; top deltaAF threshold (>= {100.0 - (100 - (100 * 0 + 0)):0.0}) used to mark windows: {thresh:.6f}\n")
    lines.append(td_df["tajimas_d"].describe().to_string())
    lines.append("\n\n")

    lines.append("## Notes for Sheki Okwayo (lso24) \n")
    lines.append("- These benchmarks used the allele-frequency table in `results/abo_expanded/allele_frequencies.csv` and the filtered delta AF table in `results/abo_expanded/delta_af_nonzero.csv`.\n")
    lines.append("- I flagged windows overlapping positions with DeltaAF >= the 95th percentile (shown as red-shaded regions in the Tajima's D plot).\n")
    lines.append("- If you want different population pairs, window sizes, or more annotation (e.g., overlapping genes), tell me and I’ll extend this.\n")

    out = OUT / "benchmark_analysis.md"
    out.write_text("\n".join(lines))
    print(f"Wrote analysis to {out}")


if __name__ == "__main__":
    pairs = [("YRI", "CEU"), ("YRI", "CHB"), ("CEU", "CHB")]
    fst_results = {}
    for a, b in pairs:
        df = compute_and_plot_fst(a, b)
        fst_results[f"{a}_vs_{b}"] = df

    td_df, thresh = tajimas_sliding_and_annotate(window_size=50, step=10, top_delta_pct=95)
    write_analysis(fst_results, td_df, thresh)
    print("Done.")
