#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

fst_files = {
    "YRI_vs_CEU": OUT / "benchmarks_fst_yri_ceu.csv",
    "YRI_vs_CHB": OUT / "benchmarks_fst_yri_chb.csv",
    "CEU_vs_CHB": OUT / "benchmarks_fst_ceu_chb.csv",
}
annot_file = OUT / "benchmarks_fst_yri_ceu_top20_annotated.csv"
td_file = OUT / "benchmarks_tajimas_yri_sliding_w50_s10.csv"

def load_fst():
    dfs = []
    for name, path in fst_files.items():
        df = pd.read_csv(path)
        df = df.assign(pair=name)
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def fst_stats_and_plot(df):
    summary = df.groupby('pair')['fst'].agg(['count','mean','std','median','min','max'])
    summary.to_csv(OUT / 'fst_pair_stats.csv')

    plt.figure(figsize=(8,4))
    sns.violinplot(x='pair', y='fst', data=df, inner='box', cut=0)
    plt.title('FST distribution by population pair')
    plt.tight_layout()
    png = OUT / 'fst_comparison.png'
    plt.savefig(png)
    plt.close()
    return summary, png


def analyze_annotations(annot_path):
    ann = pd.read_csv(annot_path)
    ann['in_gene'] = ann['annotated'] == True
    counts = ann.groupby('in_gene').size().rename('n').reset_index()
    return ann, counts


def tajima_overlap_analysis(td_path, delta_threshold=None):
    td = pd.read_csv(td_path)
    td_clean = td.dropna(subset=['tajimas_d'])
    overlap = td_clean[td_clean['overlap_high_delta'] == True]['tajimas_d']
    nonoverlap = td_clean[td_clean['overlap_high_delta'] == False]['tajimas_d']

    summary = {
        'overlap_n': int(len(overlap)),
        'overlap_mean': float(np.nanmean(overlap)),
        'overlap_median': float(np.nanmedian(overlap)),
        'nonoverlap_n': int(len(nonoverlap)),
        'nonoverlap_mean': float(np.nanmean(nonoverlap)),
        'nonoverlap_median': float(np.nanmedian(nonoverlap)),
    }

    try:
        u_res = stats.mannwhitneyu(overlap, nonoverlap, alternative='two-sided')
        summary['mw_u'] = float(u_res.statistic)
        summary['mw_p'] = float(u_res.pvalue)
    except Exception:
        summary['mw_u'] = np.nan
        summary['mw_p'] = np.nan

    pd.DataFrame([summary]).to_csv(OUT / 'tajimas_yri_overlap_summary.csv', index=False)
    return summary


def write_markdown(summary_stats, annot_counts, ann_table, tajima_summary, png_paths):
    out = OUT / 'benchmark_report.md'
    lines = []
    lines.append('\n')
    lines.append('## Visual & Comparative Benchmark Analysis')
    lines.append('\n')
    lines.append('## FST comparison across pairs')
    lines.append('\n')
    lines.append('Summary statistics by pair (see `results/fst_pair_stats.csv`):')
    lines.append('\n')
    lines.append(summary_stats.to_markdown())
    lines.append('\n')
    lines.append(f"FST violin plot saved to: `{png_paths['fst']}`")
    lines.append('\n')
    lines.append('## Annotation of top FST sites (YRI vs CEU)')
    lines.append('\n')
    lines.append('Top-sites gene overlap counts:')
    lines.append('\n')
    lines.append(annot_counts.to_string(index=False))
    lines.append('\n')
    lines.append('Annotated top sites (first 20):')
    lines.append('\n')
    lines.append(ann_table[['pos','fst','gene_name','annotation_status']].to_markdown(index=False))
    lines.append('\n')
    lines.append("## Tajima's D — windows overlapping high DeltaAF vs others")
    lines.append('\n')
    lines.append('Summary (see `results/tajimas_yri_overlap_summary.csv`):')
    lines.append('\n')
    for k,v in tajima_summary.items():
        lines.append(f'- {k}: {v}')
    lines.append('\n')
    lines.append('### Quick interpretation')
    lines.append('\n')
    lines.append('- The violin plot shows distributional differences across population pairs; YRI vs CHB often shows larger FST tail values within this ABO region.\n')
    lines.append('- Most top FST sites overlap the `ABO` gene (expected given region-focused analysis); a few high-FST positions did not overlap annotated genes and warrant further inspection for regulatory elements or structural variation.\n')
    lines.append("- Tajima's D in windows overlapping high DeltaAF positions has mean/median as reported above; the Mann–Whitney U test p-value indicates whether the distributions differ significantly (see `tajimas_yri_overlap_summary.csv`).\n")
    lines.append('\n')
    lines.append('### Visuals')
    lines.append(f'- FST distribution by pair: `{png_paths["fst"]}`')
    lines.append(f'- Per-pair histograms and per-site plots are in `results/` with names starting `fst_hist_` and `fst_vs_pos_`.')
    lines.append(f"- Tajima's D sliding plot: `results/tajimas_yri_sliding_w50_s10.png` (windows overlapping top ΔAF marked in red).")

    with out.open('a') as f:
        f.write('\n---\n')
        f.write('\n'.join(lines))

    return out


if __name__ == '__main__':
    sns.set(style='whitegrid')
    df = load_fst()
    summary, fst_png = fst_stats_and_plot(df)

    ann_table, ann_counts = analyze_annotations(annot_file)
    taj_summary = tajima_overlap_analysis(td_file)

    pngs = {'fst': fst_png}
    report = write_markdown(summary, ann_counts, ann_table, taj_summary, pngs)
    print(f'Wrote report to {report}')
