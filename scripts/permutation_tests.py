#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
OUT.mkdir(exist_ok=True)

FST_FILE = OUT / 'benchmarks_fst_yri_ceu.csv'
TOP_ANN = OUT / 'benchmarks_fst_yri_ceu_top20_annotated.csv'
TD_FILE = OUT / 'benchmarks_tajimas_yri_sliding_w50_s10.csv'

np.random.seed(12345)

def max_fst_perm(k=20, nperm=5000):
    fst = pd.read_csv(FST_FILE)
    observed_top = fst.nlargest(k, 'fst')
    observed_max = observed_top['fst'].max()

    maxs = np.empty(nperm)
    positions = fst['pos'].values
    fsts = fst['fst'].values
    start = time.time()
    for i in range(nperm):
        idx = np.random.choice(len(fst), size=k, replace=False)
        maxs[i] = fsts[idx].max()
    elapsed = time.time() - start

    pval = (np.sum(maxs >= observed_max) + 1) / (nperm + 1)

    df = pd.DataFrame({'perm_max_fst': maxs})
    df.to_csv(OUT / 'perm_max_fst_top20.csv', index=False)

    plt.figure(figsize=(6,4))
    sns.histplot(maxs, bins=50, kde=False)
    plt.axvline(observed_max, color='red', linestyle='--', label=f'observed max={observed_max:.4f}')
    plt.title(f'Permutation: max FST of {k} random sites (n={nperm})')
    plt.xlabel('max FST')
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / 'perm_max_fst_top20.png')
    plt.close()

    return observed_max, pval, elapsed


def gene_overlap_perm(k=20, nperm=5000):
    top_ann = pd.read_csv(TOP_ANN)
    abo_rows = top_ann[top_ann['gene_name'] == 'ABO']
    if not abo_rows.empty:
        gene_start = int(abo_rows['start'].min())
        gene_end = int(abo_rows['end'].max())
    else:
        ann_all = top_ann.dropna(subset=['start','end'])
        if not ann_all.empty:
            gene_start = int(ann_all['start'].min())
            gene_end = int(ann_all['end'].max())
        else:
            gene_start = None
            gene_end = None

    fst = pd.read_csv(FST_FILE)
    positions = fst['pos'].values

    observed_positions = top_ann['pos'].astype(int).values
    if gene_start is not None:
        observed_count = np.sum((observed_positions >= gene_start) & (observed_positions <= gene_end))
    else:
        observed_count = 0

    counts = np.empty(nperm, dtype=int)
    for i in range(nperm):
        sample_pos = np.random.choice(positions, size=k, replace=False)
        if gene_start is not None:
            counts[i] = np.sum((sample_pos >= gene_start) & (sample_pos <= gene_end))
        else:
            counts[i] = 0

    pval = (np.sum(counts >= observed_count) + 1) / (nperm + 1)

    pd.DataFrame({'perm_counts': counts}).to_csv(OUT / 'perm_gene_overlap_top20.csv', index=False)

    plt.figure(figsize=(6,4))
    sns.histplot(counts, bins=range(counts.min(), counts.max()+2), discrete=True)
    plt.axvline(observed_count, color='red', linestyle='--', label=f'observed count={observed_count}')
    plt.title(f'Permutation: count of {k} random sites overlapping gene region (n={nperm})')
    plt.xlabel('count overlapping gene')
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / 'perm_gene_overlap_top20.png')
    plt.close()

    return observed_count, pval, gene_start, gene_end


def tajima_label_perm(nperm=5000):
    td = pd.read_csv(TD_FILE)
    td_clean = td.dropna(subset=['tajimas_d']).reset_index(drop=True)
    obs_overlap = td_clean[td_clean['overlap_high_delta'] == True]['tajimas_d'].values
    obs_non = td_clean[td_clean['overlap_high_delta'] == False]['tajimas_d'].values

    obs_diff = np.nanmean(obs_overlap) - np.nanmean(obs_non)

    diffs = np.empty(nperm)
    n = len(td_clean)
    labels = td_clean['overlap_high_delta'].values
    for i in range(nperm):
        perm = np.random.permutation(labels)
        perm_overlap = td_clean['tajimas_d'].values[perm == True]
        perm_non = td_clean['tajimas_d'].values[perm == False]
        if len(perm_overlap) == 0 or len(perm_non) == 0:
            diffs[i] = 0.0
        else:
            diffs[i] = np.nanmean(perm_overlap) - np.nanmean(perm_non)

    pval = (np.sum(np.abs(diffs) >= np.abs(obs_diff)) + 1) / (nperm + 1)

    pd.DataFrame({'perm_diff': diffs}).to_csv(OUT / 'perm_tajimas_diff.csv', index=False)

    plt.figure(figsize=(6,4))
    sns.histplot(diffs, bins=50, kde=False)
    plt.axvline(obs_diff, color='red', linestyle='--', label=f'observed diff={obs_diff:.4f}')
    plt.title(f'Permutation: mean(Tajima D overlap) - mean(Tajima D non-overlap) (n={nperm})')
    plt.xlabel('mean difference')
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / 'perm_tajimas_diff.png')
    plt.close()

    return obs_diff, pval


def append_results_to_md(max_res, gene_res, taj_res):
    md = OUT / 'benchmark_report.md'
    now = time.strftime('%Y-%m-%d %H:%M:%S')
    lines = []
    lines.append('\n')
    lines.append('## Permutation-based significance tests')
    lines.append('\n')
    lines.append(f'_Tests run: {now}_')
    lines.append('\n')
    lines.append('### 1) Max-FST of top 20 sites')
    lines.append('\n')
    observed_max, pval_max, elapsed = max_res
    lines.append(f'- Observed max FST among top 20: **{observed_max:.4f}**')
    lines.append(f'- Permutation p-value (how often random sets of 20 produce max >= observed): **{pval_max:.4f}**')
    lines.append(f'- Plot: `results/perm_max_fst_top20.png`')
    lines.append('\n')
    lines.append('### 2) Gene-overlap of top 20 sites')
    lines.append('\n')
    obs_count, pval_gene, gstart, gend = gene_res
    if gstart is not None:
        lines.append(f'- Observed count of top-20 sites inside `ABO` gene region ({gstart}-{gend}): **{obs_count}**')
    else:
        lines.append(f'- Observed count (gene bounds not available): **{obs_count}**')
    lines.append(f'- Permutation p-value (how often random sets have ≥ observed count): **{pval_gene:.4f}**')
    lines.append(f'- Plot: `results/perm_gene_overlap_top20.png`')
    lines.append('\n')
    lines.append('### 3) Tajima\'s D: windows overlapping high ΔAF vs others')
    lines.append('\n')
    obs_diff, pval_td = taj_res
    lines.append(f'- Observed mean difference (overlap - non-overlap): **{obs_diff:.4f}**')
    lines.append(f'- Permutation p-value (two-sided): **{pval_td:.4f}**')
    lines.append(f'- Plot: `results/perm_tajimas_diff.png`')
    lines.append('\n')
    lines.append('### Interpretation')
    lines.append('\n')
    lines.append('- Small permutation p-values (e.g., < 0.05) indicate that the observed statistic is unlikely under random draws from the region and therefore may be indicative of a localized signal rather than chance.\n')
    lines.append('- Remember: these tests are **local** to the ABO region and compare against random draws from the same site set; for genome-wide significance, compare to genome-wide matched backgrounds.\n')

    with md.open('a') as f:
        f.write('\n---\n')
        f.write('\n'.join(lines))

    return md


if __name__ == '__main__':
    sns.set(style='whitegrid')
    print('Running max-FST permutation...')
    max_res = max_fst_perm(k=20, nperm=5000)
    print('Running gene-overlap permutation...')
    gene_res = gene_overlap_perm(k=20, nperm=5000)
    print('Running Tajima\'s D permutation...')
    taj_res = tajima_label_perm(nperm=5000)

    print('Appending results to benchmarkanalysis.md')
    md = append_results_to_md(max_res, gene_res, taj_res)
    print('Wrote permutation outputs and updated', md)
