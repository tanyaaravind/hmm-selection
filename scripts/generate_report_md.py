#!/usr/bin/env python3
"""Generate a concise Markdown report for graders from analysis artifacts.

Produces:
- results/benchmark_report.md

Includes:
- Short summary (Methods + Conclusion excerpt)
- Embedded plots (as relative links)
- Top annotated sites table (first 20 rows)
- Permutation test summaries (extracted from benchmarkanalysis.md)
- Links to detailed files
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
MD_SRC = OUT / 'benchmark_report.md'
OUT_MD = OUT / 'benchmark_report.md'

plots = [
    ('FST comparison (violin)', 'fst_comparison.png'),
    ('FST per-site (YRI vs CEU)', 'fst_vs_pos_yri_ceu.png'),
    ("Tajima's D (YRI) sliding windows", 'tajimas_yri_sliding_w50_s10.png'),
    ('Permutation: max FST (top 20)', 'perm_max_fst_top20.png'),
    ("Permutation: Tajima's D diff", 'perm_tajimas_diff.png'),
]

ann_csv = OUT / 'benchmarks_fst_yri_ceu_top20_annotated.csv'

# Read benchmarkanalysis.md and extract Methods, Conclusion, Permutation section
text = MD_SRC.read_text() if MD_SRC.exists() else ''

def excerpt_section(text, header):
    if header not in text:
        return ''
    part = text.split(header,1)[1]
    # stop at next '## '
    if '\n## ' in part:
        part = part.split('\n## ',1)[0]
    return part.strip()

methods = excerpt_section(text, '## Methods (concise)')
conclusion = excerpt_section(text, '## Conclusion')
perms = excerpt_section(text, '## Permutation-based significance tests')

lines = []
lines.append('# ABO Benchmark — Concise Report')
lines.append('')
lines.append('**Short summary:** This report contains key visuals, top annotated sites, and results of quick permutation-based tests to assess whether observed signals in the ABO region are unlikely by chance. The pipeline is intentionally focused (2-state HMM on ABO; complementary statistics: ΔAF, FST, Tajima\'s D).')
lines.append('')

if methods:
    lines.append('## Methods (brief)')
    lines.append('')
    lines.append(methods)
    lines.append('')

lines.append('## Key visuals')
lines.append('')
for caption, fname in plots:
    fpath = Path(fname)
    if (OUT / fname).exists():
        lines.append(f'### {caption}')
        lines.append(f'![]({fname})')
        lines.append('')
    else:
        lines.append(f'- Missing plot: `{fname}`')
        lines.append('')

lines.append('## Top annotated FST sites (YRI vs CEU)')
lines.append('')
if ann_csv.exists():
    ann = pd.read_csv(ann_csv)
    cols = [c for c in ['pos','fst','gene_name','annotation_status','known_rs'] if c in ann.columns]
    ann = ann[cols].head(20)
    # convert to markdown table
    lines.append(ann.to_markdown(index=False))
else:
    lines.append('No annotation CSV found.')
lines.append('')

lines.append('## Permutation test summaries (brief)')
lines.append('')
if perms:
    lines.append(perms)
else:
    lines.append('No permutation results found.')
lines.append('')

lines.append('## Conclusion (brief)')
lines.append('')
if conclusion:
    lines.append(conclusion)
else:
    lines.append('No conclusion found.')
lines.append('')

lines.append('## Files & Reproducibility')
lines.append('')
lines.append('- Full write-up and details: included in `results/benchmark_report.md` (extended analysis appended)')
lines.append('- All CSVs and plots are in `results/` (use scripts in `scripts/` to reproduce)')
lines.append('')
lines.append('---')
lines.append('Generated automatically from analysis artifacts.')

OUT_MD.write_text('\n\n'.join(lines))
print('Wrote concise Markdown report to', OUT_MD)
