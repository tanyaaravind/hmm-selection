### Project
Applying Hidden Markov Models to Detect Selection on ABO Blood Group Genes Across Human Populations

### Team
- Wellington Mapise (wpm44), Junior, Computer Science
- Tanya Aravind (ta374), Junior, Computer Science
- Amy Oduor (aao73), Junior, Computer Science
- Sheki Okwayo (lso24), Junior, Computer Science

### Scope
- Default analysis: 2-state HMM (neutral vs selection) on the ABO locus only. Balancing-selection, Rh locus, or full Shen et al. non-homogeneous features are stretch goals, not required.
- Emissions: one scalar population-differentiation metric per SNP. Defaults to ΔAF (e.g., |AF_YRI−AF_CEU|); FST can be swapped in. Set the Gaussian means/stds in `hmm_core.py`/scripts to match whichever metric you choose.
- Heavyweight stats (full iHS/Tajima’s D pipelines) can be pulled from existing tools; `src/stats.py` keeps only lightweight helpers. Use external software if time is short.

### Reporting reminders
- Clearly state the two hidden states (neutral, selection) and the exact per-SNP emission metric used.
- If you change emissions (ΔAF vs FST), log the parameter choices for mean/std in the report so readers can interpret the HMM outputs.

### Overview
- Two-state HMM in `src/hmm_core.py` detects neutral vs selection signals from differentiation metrics. Forward/backward, Viterbi decoding, and Baum–Welch training are implemented.
- Simulations in `src/simulation.py` generate toy SNP series with a known selection segment.
- Visualization scripts in `src/visualization.py` and `tests/test_abo_data.py` produce figures saved in `results/`.
- Data prep stub in `src/data_prep.py` extracts allele frequencies from 1000 Genomes VCF slices and computes ΔAF tables for specified populations.

### Quickstart
1) Install deps: `pip install -r requirements.txt` (optionally `pip install cyvcf2` for VCF parsing).
2) Run toy sanity check: `python tests/test_toy_example.py`.
3) Minimal fit/decoding check: `python tests/test_hmm_fit.py`.
4) Generate simulated figure: `python src/visualization.py` → `results/preliminary_results.png`.
5) ABO example figure (hardcoded ΔAF table): `python tests/test_abo_data.py` → `results/real_abo_analysis.png`.

### Training (Baum–Welch)
Use `SelectionHMM.fit(observations, positions, n_iter=10)` to re-estimate emission means/stds and a heuristic distance scale. Example:
```python
from hmm_core import SelectionHMM
hmm = SelectionHMM(emission_params, transition_params)
hmm.fit(observations, positions, n_iter=10)
```

### Data Prep
- Input: bgzipped VCF slice for ABO (Rh optional) and a sample map CSV (`sample_id,pop`).
- Run:
```bash
python src/data_prep.py --vcf data/abo_slice.vcf.gz --sample-map data/samples.csv --pop-a YRI --pop-b CEU --out results/abo_freqs
```
- Outputs: `allele_frequencies.csv` and `delta_af.csv` under the chosen directory. If VCF parsing is unavailable, supply `--af-table` with a precomputed allele frequency CSV.

### Notes
- `src/stats.py` is reserved for benchmark metrics (FST/Tajima’s D/iHS) if we add them; otherwise use external tools and just compare.
- Plots require matplotlib/seaborn; headless environments may need `MPLBACKEND=Agg`.

### Benchmark helpers (Sheki)
- `src/stats.py` now provides lightweight baselines: per-site FST (with counts or frequencies), Tajima's D from allele counts within a window, and a toy iHS that consumes precomputed EHH curves.
- Example FST on a tidy allele-frequency table: 
  `python - <<'PY'\nimport pandas as pd\nfrom src.stats import fst_from_af_table\n\naf = pd.read_csv('data/example_abo_af.csv')\nfst = fst_from_af_table(af, 'YRI', 'CEU')\nprint(fst.head())\nPY`
- Tajima's D expects allele counts (2N haplotypes) per site within a window: 
  `from src.stats import tajimas_d_from_counts`
