### Project
Applying Hidden Markov Models to Detect Selection on ABO Blood Group Genes Across Human Populations

### Team
- Wellington Mapise (wpm44), Junior, Computer Science
- Tanya Aravind (ta374), Junior, Computer Science
- Amy Oduor (aao73), Junior, Computer Science
- Sheki Okwayo (lso24), Junior, Computer Science

### Scope
- We are scoping to a 2-state HMM (neutral vs selection) on the ABO locus as the primary result. Balancing-selection state, Rh locus, and full Shen et al. non-homogeneous details are optional/future work.
- Emissions: per-SNP population differentiation metrics. Current scripts use ΔAF between two populations (e.g., YRI vs CEU); `hmm_core.py` assumes Gaussian emissions over these values. If you swap in FST or another scalar, update the means/stds accordingly.
- Complex stats (iHS/Tajima’s D/FST computation) can rely on existing tools; `src/stats.py` is a placeholder if needed.

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
