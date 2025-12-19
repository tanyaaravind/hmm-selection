### Project
Applying Hidden Markov Models to Detect Selection on ABO Blood Group Genes Across Human Populations

### Team
- Wellington Mapise (wpm44), Junior, Computer Science
- Tanya Aravind (ta374), Junior, Computer Science
- Amy Oduor (aao73), Junior, Computer Science
- Sheki Okwayo (lso24), Junior, Computer Science

### Scope
- We are scoping to a 2-state HMM (neutral vs selection) on the ABO locus as the primary result. Balancing-selection state, Rh locus, and full Shen et al. non-homogeneous details are optional/future work.
- Emissions: per-SNP population differentiation metrics. Current scripts use ΔAF (absolute difference in allele frequencies) between two populations (e.g., YRI vs CEU); `hmm_core.py` assumes Gaussian emissions over these values. If you swap in FST or another scalar, update the means/stds accordingly.
- Complex stats (iHS/Tajima's D/FST computation) can rely on existing tools; `src/stats.py` is a placeholder if needed.

### Overview
- Two-state HMM in `src/hmm_core.py` detects neutral vs selection signals from differentiation metrics. Forward/backward, Viterbi decoding, and Baum–Welch training are implemented.
- Simulations in `src/simulation.py` generate toy SNP series with a known selection segment.
- Visualization scripts:
  - `src/visualization.py`: Generates preliminary results figure from simulated data → `results/preliminary_results.png`
  - `src/visualization_expanded.py`: Main script for real ABO data analysis with HMM → `results/real_abo_analysis_expanded.png`
  - `tests/test_abo_data.py`: Test script with hardcoded small dataset → `results/real_abo_analysis.png`
- Data prep in `src/data_prep.py` extracts allele frequencies from 1000 Genomes VCF slices and computes ΔAF tables for specified populations (YRI, CEU, CHB).

### Quickstart

#### 1. Install Dependencies
```bash
pip install -r requirements.txt
# For VCF parsing (required for data extraction):
pip install cyvcf2
# Or with conda:
conda install -c bioconda cyvcf2
```

#### 2. Run Validation Tests
```bash
# Toy example with known ground truth
python tests/test_toy_example.py

# Minimal fit/decoding check
python tests/test_hmm_fit.py
```

#### 3. Generate Simulated Data Figure
```bash
python src/visualization.py
# Output: results/preliminary_results.png
```

#### 4. Run Real ABO Data Analysis (Main Workflow)
```bash
# Step 1: Extract data from VCF (if you have a VCF file)
python src/data_prep.py \
    --vcf data/abo_region.vcf.gz \
    --sample-map data/samples.csv \
    --pop-a YRI \
    --pop-b CEU \
    --out results/abo_expanded

# Step 2: Filter to non-zero DeltaAF (optional, recommended)
python filter_nonzero_deltaaf.py \
    --input results/abo_expanded/delta_af.csv \
    --output results/abo_expanded/delta_af_nonzero.csv

# Step 3: Run HMM analysis and generate visualization
python src/visualization_expanded.py
# Output: results/real_abo_analysis_expanded.png
```

### Data Preparation

#### Input Requirements
- **VCF file**: bgzipped VCF slice for ABO region (e.g., `data/abo_region.vcf.gz`)
- **Sample map**: CSV file mapping sample IDs to populations (`sample_id,pop`)
  - Can be created from 1000 Genomes panel file using `scripts/create_sample_map.py`

#### Running Data Prep
```bash
python src/data_prep.py \
    --vcf data/abo_region.vcf.gz \
    --sample-map data/samples.csv \
    --pop-a YRI \
    --pop-b CEU \
    --out results/abo_expanded
```

#### Outputs
- `allele_frequencies.csv`: Allele frequencies per population (YRI, CEU, CHB if available)
- `delta_af.csv`: DeltaAF values for all SNPs (includes zeros)
- `delta_af_nonzero.csv`: Filtered version with only non-zero DeltaAF values (recommended for analysis)

#### Filtering Non-Zero DeltaAF
Many SNPs have zero DeltaAF (no population differentiation). Filtering removes these for cleaner analysis:
```bash
python filter_nonzero_deltaaf.py \
    --input results/abo_expanded/delta_af.csv \
    --output results/abo_expanded/delta_af_nonzero.csv
```

### Main Analysis Script

**`src/visualization_expanded.py`** is the primary script for real data analysis:
- Loads expanded ABO data (prefers filtered non-zero version)
- Runs HMM with Forward-Backward algorithm
- Generates comprehensive 3-panel figure:
  1. Allele frequencies across populations (YRI, CEU, CHB)
  2. DeltaAF with highlighting for high-differentiation regions
  3. HMM posterior probabilities of selection state
- Saves to `results/real_abo_analysis_expanded.png`

### Training (Baum–Welch)
Use `SelectionHMM.fit(observations, positions, n_iter=10)` to re-estimate emission means/stds and a heuristic distance scale. Example:
```python
from hmm_core import SelectionHMM
hmm = SelectionHMM(emission_params, transition_params)
hmm.fit(observations, positions, n_iter=10)
```

### Project Structure
```
hmm-selection/
├── src/
│   ├── hmm_core.py              # Core HMM implementation
│   ├── data_prep.py             # VCF parsing and DeltaAF computation
│   ├── simulation.py            # Simulated data generation
│   ├── visualization.py         # Simulated data visualization
│   └── visualization_expanded.py # Main real data analysis script (includes visualization functions)
├── tests/
│   ├── test_toy_example.py      # Validation with known ground truth
│   ├── test_hmm_fit.py          # HMM fitting tests
│   └── test_abo_data.py         # Small real data test
├── scripts/
│   ├── create_sample_map.py     # Convert 1000 Genomes panel to sample map
│   └── check_data_size.py       # Check number of SNPs in CSV
├── data/                        # Input data (VCF files, sample maps)
├── results/                     # Output figures and processed data
│   ├── abo_expanded/            # Expanded ABO dataset
│   │   ├── allele_frequencies.csv
│   │   ├── delta_af.csv
│   │   └── delta_af_nonzero.csv
│   ├── preliminary_results.png  # Simulated data figure
│   └── real_abo_analysis_expanded.png  # Main real data figure
└── filter_nonzero_deltaaf.py    # Filter script for non-zero DeltaAF
```

### Notes
- `src/stats.py` is reserved for benchmark metrics (FST/Tajima's D/iHS) if we add them; otherwise use external tools and just compare.
- Plots require matplotlib/seaborn; headless environments may need `MPLBACKEND=Agg` or `matplotlib.use('Agg')`.
