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
  - `src/visualization.py`: Generates 2-panel preliminary results figure from simulated data → `results/preliminary_results.png`
  - `src/visualization_expanded.py`: Generates 2-panel preliminary results figure for real ABO data → `results/preliminary_results_expanded.png`
  - `tests/test_abo_data.py`: Generates comprehensive 3-panel figures for both small and expanded datasets → `results/real_abo_analysis.png` and `results/real_abo_analysis_expanded.png`
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
cd selection_hmm

# Toy example with known ground truth
python tests/test_toy_example.py

# Minimal fit/decoding check
python tests/test_hmm_fit.py
```

#### 3. Generate Preliminary Results Figures
```bash
# Simulated data (2-panel figure)
python src/visualization.py
# Output: results/preliminary_results.png

# Real ABO data (2-panel figure)
python src/visualization_expanded.py
# Output: results/preliminary_results_expanded.png
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

# Step 3: Generate comprehensive 3-panel analysis figures
python tests/test_abo_data.py
# Outputs:
#   - results/real_abo_analysis.png (small dataset, 8 SNPs)
#   - results/real_abo_analysis_expanded.png (expanded dataset, 648+ SNPs)
```


#### 5. Benchmarking & report

Run the commands below to reproduce the deliverables. All outputs go to `results/`.

- Install dependencies:
  ```bash
  pip install -r requirements.txt
  # (requests, seaborn, cyvcf2 optional)
  ```

- Quick run (FST + Tajima's D):
  ```bash
  python scripts/run_benchmarks.py
  # -> results/benchmarks_fst_yri_ceu.csv, results/benchmarks_tajimas_yri.csv
  ```

- Full run (FST for all pairs, plots, sliding Tajima's D, analysis):
  ```bash
  python scripts/benchmark_and_plot.py
  # -> CSVs, plots, and results/benchmark_analysis.md
  ```

- Annotate top sites:
  ```bash
  python scripts/annotate_top_fst.py --fst results/benchmarks_fst_yri_ceu.csv --top 20 --chrom 9
  # -> results/benchmarks_fst_yri_ceu_top20_annotated.csv
  ```

- report (already generated): `results/benchmark_report.md` (includes key plots and permutation summaries).

All scripts are in `scripts/` and primary outputs are in `results/`.

### Quick grader checklist (recommended)

Run the minimal reproducible pipeline:

```bash
make benchmarks   # runs full benchmark_and_plot pipeline
make report       # generates concise markdown report (results/benchmark_report.md)
make verify       # quick check that key outputs exist
```

All outputs are saved in `results/` and the concise report is `results/benchmark_report.md`.

Notes:
- The annotation script falls back gracefully if Ensembl REST is unreachable (`annotation_status=api-failed` in the CSV). If you prefer offline annotation, point me to a local GTF and I can add GTF-based annotation.
- All generated CSVs and images are saved under `results/`.

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

### Visualization Scripts

**`src/visualization.py`** - Simulated data visualizations:
- 2-panel preliminary figure (FST values + HMM posteriors)
  - Output: `results/preliminary_results.png`
- 3-panel synthetic results figure (FST values + true vs predicted states + HMM posteriors)
  - Output: `results/synthetic_results.png`

**`src/visualization_expanded.py`** - Real data preliminary results:
- Generates 2-panel figure (DeltaAF values + HMM posteriors)
- Uses real ABO data from 1000 Genomes
- Output: `results/preliminary_results_expanded.png`

**`tests/test_abo_data.py`** - Comprehensive real data analysis:
- Generates comprehensive 3-panel figures for both datasets:
  1. Allele frequencies across populations (YRI, CEU, CHB)
  2. DeltaAF with highlighting for high-differentiation regions
  3. HMM posterior probabilities of selection state
- Outputs:
  - `results/real_abo_analysis.png` (small dataset, 8 SNPs)
  - `results/real_abo_analysis_expanded.png` (expanded dataset, 648+ SNPs)

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
│   ├── abo_freqs/               # Initial 8-SNP dataset
│   │   ├── allele_frequencies.csv
│   │   └── delta_af.csv
│   ├── abo_expanded/            # Expanded ABO dataset
│   │   ├── allele_frequencies.csv
│   │   ├── delta_af.csv
│   │   └── delta_af_nonzero.csv
│   ├── preliminary_results.png  # Simulated data (2-panel)
│   ├── preliminary_results_expanded.png  # Real data (2-panel)
│   ├── real_abo_analysis.png  # Small dataset comprehensive (3-panel)
│   └── real_abo_analysis_expanded.png  # Expanded dataset comprehensive (3-panel)
└── filter_nonzero_deltaaf.py    # Filter script for non-zero DeltaAF
```

### Reproducibility

This repository is self-contained and all results can be reproduced following the steps below. All code, processed data files, and analysis scripts are included.

#### Dataset Explanation & Data Sources

**All allele-frequency data comes from the 1000 Genomes Project Phase 3.** We use two datasets that fit into our analysis pipeline:

**Dataset 1: Initial 8-SNP Test Dataset (Testing Phase)**
- **Source**: 1000 Genomes Project Phase 3, accessed via **Ensembl REST API** during initial testing
- **API Base URL**: https://rest.ensembl.org/
- **API Endpoints Used**:
  - `GET /overlap/region/human/{region}` - to identify SNPs in ABO region
  - `GET /variation/human/{rsID}?pops=1` - to retrieve allele frequencies for specific rsIDs
- **SNPs**: 8 key ABO blood group SNPs:
  - rs8176749, rs8176746, rs8176747, rs8176743, rs8176740, rs8176719, rs687289, rs505922
- **Region**: chr9:133,255,801 - 133,273,813 (~18 kb)
- **Populations**: YRI (African), CEU (European), CHB (East Asian)
- **Purpose**: Initial validation and testing of HMM implementation
- **Status**: Data is hardcoded in `tests/test_abo_data.py` for full reproducibility (no API calls needed)
- **Reproducibility**: The exact allele frequencies are embedded in the code, so results are 100% reproducible
- **Output**: `results/real_abo_analysis.png` (3-panel comprehensive figure)

**Dataset 2: Expanded ABO Region Dataset (Full Analysis)**
- **Source**: 1000 Genomes Project Phase 3 (20130502 release)
- **Repository**: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
- **VCF File**: `ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
- **Sample Metadata**: `integrated_call_samples_v3.20130502.ALL.panel`
- **Region Analyzed**: Chromosome 9, ABO locus (chr9:133,240,000-133,290,000, ~50 kb)
- **Populations**: YRI (African), CEU (European), CHB (East Asian)
- **SNPs**: 1831 total variants (648 non-zero DeltaAF after filtering)
- **Purpose**: Comprehensive analysis with higher SNP density
- **Processing**: Extracted from VCF using `src/data_prep.py`
- **Output**: `results/real_abo_analysis_expanded.png` (3-panel comprehensive figure)

**How Datasets Fit Into Pipeline:**
1. **Testing Phase**: Used 8-SNP dataset (Dataset 1) to validate HMM on small, curated set
2. **Full Analysis**: Expanded to 1831 SNPs (Dataset 2) for comprehensive region-wide analysis
3. Both datasets use same populations (YRI, CEU, CHB) and same analysis methods
4. Results from both are included in final figures for comparison

**Note on Data Files:**
- Large VCF files (>600MB) are excluded from the repository (see `.gitignore`)
- Processed data files are included: `results/abo_expanded/` contains all processed CSV files
- Sample map (`data/samples.csv`) is included for immediate use
- Initial 8-SNP dataset is available in `results/abo_freqs/` CSV files and also hardcoded in `tests/test_abo_data.py`
- To regenerate expanded dataset from raw VCF, follow the data preparation steps below

#### Complete Reproduction Pipeline

**Step 1: Install Dependencies**
```bash
pip install -r requirements.txt
conda install -c bioconda cyvcf2  # For VCF parsing
```

**Step 2: Validate HMM Implementation (Simulated Data)**
```bash
# Test with known ground truth
python tests/test_toy_example.py
python tests/test_hmm_fit.py

# Generate simulated data figure
python src/visualization.py
# Output: results/preliminary_results.png
```

**Step 3: Generate Real Data Analysis Figures**

The 8-SNP dataset is hardcoded and requires no data download. The expanded dataset uses processed CSV files included in the repository.

```bash
# Generate comprehensive 3-panel figures for both datasets
python tests/test_abo_data.py
# Outputs:
#   - results/real_abo_analysis.png (8-SNP dataset, hardcoded)
#   - results/real_abo_analysis_expanded.png (expanded dataset, from processed CSV files)
```

**Step 4: Regenerate Expanded Dataset (Optional - if you need to process from raw VCF)**

If you need to regenerate the expanded dataset from raw VCF files (processed files are already included):

```bash
# 4a. Download 1000 Genomes VCF (if not already present)
# Full chromosome 9 VCF (~614MB):
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

# 4b. Download sample metadata panel file
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# 4c. Extract ABO region using bcftools
bcftools view -r 9:133240000-133290000 \
    data/raw_vcf/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
    -Oz -o data/abo_region.vcf.gz
bcftools index data/abo_region.vcf.gz

# 4d. Create sample map from panel file
python scripts/create_sample_map.py \
    --panel data/integrated_call_samples_v3.20130502.ALL.panel \
    --output data/samples.csv

# 4e. Extract allele frequencies and compute DeltaAF
python src/data_prep.py \
    --vcf data/abo_region.vcf.gz \
    --sample-map data/samples.csv \
    --pop-a YRI \
    --pop-b CEU \
    --out results/abo_expanded

# 4f. Filter to non-zero DeltaAF (recommended)
python filter_nonzero_deltaaf.py \
    --input results/abo_expanded/delta_af.csv \
    --output results/abo_expanded/delta_af_nonzero.csv
```

**Step 5: Generate All Figures**

```bash
# 5a. Preliminary results (2-panel figures)
python src/visualization.py
# Output: results/preliminary_results.png (simulated data)

python src/visualization_expanded.py
# Output: results/preliminary_results_expanded.png (real data, expanded dataset)

# 5b. Comprehensive analysis (3-panel figures)
python tests/test_abo_data.py
# Outputs:
#   - results/real_abo_analysis.png (8-SNP dataset, hardcoded)
#   - results/real_abo_analysis_expanded.png (expanded dataset, 648+ SNPs from processed CSV)
```

#### Included Processed Data

The following processed data files are included in the repository for immediate analysis:

**Dataset 1 (8-SNP):**
- **Files**: `results/abo_freqs/allele_frequencies.csv` and `results/abo_freqs/delta_af.csv`
- **Also available**: `data/example_abo_af.py` 
- **Content**: Allele frequencies for 8 SNPs across YRI, CEU, CHB populations
- **Reproducibility**: CSV files included in repository, fully reproducible without data downloads

**Dataset 2 (Expanded):**
- `results/abo_expanded/allele_frequencies.csv`: Allele frequencies for all populations (YRI, CEU, CHB) - 1831 SNPs
- `results/abo_expanded/delta_af.csv`: DeltaAF values for all 1831 SNPs
- `results/abo_expanded/delta_af_nonzero.csv`: Filtered to 648 non-zero DeltaAF SNPs (recommended for analysis)
- `data/samples.csv`: Sample ID to population mapping (2506 samples from 1000 Genomes)

**Note on Excluded Files:**
- Large raw VCF files (>600MB) are excluded due to size (see `.gitignore`)
- The processed CSV files contain all necessary data for reproducing analyses
- To regenerate expanded dataset from raw VCF, use the commands in Step 4 above
- The 8-SNP dataset requires no external files - it's fully self-contained in the code

#### Expected Outputs

Running the complete pipeline will generate:
- `results/preliminary_results.png`: Simulated data validation (2-panel)
- `results/preliminary_results_expanded.png`: Real data preliminary results (2-panel)
- `results/real_abo_analysis.png`: Small dataset comprehensive analysis (3-panel, 8 SNPs)
- `results/real_abo_analysis_expanded.png`: Expanded dataset comprehensive analysis (3-panel, 648 SNPs)

All figures are saved at 300 DPI and can be reproduced exactly by following the steps above.

### Notes
- `src/stats.py` is reserved for benchmark metrics (FST/Tajima's D/iHS) if we add them; otherwise use external tools and just compare.
- Plots require matplotlib/seaborn; headless environments may need `MPLBACKEND=Agg` or `matplotlib.use('Agg')`.
