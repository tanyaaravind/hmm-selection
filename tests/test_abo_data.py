
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
src_path = os.path.join(project_root, 'src')
sys.path.append(src_path)

from hmm_core import SelectionHMM

sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['font.size'] = 10


def analyze_real_abo_data():
    print("="*70)
    print("REAL DATA ANALYSIS: ABO Blood Type Locus")
    print("="*70)

    snp_ids = [
        'rs8176749', 'rs8176746', 'rs8176747', 'rs8176743',
        'rs8176740', 'rs8176719', 'rs687289', 'rs505922'
    ]
    
    positions = np.array([
        133255801, 133255935, 133255936, 133256028,
        133256085, 133257521, 133261703, 133273813
    ])
    
    af_yri = np.array([0.829, 0.829, 0.829, 0.829, 0.806, 0.769, 0.606, 0.690])
    af_ceu = np.array([0.929, 0.929, 0.929, 0.929, 0.722, 0.621, 0.621, 0.631])
    af_chb = np.array([0.772, 0.772, 0.772, 0.772, 0.767, 0.607, 0.602, 0.602])
    
    observations = np.array([0.101, 0.101, 0.101, 0.101, 0.083, 0.147, 0.015, 0.059])
    
    print("\nDataset Summary:")
    print(f"  Number of SNPs: {len(snp_ids)}")
    print(f"  Genomic region: chr9:{positions[0]:,} - {positions[-1]:,}")
    print(f"  Region size: {positions[-1] - positions[0]:,} bp (~18 kb)")
    print(f"  Populations: YRI (African), CEU (European), CHB (East Asian)")
    
    print("\n" + "-"*70)
    print("SNP Data:")
    print("-"*70)
    print(f"{'SNP ID':<15} {'Position':>12} {'AF_YRI':>8} {'AF_CEU':>8} {'AF_CHB':>8} {'DeltaAF':>8}")
    print("-"*70)
    for i in range(len(snp_ids)):
        print(f"{snp_ids[i]:<15} {positions[i]:>12,} {af_yri[i]:>8.3f} {af_ceu[i]:>8.3f} "
              f"{af_chb[i]:>8.3f} {observations[i]:>8.3f}")
    print("-"*70)
    
    print("\n" + "="*70)
    print("BIOLOGICAL PATTERN:")
    print("="*70)
    print("Notice the DeltaAF pattern:")
    print("  • SNPs 1-4 (rs8176749-rs8176743): DeltaAF = 0.101 (consistently high)")
    print("  • SNP 5 (rs8176740): DeltaAF = 0.083 (moderate)")
    print("  • SNP 6 (rs8176719): DeltaAF = 0.147 ⭐ PEAK - strongest signal!")
    print("  • SNP 7 (rs687289): DeltaAF = 0.015 (very low - neutral-like)")
    print("  • SNP 8 (rs505922): DeltaAF = 0.059 (moderate)")
    print("\nThis pattern suggests SELECTION in the first part of the region,")
    print("with a peak at rs8176719, followed by return to neutrality.")
    
    print("\n" + "="*70)
    print("RUNNING HMM ANALYSIS")
    print("="*70)
    

    emission_params = {
        'neutral': {'mean': 0.04, 'std': 0.03},
        'selection': {'mean': 0.12, 'std': 0.03}
    }
    
    transition_params = {
        'distance_scale': 2000
    }
    
    hmm = SelectionHMM(emission_params, transition_params)
    hmm.fit(observations, positions, n_iter=10, verbose=False)
    
    posteriors = hmm.posterior_probabilities(observations, positions)
    
    print("\n" + "="*70)
    print("HMM RESULTS: Posterior Probabilities")
    print("="*70)
    print(f"{'SNP ID':<15} {'Position':>12} {'DeltaAF':>8} {'P(Neutral)':>12} {'P(Selection)':>13} {'State':>12}")
    print("-"*70)
    
    for i in range(len(snp_ids)):
        pred_state = 'Selection' if posteriors[i, 1] > 0.5 else 'Neutral'
        print(f"{snp_ids[i]:<15} {positions[i]:>12,} {observations[i]:>8.3f} "
              f"{posteriors[i, 0]:>12.3f} {posteriors[i, 1]:>13.3f} {pred_state:>12}")
    
    print("="*70)
    
    print("\nGenerating visualization...")
    fig = create_real_data_figure(snp_ids, positions, observations, posteriors, af_yri, af_ceu, af_chb)
    
    save_path = os.path.join(project_root, 'results', 'real_abo_analysis.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to: {save_path}")
    
    return posteriors


def create_real_data_figure(snp_ids, positions, deltaaf, posteriors, af_yri, af_ceu, af_chb):

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1.2], hspace=0.35)
    
    pos_kb = (positions - positions[0]) / 1000
    ax1 = fig.add_subplot(gs[0])
    
    ax1.plot(pos_kb, af_yri, 'o-', color='#E64B35', linewidth=2, markersize=8, 
             label='YRI (African)', alpha=0.8)
    ax1.plot(pos_kb, af_ceu, 's-', color='#4DBBD5', linewidth=2, markersize=8, 
             label='CEU (European)', alpha=0.8)
    ax1.plot(pos_kb, af_chb, '^-', color='#00A087', linewidth=2, markersize=8, 
             label='CHB (East Asian)', alpha=0.8)
    
    ax1.set_ylabel('Allele Frequency', fontsize=11, fontweight='bold')
    ax1.set_title('Panel A: Allele Frequencies Across Three Populations', 
                  fontsize=12, fontweight='bold')
    ax1.legend(loc='best', fontsize=10, framealpha=0.9)
    ax1.grid(alpha=0.3)
    ax1.set_ylim(0.5, 1.0)
    
    ax2 = fig.add_subplot(gs[1])
    
    ax2.plot(pos_kb, deltaaf, 'o-', color='purple', linewidth=2.5, 
             markersize=10, label='ΔAF = |AF_YRI - AF_CEU|', alpha=0.8)
    
    peak_idx = 5
    ax2.plot(pos_kb[peak_idx], deltaaf[peak_idx], 'o', color='red', 
             markersize=15, alpha=0.5, label='rs8176719 (peak)')
    
    ax2.axhline(0.04, color='green', linestyle='--', linewidth=1.5, 
                alpha=0.6, label='Neutral-like (<0.06)')
    ax2.axhline(0.10, color='red', linestyle='--', linewidth=1.5, 
                alpha=0.6, label='Selection-like (>0.10)')
    
    ax2.set_ylabel('ΔAF (Population Differentiation)', fontsize=11, fontweight='bold')
    ax2.set_title('Panel B: Population Differentiation Signal', 
                  fontsize=12, fontweight='bold')
    ax2.legend(loc='best', fontsize=9, framealpha=0.9)
    ax2.grid(alpha=0.3)
    ax2.set_ylim(0, 0.18)
    
    ax3 = fig.add_subplot(gs[2])
    ax3.plot(pos_kb, posteriors[:, 1], '-', color='red', linewidth=3, 
             label='P(Selection | Data)')
    ax3.fill_between(pos_kb, 0, posteriors[:, 1], alpha=0.3, color='red')

    ax3.plot(pos_kb, posteriors[:, 1], 'o', color='darkred', markersize=10, alpha=0.7)
    
    ax3.axhline(0.5, color='gray', linestyle='--', linewidth=1.5, 
                alpha=0.7, label='Decision threshold (0.5)')
    
    for i, (x, y, snp) in enumerate(zip(pos_kb, posteriors[:, 1], snp_ids)):
        if y > 0.5:
            ax3.annotate(snp, xy=(x, y), xytext=(0, 10), 
                        textcoords='offset points', fontsize=8,
                        ha='center', rotation=45, alpha=0.7)
    
    ax3.set_xlabel('Position Relative to First SNP (kb)', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Posterior Probability', fontsize=11, fontweight='bold')
    ax3.set_title('Panel C: HMM Detection of Selection Signal', 
                  fontsize=12, fontweight='bold')
    ax3.legend(loc='best', fontsize=10, framealpha=0.9)
    ax3.grid(alpha=0.3)
    ax3.set_ylim(-0.05, 1.05)
    
    n_selection = sum(posteriors[:, 1] > 0.5)
    summary_text = f'SNPs identified as under selection: {n_selection}/8'
    ax3.text(0.02, 0.95, summary_text, transform=ax3.transAxes,
             fontsize=11, fontweight='bold', verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.suptitle('Real Data Analysis: ABO Blood Type Locus (1000 Genomes)', 
                 fontsize=14, fontweight='bold', y=0.995)
    
    return fig


if __name__ == "__main__":
    posteriors = analyze_real_abo_data()
    
    print("\n" + "="*70)
    print("GENERATING EXPANDED DATASET ANALYSIS")
    print("="*70)
    
    from visualization_expanded import plot_real_abo_analysis_expanded
    
    expanded_fig, expanded_posteriors = plot_real_abo_analysis_expanded(
        save_path=os.path.join(project_root, 'results', 'real_abo_analysis_expanded.png')
    )
