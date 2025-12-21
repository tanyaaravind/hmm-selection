
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.insert(0, current_dir)
results_path = os.path.join(project_root, 'results')

from hmm_core import SelectionHMM
from data_prep import load_hmm_data

sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['font.size'] = 10


def load_expanded_data(delta_af_path, allele_freq_path=None):
    
    delta_df = pd.read_csv(delta_af_path)
    delta_df = delta_df.sort_values('pos').reset_index(drop=True)
    
    positions = delta_df['pos'].values
    delta_af = delta_df['delta_af'].values
    af_yri = delta_df['af_YRI'].values
    af_ceu = delta_df['af_CEU'].values
    
    af_chb = None
    if allele_freq_path and os.path.exists(allele_freq_path):
        af_df = pd.read_csv(allele_freq_path)
        chb_df = af_df[af_df['pop'] == 'CHB'][['pos', 'af']].rename(columns={'af': 'af_CHB'})
        merged = delta_df[['pos']].merge(chb_df, on='pos', how='left')
        af_chb = merged['af_CHB'].values
    
    metadata = {
        'n_snps': len(positions),
        'region_start': int(positions.min()),
        'region_end': int(positions.max()),
        'region_size': int(positions.max() - positions.min()),
        'n_zeros': int(np.sum(delta_af == 0)),
        'n_nonzero': int(np.sum(delta_af > 0)),
        'n_high': int(np.sum(delta_af > 0.2))
    }
    
    return {
        'positions': positions,
        'delta_af': delta_af,
        'af_yri': af_yri,
        'af_ceu': af_ceu,
        'af_chb': af_chb,
        'metadata': metadata
    }


def create_expanded_figure(data, posteriors=None, save_path='results/abo_expanded_visualization.png'):
    
    positions = data['positions']
    delta_af = data['delta_af']
    af_yri = data['af_yri']
    af_ceu = data['af_ceu']
    af_chb = data['af_chb']
    metadata = data['metadata']
    
    pos_kb = (positions - positions[0]) / 1000 
    
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1.2], hspace=0.35)
    
    ax1 = fig.add_subplot(gs[0])
    
    marker_size = 4 if len(positions) > 500 else 8
    line_width = 1.5 if len(positions) > 500 else 2
    
    ax1.plot(pos_kb, af_yri, 'o-', color='#E64B35', linewidth=line_width, 
             markersize=marker_size, label='YRI (African)', alpha=0.8)
    ax1.plot(pos_kb, af_ceu, 's-', color='#4DBBD5', linewidth=line_width, 
             markersize=marker_size, label='CEU (European)', alpha=0.8)
    
    if af_chb is not None:
        chb_mask = ~np.isnan(af_chb)
        if np.any(chb_mask):
            ax1.plot(pos_kb[chb_mask], af_chb[chb_mask], '^-', 
                    color='#00A087', linewidth=line_width, 
                    markersize=marker_size, label='CHB (East Asian)', alpha=0.8)
    
    ax1.set_ylabel('Allele Frequency', fontsize=11, fontweight='bold')
    ax1.set_title('Panel A: Allele Frequencies Across Three Populations', 
                  fontsize=12, fontweight='bold')
    ax1.legend(loc='best', fontsize=10, framealpha=0.9)
    ax1.grid(alpha=0.3)
    ax1.set_ylim(-0.05, 1.05)
    
    ax2 = fig.add_subplot(gs[1])
    
    ax2.plot(pos_kb, delta_af, 'o-', color='purple', linewidth=2.5, 
             markersize=marker_size, label='ΔAF = |AF_YRI - AF_CEU|', alpha=0.8)
    
    high_mask = delta_af > 0.2
    if np.any(high_mask):
        ax2.plot(pos_kb[high_mask], delta_af[high_mask], 'o', color='red', 
                markersize=12, alpha=0.7, label=f'High differentiation (ΔAF > 0.2, {np.sum(high_mask)} SNPs)',
                zorder=10)
    
    moderate_mask = (delta_af > 0.1) & (delta_af <= 0.2)
    if np.any(moderate_mask):
        ax2.plot(pos_kb[moderate_mask], delta_af[moderate_mask], 's', color='orange', 
                markersize=8, alpha=0.6, label=f'Moderate differentiation (0.1 < ΔAF ≤ 0.2, {np.sum(moderate_mask)} SNPs)',
                zorder=9)
    
    ax2.axhline(0.04, color='green', linestyle='--', linewidth=1.5, 
                alpha=0.6, label='Neutral-like (<0.06)')
    ax2.axhline(0.10, color='red', linestyle='--', linewidth=1.5, 
                alpha=0.6, label='Selection-like (>0.10)')
    
    ax2.set_ylabel('ΔAF (Population Differentiation)', fontsize=11, fontweight='bold')
    ax2.set_title('Panel B: Population Differentiation Signal', 
                  fontsize=12, fontweight='bold')
    ax2.legend(loc='best', fontsize=9, framealpha=0.9)
    ax2.grid(alpha=0.3)
    ax2.set_ylim(0, max(0.7, delta_af.max() * 1.1))
    
    ax3 = fig.add_subplot(gs[2])
    
    if posteriors is not None:
        ax3.plot(pos_kb, posteriors[:, 1], '-', color='red', linewidth=3, 
                label='P(Selection | Data)')
        ax3.fill_between(pos_kb, 0, posteriors[:, 1], alpha=0.3, color='red')
        
        ax3.plot(pos_kb, posteriors[:, 1], 'o', color='darkred', 
                markersize=marker_size, alpha=0.7)
        
        ax3.axhline(0.5, color='gray', linestyle='--', linewidth=1.5, 
                   alpha=0.7, label='Decision threshold (0.5)')
        
        ax3.set_ylabel('Posterior Probability', fontsize=11, fontweight='bold')
        ax3.set_title('Panel C: HMM Detection of Selection Signal', 
                      fontsize=12, fontweight='bold')
        
        n_selection = np.sum(posteriors[:, 1] > 0.5)
        summary_text = f'SNPs identified as under selection: {n_selection}/{len(positions)}'
        ax3.text(0.02, 0.95, summary_text, transform=ax3.transAxes,
                fontsize=11, fontweight='bold', verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        ax3.set_ylim(-0.05, 1.05)
    else:
        ax3.axis('off')
        summary_text = f"""
Dataset Summary:
• Total SNPs: {metadata['n_snps']:,}
• Genomic region: {metadata['region_start']:,} - {metadata['region_end']:,} bp
• Region size: {metadata['region_size']:,} bp (~{metadata['region_size']/1000:.1f} kb)

DeltaAF Distribution:
• Zero differentiation: {metadata['n_zeros']:,} SNPs ({metadata['n_zeros']/metadata['n_snps']*100:.1f}%)
• Non-zero differentiation: {metadata['n_nonzero']:,} SNPs ({metadata['n_nonzero']/metadata['n_snps']*100:.1f}%)
• High differentiation (ΔAF > 0.2): {metadata['n_high']:,} SNPs ({metadata['n_high']/metadata['n_snps']*100:.1f}%)

Next Steps:
• Run HMM analysis to detect selection signals
• Panel C will show posterior probabilities when HMM results are available
        """
        ax3.text(0.05, 0.95, summary_text, transform=ax3.transAxes,
                fontsize=11, verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        ax3.set_title('Panel C: HMM Results (Pending)', 
                      fontsize=12, fontweight='bold')
    
    ax3.set_xlabel('Position Relative to First SNP (kb)', fontsize=11, fontweight='bold')
    if posteriors is not None:
        ax3.legend(loc='best', fontsize=10, framealpha=0.9)
    ax3.grid(alpha=0.3)
    
    plt.suptitle('Real Data Analysis: ABO Blood Type Locus (1000 Genomes)', 
                 fontsize=14, fontweight='bold', y=0.995)
    
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"✅ Saved figure to: {save_path}")
    
    return fig


def plot_real_abo_analysis_expanded(save_path=None):
    
    if save_path is None:
        save_path = os.path.join(results_path, 'real_abo_analysis_expanded.png')
    
    print("="*70)
    print("REAL ABO ANALYSIS: Expanded Data with HMM")
    print("="*70)
    
    print("\n1. Loading data for HMM analysis...")
    delta_af_path = os.path.join(project_root, 'results', 'abo_expanded', 'delta_af_nonzero.csv')
    if not os.path.exists(delta_af_path):
        delta_af_path = os.path.join(project_root, 'results', 'abo_expanded', 'delta_af.csv')
        print("   (Using full dataset - consider filtering zeros for cleaner analysis)")
    else:
        print("   (Using filtered non-zero dataset)")
    allele_freq_path = os.path.join(project_root, 'results', 'abo_expanded', 'allele_frequencies.csv')
    
    if not os.path.exists(delta_af_path):
        print(f"❌ Error: Could not find data at {delta_af_path}")
        print("   Make sure you've run data_prep.py first!")
        return None
    
    observations, positions, metadata = load_hmm_data(
        delta_af_path,
        pop_a='YRI',
        pop_b='CEU'
    )
    
    print(f"✅ Loaded {metadata['n_snps']} SNPs")
    print(f"   Region: {metadata['region_start']:,} - {metadata['region_end']:,} bp")
    print(f"   Region size: {metadata['region_size']:,} bp (~{metadata['region_size']/1000:.1f} kb)")
    print(f"   DeltaAF range: {observations.min():.3f} - {observations.max():.3f}")
    
    print("\n2. Setting up HMM...")
    
    emission_params = {
        'neutral': {'mean': 0.04, 'std': 0.03}, 
        'selection': {'mean': 0.12, 'std': 0.03}
    }
    
    transition_params = {
        'distance_scale': 2000
    }
    
    hmm = SelectionHMM(emission_params, transition_params)
    
    print("2b. Running Baum–Welch (EM) to refine parameters...")
    hmm.fit(observations, positions, n_iter=10, verbose=False)

    print("\n3. Running Forward-Backward algorithm...")
    print("   (This may take a moment...)")
    
    posteriors = hmm.posterior_probabilities(observations, positions)
    
    print("\n4. Analyzing results...")
    
    predicted_states = (posteriors[:, 1] > 0.5).astype(int)
    n_selection = np.sum(predicted_states)
    n_neutral = len(predicted_states) - n_selection
    
    print(f"\n   Results Summary:")
    print(f"   • SNPs predicted as Neutral: {n_neutral} ({n_neutral/len(predicted_states)*100:.1f}%)")
    print(f"   • SNPs predicted as Selection: {n_selection} ({n_selection/len(predicted_states)*100:.1f}%)")
    
    high_selection_mask = posteriors[:, 1] > 0.8
    if np.any(high_selection_mask):
        high_selection_positions = positions[high_selection_mask]
        print(f"\n   High-confidence selection regions (P > 0.8):")
        print(f"   • {np.sum(high_selection_mask)} SNPs")
        print(f"   • Position range: {high_selection_positions.min():,} - {high_selection_positions.max():,} bp")
    
    print("\n5. Loading full data for visualization...")
    
    viz_data = load_expanded_data(delta_af_path, allele_freq_path)
    
    print("\n6. Creating comprehensive figure...")
    
    fig = create_expanded_figure(viz_data, posteriors=posteriors, save_path=save_path)
    
    print(f"\n✅ Figure saved to: {save_path}")
    
    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)
    print(f"Total SNPs: {len(observations)}")
    print(f"Neutral SNPs: {n_neutral}")
    print(f"Selection SNPs: {n_selection}")
    print(f"\nDeltaAF Statistics:")
    print(f"  Mean: {observations.mean():.3f}")
    print(f"  Median: {np.median(observations):.3f}")
    print(f"  Std: {observations.std():.3f}")
    print(f"  Range: {observations.min():.3f} - {observations.max():.3f}")
    print(f"\nHMM Results:")
    print(f"  SNPs with P(Selection) > 0.5: {n_selection} ({n_selection/len(observations)*100:.1f}%)")
    print(f"  SNPs with P(Selection) > 0.8: {np.sum(posteriors[:, 1] > 0.8)} ({np.sum(posteriors[:, 1] > 0.8)/len(observations)*100:.1f}%)")
    print("="*70)
    
    return fig, posteriors


def plot_preliminary_results_expanded(save_path=None):

    if save_path is None:
        save_path = os.path.join(results_path, 'preliminary_results_expanded.png')
    
    print("="*60)
    print("GENERATING PRELIMINARY RESULTS FIGURE (EXPANDED DATA)")
    print("="*60)
    
    print("\n1. Loading expanded ABO data...")
    delta_af_path = os.path.join(project_root, 'results', 'abo_expanded', 'delta_af_nonzero.csv')
    if not os.path.exists(delta_af_path):
        delta_af_path = os.path.join(project_root, 'results', 'abo_expanded', 'delta_af.csv')
    
    if not os.path.exists(delta_af_path):
        print(f"❌ Error: Could not find data at {delta_af_path}")
        print("   Make sure you've run data_prep.py first!")
        return None
    
    observations, positions, metadata = load_hmm_data(delta_af_path, pop_a='YRI', pop_b='CEU')
    
    print(f"   Loaded {metadata['n_snps']} SNPs")
    print(f"   Region: {metadata['region_start']:,} - {metadata['region_end']:,} bp")
    print(f"   DeltaAF range: {observations.min():.3f} - {observations.max():.3f}")
    
    print("2. Initializing HMM...")
    emission_params = {
        'neutral': {'mean': 0.05, 'std': 0.02},
        'selection': {'mean': 0.25, 'std': 0.05}
    }
    
    transition_params = {
        'distance_scale': 2000
    }
    
    hmm = SelectionHMM(emission_params, transition_params)
    
    print("2b. Running Baum–Welch (EM) to refine parameters...")
    hmm.fit(observations, positions, n_iter=10, verbose=False)

    print("3. Running Forward-Backward algorithm...")
    print("   (This may take a moment...)")
    posteriors = hmm.posterior_probabilities(observations, positions)
    
    predicted_states = (posteriors[:, 1] > 0.5).astype(int)
    n_selection = np.sum(predicted_states)
    n_neutral = len(predicted_states) - n_selection
    
    print(f"4. Results:")
    print(f"   SNPs predicted as Neutral: {n_neutral} ({n_neutral/len(predicted_states)*100:.1f}%)")
    print(f"   SNPs predicted as Selection: {n_selection} ({n_selection/len(predicted_states)*100:.1f}%)")
    
    print("5. Creating figure...")
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    
    positions_kb = (positions - positions[0]) / 1000
    
    ax1 = axes[0]
    
    ax1.plot(positions_kb, observations, 'o-', color='steelblue', 
             alpha=0.7, linewidth=1.5, markersize=6, label='Observed ΔAF')
    
    ax1.axhline(0.05, color='green', linestyle='--', linewidth=2, 
                alpha=0.7, label='Neutral mean (0.05)')
    ax1.axhline(0.25, color='red', linestyle='--', linewidth=2, 
                alpha=0.7, label='Selection mean (0.25)')
    
    ax1.set_ylabel('ΔAF (Population Differentiation)', fontsize=11, fontweight='bold')
    ax1.set_title('Panel A: Observed ΔAF Values Along Chromosome', 
                  fontsize=12, fontweight='bold')
    ax1.legend(loc='upper left', fontsize=9, framealpha=0.9)
    ax1.grid(alpha=0.3)
    ax1.set_ylim(-0.05, 0.4)
    
    ax2 = axes[1]
    
    ax2.plot(positions_kb, posteriors[:, 1], '-', color='red', 
             linewidth=3, label='P(Selection | Data)')
    ax2.fill_between(positions_kb, 0, posteriors[:, 1], 
                     alpha=0.3, color='red')
    
    ax2.axhline(0.5, color='gray', linestyle='--', linewidth=1.5, 
                alpha=0.7, label='Decision threshold')
    
    ax2.set_xlabel('Genomic Position (kb)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Posterior Probability', fontsize=11, fontweight='bold')
    ax2.set_title('Panel B: HMM Posterior Probability of Selection State', 
                  fontsize=12, fontweight='bold')
    ax2.legend(loc='upper left', fontsize=9, framealpha=0.9)
    ax2.grid(alpha=0.3)
    ax2.set_ylim(-0.05, 1.05)
    
    summary_text = f'SNPs under selection: {n_selection}/{len(observations)} ({n_selection/len(observations)*100:.1f}%)'
    ax2.text(0.98, 0.05, summary_text, 
             transform=ax2.transAxes,
             fontsize=12, fontweight='bold',
             verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"6. Figure saved to: {save_path}")
    
    plt.close()
    
    return fig, posteriors


if __name__ == "__main__":
    fig, posteriors = plot_preliminary_results_expanded()
    
    
