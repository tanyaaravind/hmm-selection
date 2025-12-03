"""
Visualization: Generate Preliminary Results Figure

This creates the figure you'll show in your presentation!
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

sys.path.append('/home/claude/selection_hmm/src')
from hmm_core import SelectionHMM
from simulation import simulate_selection_region

# Set style for publication-quality figures
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['font.size'] = 10


def plot_preliminary_results(save_path='/home/claude/selection_hmm/results/preliminary_results.png'):
    """
    Generate the preliminary results figure for presentation
    
    Creates a 2-panel figure showing:
    1. Observed FST values along chromosome
    2. HMM posterior probability of selection state
    """
    print("="*60)
    print("GENERATING PRELIMINARY RESULTS FIGURE")
    print("="*60)
    
    # ===== GENERATE SIMULATED DATA =====
    print("\n1. Generating simulated data...")
    observations, positions, true_states = simulate_selection_region(n_snps=50, seed=42)
    
    # ===== SET UP HMM =====
    print("2. Initializing HMM...")
    emission_params = {
        'neutral': {'mean': 0.05, 'std': 0.02},
        'selection': {'mean': 0.25, 'std': 0.05}
    }
    
    transition_params = {
        'distance_scale': 2000  # 2kb correlation distance
    }
    
    hmm = SelectionHMM(emission_params, transition_params)
    
    # ===== RUN HMM =====
    print("3. Running Forward-Backward algorithm...")
    posteriors = hmm.posterior_probabilities(observations, positions)
    
    # ===== CALCULATE ACCURACY =====
    predicted_states = (posteriors[:, 1] > 0.5).astype(int)
    accuracy = (predicted_states == true_states).mean()
    print(f"4. Model accuracy: {accuracy:.1%}")
    
    # ===== CREATE FIGURE =====
    print("5. Creating figure...")
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    
    # Convert positions to kb for cleaner axis
    positions_kb = positions / 1000
    
    # ===== PANEL A: OBSERVED FST VALUES =====
    ax1 = axes[0]
    
    # Plot FST values
    ax1.plot(positions_kb, observations, 'o-', color='steelblue', 
             alpha=0.7, linewidth=1.5, markersize=6, label='Observed FST')
    
    # Add reference lines for expected FST values
    ax1.axhline(0.05, color='green', linestyle='--', linewidth=2, 
                alpha=0.7, label='Neutral mean (0.05)')
    ax1.axhline(0.25, color='red', linestyle='--', linewidth=2, 
                alpha=0.7, label='Selection mean (0.25)')
    
    # Mark true selection boundary
    boundary_position = positions_kb[len(positions)//2]
    ax1.axvline(boundary_position, color='black', linestyle=':', 
                linewidth=2, alpha=0.5, label='True selection boundary')
    
    # Shade regions
    ax1.axvspan(0, boundary_position, alpha=0.1, color='green', label='True neutral region')
    ax1.axvspan(boundary_position, positions_kb[-1], alpha=0.1, 
                color='red', label='True selection region')
    
    ax1.set_ylabel('FST (Population Differentiation)', fontsize=11, fontweight='bold')
    ax1.set_title('Panel A: Observed FST Values Along Chromosome', 
                  fontsize=12, fontweight='bold')
    ax1.legend(loc='upper left', fontsize=9, framealpha=0.9)
    ax1.grid(alpha=0.3)
    ax1.set_ylim(-0.05, 0.4)
    
    # ===== PANEL B: HMM POSTERIOR PROBABILITIES =====
    ax2 = axes[1]
    
    # Plot posterior probability of selection state
    ax2.plot(positions_kb, posteriors[:, 1], '-', color='red', 
             linewidth=3, label='P(Selection | Data)')
    ax2.fill_between(positions_kb, 0, posteriors[:, 1], 
                     alpha=0.3, color='red')
    
    # Add reference line at 0.5 (decision boundary)
    ax2.axhline(0.5, color='gray', linestyle='--', linewidth=1.5, 
                alpha=0.7, label='Decision threshold')
    
    # Mark true selection boundary
    ax2.axvline(boundary_position, color='black', linestyle=':', 
                linewidth=2, alpha=0.5, label='True boundary')
    
    ax2.set_xlabel('Genomic Position (kb)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Posterior Probability', fontsize=11, fontweight='bold')
    ax2.set_title('Panel B: HMM Posterior Probability of Selection State', 
                  fontsize=12, fontweight='bold')
    ax2.legend(loc='upper left', fontsize=9, framealpha=0.9)
    ax2.grid(alpha=0.3)
    ax2.set_ylim(-0.05, 1.05)
    
    # ===== ADD ACCURACY TEXT =====
    ax2.text(0.98, 0.05, f'Accuracy: {accuracy:.1%}', 
             transform=ax2.transAxes,
             fontsize=12, fontweight='bold',
             verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    # ===== SAVE FIGURE =====
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"6. Figure saved to: {save_path}")
    
    plt.show()
    
    # ===== PRINT SUMMARY STATISTICS =====
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(f"Total SNPs: {len(observations)}")
    print(f"Neutral SNPs: {sum(true_states == 0)}")
    print(f"Selection SNPs: {sum(true_states == 1)}")
    print(f"\nFST Statistics:")
    print(f"  Neutral region: mean={observations[true_states==0].mean():.3f}, "
          f"std={observations[true_states==0].std():.3f}")
    print(f"  Selection region: mean={observations[true_states==1].mean():.3f}, "
          f"std={observations[true_states==1].std():.3f}")
    print(f"\nHMM Performance:")
    print(f"  Accuracy: {accuracy:.1%}")
    print(f"  Correctly identified neutral: {sum((predicted_states == 0) & (true_states == 0))}/{sum(true_states == 0)}")
    print(f"  Correctly identified selection: {sum((predicted_states == 1) & (true_states == 1))}/{sum(true_states == 1)}")
    
    # Calculate sensitivity and specificity
    true_positives = sum((predicted_states == 1) & (true_states == 1))
    false_positives = sum((predicted_states == 1) & (true_states == 0))
    true_negatives = sum((predicted_states == 0) & (true_states == 0))
    false_negatives = sum((predicted_states == 0) & (true_states == 1))
    
    sensitivity = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    specificity = true_negatives / (true_negatives + false_positives) if (true_negatives + false_positives) > 0 else 0
    
    print(f"  Sensitivity (recall): {sensitivity:.1%}")
    print(f"  Specificity: {specificity:.1%}")
    print("="*60)
    
    return fig, posteriors


if __name__ == "__main__":
    fig, posteriors = plot_preliminary_results()