"""
Simple Simulation: Generate Fake Genomic Data with Known Selection Signal

This creates realistic-looking FST data where:
- First half: Neutral region (low FST)
- Second half: Selection region (high FST)

We'll use this to validate that our HMM can detect selection
when we KNOW where it is.
"""

import numpy as np


def simulate_selection_region(n_snps=50, neutral_fst_mean=0.05, 
                              selection_fst_mean=0.25, seed=42):
    """
    Generate simulated FST data with selection in second half
    
    Parameters:
    -----------
    n_snps : int
        Total number of SNPs to simulate
    neutral_fst_mean : float
        Mean FST for neutral regions
    selection_fst_mean : float
        Mean FST for selection regions
    seed : int
        Random seed for reproducibility
    
    Returns:
    --------
    observations : numpy array, shape (n_snps,)
        FST values
    positions : numpy array, shape (n_snps,)
        Genomic positions in bp
    true_states : numpy array, shape (n_snps,)
        Ground truth labels (0=neutral, 1=selection)
    """
    np.random.seed(seed)
    
    # Split into two halves
    n_neutral = n_snps // 2
    n_selection = n_snps - n_neutral
    
    # Generate neutral region FST values
    # FST ~ N(0.05, 0.02) - small differences between populations
    neutral_fst = np.random.normal(neutral_fst_mean, 0.02, n_neutral)
    
    # Generate selection region FST values  
    # FST ~ N(0.25, 0.05) - large differences between populations
    selection_fst = np.random.normal(selection_fst_mean, 0.05, n_selection)
    
    # Clip to valid FST range [0, 1]
    neutral_fst = np.clip(neutral_fst, 0, 1)
    selection_fst = np.clip(selection_fst, 0, 1)
    
    # Combine
    observations = np.concatenate([neutral_fst, selection_fst])
    
    # Generate genomic positions
    # Assume SNPs span 100kb region (100,000 bp)
    # Evenly spaced for simplicity
    positions = np.linspace(0, 100000, n_snps, dtype=int)
    
    # Ground truth labels
    true_states = np.concatenate([
        np.zeros(n_neutral, dtype=int),      # 0 = neutral
        np.ones(n_selection, dtype=int)       # 1 = selection
    ])
    
    print(f"Simulated {n_snps} SNPs across 100kb region:")
    print(f"  Neutral region: positions 0-{positions[n_neutral-1]} bp ({n_neutral} SNPs)")
    print(f"    Mean FST = {neutral_fst.mean():.3f} (expected {neutral_fst_mean})")
    print(f"  Selection region: positions {positions[n_neutral]}-100000 bp ({n_selection} SNPs)")
    print(f"    Mean FST = {selection_fst.mean():.3f} (expected {selection_fst_mean})")
    
    return observations, positions, true_states


def simulate_multiple_regions(n_regions=3, snps_per_region=100, 
                               selection_probability=0.3, seed=42):
    """
    Generate more complex simulation with multiple neutral/selection regions
    
    This is for later use - creates a more realistic genome with
    scattered selection regions mixed with neutral regions.
    
    Parameters:
    -----------
    n_regions : int
        Number of genomic regions to simulate
    snps_per_region : int
        SNPs per region
    selection_probability : float
        Probability that a region is under selection
    seed : int
        Random seed
    
    Returns:
    --------
    observations, positions, true_states (same as above)
    """
    np.random.seed(seed)
    
    all_observations = []
    all_positions = []
    all_states = []
    
    current_position = 0
    
    for region_idx in range(n_regions):
        # Randomly decide if this region has selection
        has_selection = np.random.random() < selection_probability
        
        if has_selection:
            # Selection region: high FST
            fst_values = np.random.normal(0.25, 0.05, snps_per_region)
            fst_values = np.clip(fst_values, 0, 1)
            states = np.ones(snps_per_region, dtype=int)
        else:
            # Neutral region: low FST
            fst_values = np.random.normal(0.05, 0.02, snps_per_region)
            fst_values = np.clip(fst_values, 0, 1)
            states = np.zeros(snps_per_region, dtype=int)
        
        # Positions for this region (50kb per region)
        positions = np.linspace(current_position, 
                               current_position + 50000, 
                               snps_per_region, dtype=int)
        
        all_observations.append(fst_values)
        all_positions.append(positions)
        all_states.append(states)
        
        current_position += 50000
    
    observations = np.concatenate(all_observations)
    positions = np.concatenate(all_positions)
    true_states = np.concatenate(all_states)
    
    n_selection_regions = sum(has_selection for region_idx in range(n_regions))
    print(f"Simulated {n_regions} regions ({len(observations)} total SNPs):")
    print(f"  {n_selection_regions} selection regions")
    print(f"  {n_regions - n_selection_regions} neutral regions")
    
    return observations, positions, true_states


if __name__ == "__main__":
    print("="*60)
    print("SIMPLE SIMULATION TEST")
    print("="*60)
    
    # Test simple simulation
    obs, pos, states = simulate_selection_region(n_snps=50)
    
    print(f"\nFirst 5 SNPs (neutral):")
    for i in range(5):
        print(f"  Position {pos[i]:6d} bp: FST = {obs[i]:.3f}, True state = {states[i]}")
    
    print(f"\nLast 5 SNPs (selection):")
    for i in range(-5, 0):
        print(f"  Position {pos[i]:6d} bp: FST = {obs[i]:.3f}, True state = {states[i]}")