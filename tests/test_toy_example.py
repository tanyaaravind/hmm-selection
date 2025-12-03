"""
Toy Example: Test Forward-Backward Algorithm on Small Dataset

This validates that our HMM implementation works correctly
before applying it to real data.

Expected behavior:
- SNPs with FST ~0.05 should have high P(neutral)
- SNPs with FST ~0.25 should have high P(selection)
"""

import numpy as np
import sys
import os

# Get the directory where THIS file is located
current_dir = os.path.dirname(os.path.abspath(__file__))
# Go up one level to project root, then into src
project_root = os.path.dirname(current_dir)
src_path = os.path.join(project_root, 'src')
sys.path.append(src_path)


def toy_example():
    """
    Create a tiny 5-SNP example with known pattern:
    - Positions 1-2: Low FST (neutral)
    - Positions 3-4: High FST (selection)
    - Position 5: Low FST (neutral)
    """
    print("="*60)
    print("TOY EXAMPLE: 5 SNPs with Known Pattern")
    print("="*60)
    
    # ===== CREATE TOY DATA =====
    observations = np.array([0.05, 0.06, 0.24, 0.26, 0.05])
    positions = np.array([0, 1000, 2000, 3000, 4000])
    
    print("\nInput Data:")
    print("-" * 40)
    print("Position (bp) | FST Value | Expected State")
    print("-" * 40)
    expected = ['Neutral', 'Neutral', 'Selection', 'Selection', 'Neutral']
    for i in range(len(observations)):
        print(f"{positions[i]:13d} | {observations[i]:9.3f} | {expected[i]}")
    
    # ===== SET UP HMM PARAMETERS =====
    emission_params = {
        'neutral': {'mean': 0.05, 'std': 0.02},
        'selection': {'mean': 0.25, 'std': 0.03}
    }
    
    transition_params = {
        'distance_scale': 1000  # 1kb distance scale
    }
    
    print("\n" + "="*60)
    print("RUNNING HMM")
    print("="*60)
    
    # ===== CREATE AND RUN HMM =====
    hmm = SelectionHMM(emission_params, transition_params)
    posteriors = hmm.posterior_probabilities(observations, positions)
    
    # ===== DISPLAY RESULTS =====
    print("\n" + "="*60)
    print("RESULTS: Posterior Probabilities")
    print("="*60)
    print("\nPosition (bp) | FST   | P(Neutral) | P(Selection) | Predicted State")
    print("-" * 75)
    
    for i in range(len(observations)):
        pred_state = 'Selection' if posteriors[i, 1] > 0.5 else 'Neutral'
        correct = '✓' if pred_state == expected[i] else '✗'
        
        print(f"{positions[i]:13d} | {observations[i]:.3f} | "
              f"{posteriors[i, 0]:10.3f} | {posteriors[i, 1]:12.3f} | "
              f"{pred_state:9s} {correct}")
    
    # ===== ACCURACY CHECK =====
    predicted = (posteriors[:, 1] > 0.5).astype(int)
    true_labels = np.array([0, 0, 1, 1, 0])  # 0=neutral, 1=selection
    accuracy = (predicted == true_labels).mean()
    
    print("\n" + "="*60)
    print(f"ACCURACY: {accuracy:.1%}")
    print("="*60)
    
    if accuracy >= 0.8:
        print("✓ SUCCESS: HMM correctly identifies states!")
    else:
        print("⚠ WARNING: HMM may need parameter tuning")
    
    return posteriors


if __name__ == "__main__":
    posteriors = toy_example()