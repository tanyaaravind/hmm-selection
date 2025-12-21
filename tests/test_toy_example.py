
import numpy as np
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
src_path = os.path.join(project_root, 'src')
sys.path.insert(0, src_path)

from hmm_core import SelectionHMM


def toy_example():

    print("="*60)
    print("TOY EXAMPLE: 5 SNPs with Known Pattern")
    print("="*60)
    observations = np.array([0.05, 0.06, 0.24, 0.26, 0.05])
    positions = np.array([0, 1000, 2000, 3000, 4000])
    
    print("\nInput Data:")
    print("-" * 40)
    print("Position (bp) | FST Value | Expected State")
    print("-" * 40)
    expected = ['Neutral', 'Neutral', 'Selection', 'Selection', 'Neutral']
    for i in range(len(observations)):
        print(f"{positions[i]:13d} | {observations[i]:9.3f} | {expected[i]}")
    
    emission_params = {
        'neutral': {'mean': 0.05, 'std': 0.02},
        'selection': {'mean': 0.25, 'std': 0.03}
    }
    
    transition_params = {
        'distance_scale': 1000
    }
    
    print("\n" + "="*60)
    print("RUNNING HMM")
    print("="*60)
    
    hmm = SelectionHMM(emission_params, transition_params)
    posteriors = hmm.posterior_probabilities(observations, positions)
    
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
    
    predicted = (posteriors[:, 1] > 0.5).astype(int)
    true_labels = np.array([0, 0, 1, 1, 0])
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