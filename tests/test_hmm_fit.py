"""
Minimal checks for training (Baum–Welch) and decoding on toy data.

This is a sanity test for Wellington's tasks:
- Confirm Viterbi path matches expected neutral/selection pattern.
- Confirm Baum–Welch nudges emission means in the right order.
"""

import numpy as np
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
src_path = os.path.join(project_root, "src")
sys.path.append(src_path)

from hmm_core import SelectionHMM  # noqa: E402


def test_fit_and_viterbi_toy():
    observations = np.array([0.05, 0.06, 0.24, 0.26, 0.05])
    positions = np.array([0, 1000, 2000, 3000, 4000])

    emission_params = {
        "neutral": {"mean": 0.05, "std": 0.02},
        "selection": {"mean": 0.25, "std": 0.03},
    }
    transition_params = {"distance_scale": 1000}

    hmm = SelectionHMM(emission_params, transition_params)

    # Viterbi should follow the obvious low-low-high-high-low pattern
    path = hmm.viterbi(observations, positions)
    assert path.tolist() == [0, 0, 1, 1, 0]

    # Fit should keep selection mean above neutral mean
    hmm.fit(observations, positions, n_iter=5, verbose=False)
    neutral_mean = hmm.emission_params["neutral"]["mean"]
    selection_mean = hmm.emission_params["selection"]["mean"]
    assert selection_mean > neutral_mean
    assert neutral_mean > 0  # sanity: no degenerate collapse
