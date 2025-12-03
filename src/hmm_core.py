"""
Core Hidden Markov Model for Selection Detection
Using FST values as observations to detect selection vs neutral evolution

Author: [Your Team]
Date: December 2024
"""

import numpy as np
from scipy.stats import norm


class SelectionHMM:
    """
    2-State Hidden Markov Model for detecting natural selection
    
    States:
        0 = Neutral: Random genetic drift only
        1 = Selection: Adaptive evolution
    
    Observations:
        FST values (population differentiation statistic)
    
    Key Features:
        - Distance-dependent transition probabilities
        - Log-space computation for numerical stability
        - Emission model: Normal distributions for FST
    """
    
    def __init__(self, emission_params, transition_params):
        """
        Initialize HMM with parameters
        
        Parameters:
        -----------
        emission_params : dict
            Parameters for emission distributions
            Format: {
                'neutral': {'mean': float, 'std': float},
                'selection': {'mean': float, 'std': float}
            }
            Example: 
                neutral state expects low FST (mean=0.05)
                selection state expects high FST (mean=0.25)
        
        transition_params : dict
            Parameters for distance-dependent transitions
            Format: {
                'distance_scale': float  # Distance in bp where correlation decays
            }
            Example: distance_scale=1000 means SNPs 1kb apart are ~independent
        """
        self.n_states = 2
        self.state_names = ['neutral', 'selection']
        
        # Store emission parameters
        self.emission_params = emission_params
        
        # Store transition parameters
        self.transition_params = transition_params
        self.distance_scale = transition_params['distance_scale']
        
        # Initial state probabilities (equal for now)
        self.log_initial = np.log([0.5, 0.5])
        
        print("SelectionHMM initialized")
        print(f"  Neutral state: FST ~ N({emission_params['neutral']['mean']}, "
              f"{emission_params['neutral']['std']})")
        print(f"  Selection state: FST ~ N({emission_params['selection']['mean']}, "
              f"{emission_params['selection']['std']})")
        print(f"  Transition distance scale: {self.distance_scale} bp")
    
    
    def emission_log_prob(self, observation, state):
        """
        Calculate log P(FST | state)
        
        Uses normal distribution - biologically, FST values tend to be
        normally distributed around different means for neutral vs selection
        
        Parameters:
        -----------
        observation : float
            FST value (must be between 0 and 1)
        state : int
            Hidden state (0=neutral, 1=selection)
        
        Returns:
        --------
        log_prob : float
            Log probability of observing this FST in this state
        
        Mathematical Formula:
        --------------------
        If FST ~ N(μ, σ²), then:
        log P(FST) = -log(σ√(2π)) - (FST - μ)² / (2σ²)
        """
        state_name = self.state_names[state]
        params = self.emission_params[state_name]
        mean = params['mean']
        std = params['std']
        
        # Use scipy's norm for numerical stability
        log_prob = norm.logpdf(observation, loc=mean, scale=std)
        
        return log_prob
    
    
    def transition_log_prob(self, from_state, to_state, distance):
        """
        Calculate log P(to_state | from_state, distance)
        
        BIOLOGICAL INTUITION:
        ---------------------
        Selection doesn't affect individual SNPs - it affects REGIONS.
        If position t is under selection, nearby positions are also likely
        to be under selection. The probability of staying in the same state
        decreases exponentially with distance.
        
        Distance-dependent model:
        - Short distance (100 bp): Very likely to stay in same state
        - Medium distance (1 kb): Moderate correlation  
        - Long distance (10+ kb): States become independent
        
        Parameters:
        -----------
        from_state : int
            Previous state (0 or 1)
        to_state : int
            Current state (0 or 1)
        distance : float
            Genomic distance in base pairs between positions
        
        Returns:
        --------
        log_prob : float
            Log transition probability
        
        Mathematical Formula:
        --------------------
        P(stay in state) = 0.5 + 0.5 * exp(-distance / λ)
        where λ = distance_scale parameter
        
        At distance = 0: P(stay) = 1.0 (certain)
        At distance = λ: P(stay) = 0.68 (moderately correlated)
        At distance >> λ: P(stay) = 0.5 (independent)
        """
        # Base probability of staying in same state
        # Decays exponentially with distance
        p_stay = 0.5 + 0.5 * np.exp(-distance / self.distance_scale)
        
        if from_state == to_state:
            # Staying in same state
            return np.log(p_stay)
        else:
            # Switching to different state
            p_switch = 1.0 - p_stay
            return np.log(p_switch)
    
    
    def forward_algorithm(self, observations, positions):
        """
        Forward algorithm: Computes P(observations up to t, state at t)
        
        WHAT IT DOES:
        -------------
        For each position t and each state k, computes:
        α_t(k) = P(x_1, x_2, ..., x_t, state_t = k)
        
        This is the probability of:
        1. Seeing all observations from position 1 to t, AND
        2. Being in state k at position t
        
        WHY WE NEED IT:
        ---------------
        Forward probabilities are half of what we need to compute posteriors.
        Combined with backward probabilities, they give us:
        P(state_t = k | all observations)
        
        ALGORITHM:
        ----------
        1. Initialization: α_1(k) = P(state_1 = k) * P(x_1 | state_1 = k)
        2. Recursion: α_t(k) = P(x_t | state_t = k) * Σ_j [α_{t-1}(j) * P(state_t=k | state_{t-1}=j)]
        3. All computations in log-space to avoid underflow
        
        Parameters:
        -----------
        observations : numpy array, shape (T,)
            FST values at each SNP position
        positions : numpy array, shape (T,)
            Genomic positions (in bp) of each SNP
        
        Returns:
        --------
        log_forward : numpy array, shape (T, n_states)
            Log forward probabilities
            log_forward[t, k] = log α_t(k)
        """
        T = len(observations)
        log_forward = np.zeros((T, self.n_states))
        
        # ===== INITIALIZATION (t=0) =====
        # At first position, probability = P(initial state) * P(observation | state)
        for state in range(self.n_states):
            log_forward[0, state] = (
                self.log_initial[state] +
                self.emission_log_prob(observations[0], state)
            )
        
        # ===== RECURSION (t=1 to T-1) =====
        for t in range(1, T):
            # Calculate distance from previous SNP
            distance = positions[t] - positions[t-1]
            
            # For each possible current state
            for to_state in range(self.n_states):
                # Sum over all possible previous states
                # Using log-sum-exp trick for numerical stability
                log_probs = []
                
                for from_state in range(self.n_states):
                    # log P(current state | previous state, distance)
                    log_transition = self.transition_log_prob(
                        from_state, to_state, distance
                    )
                    
                    # log[α_{t-1}(j) * P(state_t | state_{t-1})]
                    log_prob = log_forward[t-1, from_state] + log_transition
                    log_probs.append(log_prob)
                
                # Log-sum-exp: log(Σ exp(x_i)) = max(x) + log(Σ exp(x_i - max(x)))
                # This prevents numerical underflow
                max_log_prob = max(log_probs)
                sum_exp = sum(np.exp(lp - max_log_prob) for lp in log_probs)
                log_sum = max_log_prob + np.log(sum_exp)
                
                # Add emission probability
                log_forward[t, to_state] = (
                    log_sum + 
                    self.emission_log_prob(observations[t], to_state)
                )
        
        return log_forward
    
    
    def backward_algorithm(self, observations, positions):
        """
        Backward algorithm: Computes P(future observations | state at t)
        
        WHAT IT DOES:
        -------------
        For each position t and each state k, computes:
        β_t(k) = P(x_{t+1}, x_{t+2}, ..., x_T | state_t = k)
        
        This is the probability of seeing all FUTURE observations
        given that we're in state k at position t
        
        WHY WE NEED IT:
        ---------------
        Combined with forward probabilities:
        P(state_t = k | all data) ∝ α_t(k) * β_t(k)
        
        ALGORITHM:
        ----------
        1. Initialization: β_T(k) = 1 for all k (log β_T(k) = 0)
        2. Recursion (going backwards): 
           β_t(k) = Σ_j [P(state_{t+1}=j | state_t=k) * P(x_{t+1} | state_{t+1}=j) * β_{t+1}(j)]
        
        Parameters:
        -----------
        observations : numpy array, shape (T,)
            FST values at each SNP position
        positions : numpy array, shape (T,)
            Genomic positions (in bp) of each SNP
        
        Returns:
        --------
        log_backward : numpy array, shape (T, n_states)
            Log backward probabilities
            log_backward[t, k] = log β_t(k)
        """
        T = len(observations)
        log_backward = np.zeros((T, self.n_states))
        
        # ===== INITIALIZATION (t=T-1) =====
        # At last position, β_T(k) = 1 for all states
        # log(1) = 0, so already initialized correctly
        
        # ===== RECURSION (t=T-2 down to 0) =====
        for t in range(T-2, -1, -1):
            # Calculate distance to next SNP
            distance = positions[t+1] - positions[t]
            
            # For each possible current state
            for from_state in range(self.n_states):
                # Sum over all possible next states
                log_probs = []
                
                for to_state in range(self.n_states):
                    # log P(next state | current state, distance)
                    log_transition = self.transition_log_prob(
                        from_state, to_state, distance
                    )
                    
                    # log P(next observation | next state)
                    log_emission = self.emission_log_prob(
                        observations[t+1], to_state
                    )
                    
                    # Combine: transition * emission * future probability
                    log_prob = (
                        log_transition + 
                        log_emission + 
                        log_backward[t+1, to_state]
                    )
                    log_probs.append(log_prob)
                
                # Log-sum-exp
                max_log_prob = max(log_probs)
                sum_exp = sum(np.exp(lp - max_log_prob) for lp in log_probs)
                log_backward[t, from_state] = max_log_prob + np.log(sum_exp)
        
        return log_backward
    
    
    def posterior_probabilities(self, observations, positions):
        """
        Compute posterior probabilities: P(state_t | all observations)
        
        THIS IS THE KEY OUTPUT FOR YOUR PRELIMINARY RESULTS!
        
        WHAT IT COMPUTES:
        -----------------
        For each SNP position t and each state k:
        P(state_t = k | x_1, x_2, ..., x_T)
        
        This tells us: "Given ALL the FST values we observed,
        what's the probability that position t is under selection?"
        
        BIOLOGICAL INTERPRETATION:
        --------------------------
        - If P(selection | data) > 0.9: Strong evidence for selection at this SNP
        - If P(selection | data) ≈ 0.5: Ambiguous, could be either
        - If P(selection | data) < 0.1: Strong evidence for neutral evolution
        
        MATHEMATICAL FORMULA:
        --------------------
        P(state_t = k | all data) = [α_t(k) * β_t(k)] / Σ_j [α_t(j) * β_t(j)]
        
        In log-space:
        log P(state_t = k | all data) = log α_t(k) + log β_t(k) - log[Σ_j exp(log α_t(j) + log β_t(j))]
        
        Parameters:
        -----------
        observations : numpy array, shape (T,)
            FST values
        positions : numpy array, shape (T,)
            Genomic positions
        
        Returns:
        --------
        posteriors : numpy array, shape (T, n_states)
            Posterior probabilities
            posteriors[t, k] = P(state_t = k | all observations)
            Each row sums to 1.0
        """
        # Run forward and backward algorithms
        print("Running forward algorithm...")
        log_forward = self.forward_algorithm(observations, positions)
        
        print("Running backward algorithm...")
        log_backward = self.backward_algorithm(observations, positions)
        
        # Compute posteriors
        print("Computing posterior probabilities...")
        T = len(observations)
        posteriors = np.zeros((T, self.n_states))
        
        for t in range(T):
            # Combine forward and backward probabilities
            log_post = log_forward[t, :] + log_backward[t, :]
            
            # Normalize to get probabilities (convert from log-space)
            # Use log-sum-exp for numerical stability
            max_log = np.max(log_post)
            exp_log = np.exp(log_post - max_log)
            posteriors[t, :] = exp_log / exp_log.sum()
        
        return posteriors


# ===== UTILITY FUNCTION FOR LOG-SUM-EXP =====
def log_sum_exp(log_values):
    """
    Numerically stable computation of log(sum(exp(log_values)))
    
    Used throughout HMM to avoid underflow when summing probabilities
    
    Standard approach would be:
        sum_prob = sum(exp(log_p) for log_p in log_values)
        result = log(sum_prob)
    
    Problem: exp(log_p) can underflow if log_p is very negative
    
    Solution:
        log(Σ exp(x_i)) = max(x) + log(Σ exp(x_i - max(x)))
    
    This shifts all values by subtracting the maximum before exp(),
    ensuring at least one value is exp(0) = 1, preventing underflow.
    """
    max_log = np.max(log_values)
    return max_log + np.log(np.sum(np.exp(log_values - max_log)))