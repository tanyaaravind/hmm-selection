
import numpy as np
from scipy.stats import norm


class SelectionHMM:
    
    def __init__(self, emission_params, transition_params):
        self.n_states = 2
        self.state_names = ['neutral', 'selection']
        
        self.emission_params = emission_params
        
        self.transition_params = transition_params
        self.distance_scale = transition_params['distance_scale']
        
        self.log_initial = np.log([0.5, 0.5])
        
        print("SelectionHMM initialized")
        print(f"  Neutral state: FST ~ N({emission_params['neutral']['mean']}, "
              f"{emission_params['neutral']['std']})")
        print(f"  Selection state: FST ~ N({emission_params['selection']['mean']}, "
              f"{emission_params['selection']['std']})")
        print(f"  Transition distance scale: {self.distance_scale} bp")
    
    
    def emission_log_prob(self, observation, state):
        
        state_name = self.state_names[state]
        params = self.emission_params[state_name]
        mean = params['mean']
        std = params['std']
        
        log_prob = norm.logpdf(observation, loc=mean, scale=std)
        
        return log_prob
    
    
    def transition_log_prob(self, from_state, to_state, distance):
    
        p_stay = 0.5 + 0.5 * np.exp(-distance / self.distance_scale)
        
        if from_state == to_state:
            return np.log(p_stay)
        else:
            p_switch = 1.0 - p_stay
            return np.log(p_switch)
    
    
    def forward_algorithm(self, observations, positions):
        
        T = len(observations)
        log_forward = np.zeros((T, self.n_states))
        
        for state in range(self.n_states):
            log_forward[0, state] = (
                self.log_initial[state] +
                self.emission_log_prob(observations[0], state)
            )
        
        for t in range(1, T):
            distance = positions[t] - positions[t-1]
            
            for to_state in range(self.n_states):
                log_probs = []
                
                for from_state in range(self.n_states):
                    log_transition = self.transition_log_prob(
                        from_state, to_state, distance
                    )
                    
                    log_prob = log_forward[t-1, from_state] + log_transition
                    log_probs.append(log_prob)
                
                max_log_prob = max(log_probs)
                sum_exp = sum(np.exp(lp - max_log_prob) for lp in log_probs)
                log_sum = max_log_prob + np.log(sum_exp)
                
                log_forward[t, to_state] = (
                    log_sum + 
                    self.emission_log_prob(observations[t], to_state)
                )
        
        return log_forward
    
    
    def backward_algorithm(self, observations, positions):
        
        T = len(observations)
        log_backward = np.zeros((T, self.n_states))
        
        for t in range(T-2, -1, -1):
            distance = positions[t+1] - positions[t]
            
            for from_state in range(self.n_states):
                log_probs = []
                
                for to_state in range(self.n_states):

                    log_transition = self.transition_log_prob(
                        from_state, to_state, distance
                    )
                    
                    log_emission = self.emission_log_prob(
                        observations[t+1], to_state
                    )
                    
                    log_prob = (
                        log_transition + 
                        log_emission + 
                        log_backward[t+1, to_state]
                    )
                    log_probs.append(log_prob)
                
                max_log_prob = max(log_probs)
                sum_exp = sum(np.exp(lp - max_log_prob) for lp in log_probs)
                log_backward[t, from_state] = max_log_prob + np.log(sum_exp)
        
        return log_backward
    
    
    def posterior_probabilities(self, observations, positions):
        
        print("Running forward algorithm...")
        log_forward = self.forward_algorithm(observations, positions)
        
        print("Running backward algorithm...")
        log_backward = self.backward_algorithm(observations, positions)
        
        print("Computing posterior probabilities...")
        T = len(observations)
        posteriors = np.zeros((T, self.n_states))
        
        for t in range(T):
            log_post = log_forward[t, :] + log_backward[t, :]
            
            max_log = np.max(log_post)
            exp_log = np.exp(log_post - max_log)
            posteriors[t, :] = exp_log / exp_log.sum()
        
        return posteriors

    def viterbi(self, observations, positions):

        T = len(observations)
        log_delta = np.zeros((T, self.n_states))
        psi = np.zeros((T, self.n_states), dtype=int)

        for state in range(self.n_states):
            log_delta[0, state] = (
                self.log_initial[state] +
                self.emission_log_prob(observations[0], state)
            )

        for t in range(1, T):
            distance = positions[t] - positions[t-1]
            for to_state in range(self.n_states):
                candidates = []
                for from_state in range(self.n_states):
                    log_transition = self.transition_log_prob(
                        from_state, to_state, distance
                    )
                    candidates.append(log_delta[t-1, from_state] + log_transition)
                best_prev = np.argmax(candidates)
                log_delta[t, to_state] = candidates[best_prev] + self.emission_log_prob(
                    observations[t], to_state
                )
                psi[t, to_state] = best_prev

        path = np.zeros(T, dtype=int)
        path[T-1] = np.argmax(log_delta[T-1, :])
        for t in range(T-2, -1, -1):
            path[t] = psi[t+1, path[t+1]]
        return path

    def fit(self, observations, positions, n_iter=10, tol=1e-4, verbose=True):

        observations = np.asarray(observations)
        positions = np.asarray(positions)
        T = len(observations)

        prev_loglik = -np.inf
        for iteration in range(n_iter):
            log_forward = self.forward_algorithm(observations, positions)
            log_backward = self.backward_algorithm(observations, positions)

            loglik = log_sum_exp(log_forward[-1, :])

            log_gamma = log_forward + log_backward
            max_log = np.max(log_gamma, axis=1, keepdims=True)
            gamma = np.exp(log_gamma - max_log)
            gamma /= gamma.sum(axis=1, keepdims=True)

            xi = np.zeros((T - 1, self.n_states, self.n_states))
            for t in range(T - 1):
                distance = positions[t+1] - positions[t]
                for i in range(self.n_states):
                    for j in range(self.n_states):
                        log_term = (
                            log_forward[t, i]
                            + self.transition_log_prob(i, j, distance)
                            + self.emission_log_prob(observations[t+1], j)
                            + log_backward[t+1, j]
                        )
                        xi[t, i, j] = np.exp(log_term - loglik)

                xi_sum = xi[t].sum()
                if xi_sum > 0:
                    xi[t] /= xi_sum

            for state_name, state_idx in zip(self.state_names, range(self.n_states)):
                weight = gamma[:, state_idx].sum()
                if weight == 0:
                    continue
                mean = np.sum(gamma[:, state_idx] * observations) / weight
                var = np.sum(gamma[:, state_idx] * (observations - mean) ** 2) / weight
                std = np.sqrt(max(var, 1e-6))
                self.emission_params[state_name]['mean'] = mean
                self.emission_params[state_name]['std'] = std

            stay_prob = np.sum(xi[:, 0, 0] + xi[:, 1, 1]) / np.sum(xi)
            avg_distance = np.mean(np.diff(positions))
            stay_prob = min(max(stay_prob, 0.51), 0.99) 
            new_lambda = -avg_distance / np.log(2 * stay_prob - 1)
            if np.isfinite(new_lambda) and new_lambda > 1:
                self.distance_scale = new_lambda
                self.transition_params['distance_scale'] = new_lambda

            if verbose:
                print(f"[EM] Iter {iteration+1}: loglik={loglik:.3f}, "
                      f"lambda={self.distance_scale:.1f}, "
                      f"neutral mean={self.emission_params['neutral']['mean']:.3f}, "
                      f"selection mean={self.emission_params['selection']['mean']:.3f}")

            if abs(loglik - prev_loglik) < tol:
                break
            prev_loglik = loglik

            print(f"[EM] Iter {iteration+1}: loglik={loglik:.3f}, "
                f"lambda={self.distance_scale:.1f}, "
                f"neutral mean={self.emission_params['neutral']['mean']:.3f}, "
                f"selection mean={self.emission_params['selection']['mean']:.3f}")

        return self


def log_sum_exp(log_values):
    max_log = np.max(log_values)
    return max_log + np.log(np.sum(np.exp(log_values - max_log)))
