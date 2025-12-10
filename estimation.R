###############################################################################
# Topic: Estimating incomplete information dynamic games
# Goal: To generate functions to estimate data from a dynamic game à la Pesendorfer and Schmidt-Dengler (2008)
# Keywords: Empirical IO, Dynamic Games, Incomplete Information, ALS, Structural Estimation 
# Autor: Nícolas de Moura (<nicolasgoulartdemoura@gmail.com>) and Michel Wachsmann (<michel@wachsmann@gmail.com>)
# Date: 2025-11-30
###############################################################################

# Function to get p_hat using frequentist approach (optimized)
get_p_hat <- function(simulated_data, game) {
    
    # Initialize p_hat array
    p_hat <- array(0, dim = c(simulated_data$draws, game$num_states, game$num_actions, game$num_players))
    
    # For each draw estimate the conditional action probabilities using frequency counts
    for (d in 1:simulated_data$draws) {
        states_d <- simulated_data$states[d, ]
        actions_d <- simulated_data$actions[d, , ]
        
        for (s in 1:game$num_states) {
            # Vectorized state matching - much faster than which()
            state_mask <- (states_d == s)
            total_count <- sum(state_mask)
            
            if (total_count > 0) {
                for (p in 1:game$num_players) {
                    # Vectorized action counting
                    actions_p <- actions_d[state_mask, p]
                    for (a in 1:game$num_actions) {
                        p_hat[d, s, a, p] <- sum(actions_p == a) / total_count
                    }
                }
            } else {
                # If state s was never visited, assign uniform probabilities
                p_hat[d, s, , ] <- 1 / game$num_actions
            }
        }
    }
    
    return(p_hat)
}

# Function to compute the ALS objective function
# This implements the Asymptotic Least Squares (ALS) estimator from
# Pesendorfer and Schmidt-Dengler (2008, Econometrica)
#
# ALS Objective: Q(θ) = Σ_s Σ_a Σ_i [BR_i(a|s; p̂, θ) - p̂_i(a|s)]² W(s,a)
#
# Key insight: Find parameters θ such that observed choice probabilities p̂
# satisfy the fixed-point condition: p̂ = BR(p̂; θ)
#
# This avoids solving for equilibrium at each θ, only computing best responses
als_objective <- function(params, simulated_data, weights = NULL, estimate_beta = FALSE, beta_fixed = 0.9) {
    # Extract structural parameters
    if(estimate_beta) {
        pi0 <- params[1]  # Payoff when opponent is inactive
        pi1 <- params[2]  # Payoff when opponent is active
        F <- params[3]    # Fixed cost of being active
        beta <- params[4] # Discount factor
    } else {
        pi0 <- params[1]  # Payoff when opponent is inactive
        pi1 <- params[2]  # Payoff when opponent is active
        F <- params[3]    # Fixed cost of being active
        beta <- beta_fixed # Use calibrated beta
    }
    
    # Construct game with candidate parameters
    game <- get_game_parameters(pi_0 = pi0, pi_1 = pi1, F = F, beta = beta)
    
    # Set weighting matrix (identity by default)
    if (is.null(weights)) {
        weights <- rep(1, game$num_states * game$num_actions * game$num_players)
    }
    
    # Step 1 (PSD Eq. 3): Estimate choice probabilities from data
    # p̂_i(a|s) = frequency of action a by player i in state s
    p_hat <- get_p_hat(simulated_data, game)
    p_hat_single <- p_hat[1, , , ]
    
    # Step 2: Form beliefs consistent with p̂
    # σ_{-i}(a_{-i}|s) = product of marginal probabilities
    beliefs <- make_consistent_beliefs(game, p_hat_single)
    
    # Step 3 (PSD Eq. 4-5): Compute best response probabilities
    # BR_i(a|s; p̂, θ) given value functions V_i(s; p̂, θ)
    # For logit shocks: BR_i(a|s) = Φ(V_i(s,a) - V_i(s,a'))
    best_response <- get_conditional_action_probabilities(game, beliefs)
    
    # Step 4 (PSD Eq. 6): Compute ALS objective
    # Q(θ) = ||BR(p̂; θ) - p̂||²_W
    als_value <- 0
    weight_idx <- 1
    for (s in 1:game$num_states) {
        for (p in 1:game$num_players) {
            for (a in 1:game$num_actions) {
                diff <- best_response[s, a, p] - p_hat[1, s, a, p]
                als_value <- als_value + diff^2 * weights[weight_idx]
                weight_idx <- weight_idx + 1
            }
        }
    }
    
    return(als_value)
}

# Function to estimate parameters using ALS (PSD 2008)
# Minimizes: θ̂ = argmin_θ Q(θ) = argmin_θ ||BR(p̂; θ) - p̂||²
estimate_parameters_als <- function(simulated_data, initial_params, weights = NULL, estimate_beta = FALSE, beta_fixed = 0.9) {
    result <- optim(par = initial_params, 
                    fn = als_objective, 
                    simulated_data = simulated_data, 
                    weights = weights,
                    estimate_beta = estimate_beta,
                    beta_fixed = beta_fixed,
                    method = "BFGS",
                    control = list(maxit = 500, trace = 0))
    
    return(list(par = result$par, value = result$value, convergence = result$convergence))
}

# Function to estimate parameters for each draw
estimate_parameters_by_draw <- function(simulated_data, initial_params, weights = NULL, estimate_beta = FALSE, beta_fixed = 0.9) {
    
    num_params <- length(initial_params)
    estimates <- matrix(NA, nrow = simulated_data$draws, ncol = num_params)
    convergence_codes <- rep(NA, simulated_data$draws)
    objective_values <- rep(NA, simulated_data$draws)
    
    for (d in 1:simulated_data$draws) {
        cat(sprintf("Draw %d/%d: ", d, simulated_data$draws))
        
        # Create data for this draw only
        single_draw_data <- list(
            draws = 1,
            states = matrix(simulated_data$states[d, ], nrow = 1),
            actions = array(simulated_data$actions[d, , ], dim = c(1, dim(simulated_data$actions)[2], dim(simulated_data$actions)[3]))
        )
        
        # Estimate using this draw
        result <- estimate_parameters_als(single_draw_data, initial_params = initial_params, weights = weights, 
                                         estimate_beta = estimate_beta, beta_fixed = beta_fixed)
        
        # Store results
        estimates[d, ] <- result$par
        convergence_codes[d] <- result$convergence
        objective_values[d] <- result$value
        
     }
    
    return(list(
        estimates = estimates,
        convergence = convergence_codes,
        objectives = objective_values
    ))
}

###############################################################################
# Hotz-Miller Two-Step Estimator (HMSS 1994) - Optimized Version
###############################################################################

# Pre-compute continuation value multipliers (done once per dataset)
# Returns: array[state, action, player, next_state] of discounted visit probabilities
compute_continuation_multipliers <- function(game, p_hat, H, num_sim_draws = 500) {
    cont_mult <- array(0, dim = c(game$num_states, game$num_actions, game$num_players, game$num_states))
    
    for (s in 1:game$num_states) {
        for (p in 1:game$num_players) {
            opponent_idx <- ifelse(p == 1, 2, 1)
            
            for (a_i in game$actions) {
                # Simulate forward H periods for this (state, action, player)
                state_visit_counts <- matrix(0, nrow = num_sim_draws, ncol = game$num_states)
                
                for (r in 1:num_sim_draws) {
                    # Start from state s, take action a_i
                    # First, determine next state from this initial action
                    
                    # Sample opponent's action in initial period
                    a_opp_init <- sample(game$actions, size = 1, prob = p_hat[s, , opponent_idx])
                    
                    # Construct action profile for initial transition
                    actions_init <- numeric(game$num_players)
                    actions_init[p] <- a_i
                    actions_init[opponent_idx] <- a_opp_init
                    
                    # Transition to first next state
                    action_profile_idx_init <- which(apply(game$action_profiles, 1, 
                                                     function(x) all(x == actions_init)))
                    transition_row_init <- (s - 1) * game$num_action_profiles + action_profile_idx_init
                    next_state_probs_init <- game$transition_matrix[transition_row_init, ]
                    current_state <- which.max(next_state_probs_init)
                    
                    # Now simulate forward H-1 periods from this next state
                    # Using CCPs for BOTH players (no longer conditioning on a_i)
                    for (t in 1:(H - 1)) {
                        # Record discounted visit to current state
                        # Discount is β^t because this is t periods ahead
                        state_visit_counts[r, current_state] <- state_visit_counts[r, current_state] + game$beta^t
                        
                        # Sample both players' actions according to CCPs
                        a_p_next <- sample(game$actions, size = 1, prob = p_hat[current_state, , p])
                        a_opp_next <- sample(game$actions, size = 1, prob = p_hat[current_state, , opponent_idx])
                        
                        # Construct action profile
                        actions_next <- numeric(game$num_players)
                        actions_next[p] <- a_p_next
                        actions_next[opponent_idx] <- a_opp_next
                        
                        # Transition to next state
                        action_profile_idx <- which(apply(game$action_profiles, 1, 
                                                         function(x) all(x == actions_next)))
                        transition_row <- (current_state - 1) * game$num_action_profiles + action_profile_idx
                        next_state_probs <- game$transition_matrix[transition_row, ]
                        current_state <- which.max(next_state_probs)
                    }
                }
                
                # Average across simulation draws
                cont_mult[s, a_i + 1, p, ] <- colMeans(state_visit_counts)
            }
        }
    }
    
    return(cont_mult)
}

# Hotz-Miller objective (optimized - uses pre-computed multipliers)
hotz_miller_objective_fast <- function(params, p_hat, cont_mult, game, beta_fixed = 0.9) {
    pi0 <- params[1]
    pi1 <- params[2]
    F <- params[3]
    W <- 0.1  # Scrap value (fixed)
    
    choice_values <- array(0, dim = c(game$num_states, game$num_actions, game$num_players))
    
    for (s in 1:game$num_states) {
        state_vector <- as.numeric(game$states[s, ])
        
        for (p in 1:game$num_players) {
            opponent_idx <- ifelse(p == 1, 2, 1)
            
            for (a_i in game$actions) {
                # Immediate payoff E_{a_{-i}}[u(s, a_i, a_{-i})]
                immediate_payoff <- 0
                
                for (a_opp in game$actions) {
                    actions <- numeric(game$num_players)
                    actions[p] <- a_i
                    actions[opponent_idx] <- a_opp
                    
                    # Flow payoff with current parameters
                    # Note: The payoff function is:
                    # u_p = a_p * (pi0*(1-a_{-p}) + pi1*a_{-p}) + a_p*(1-s_p)*F + (1-a_p)*s_p*W
                    u <- a_i * (pi0 * (1 - a_opp) + pi1 * a_opp) + 
                         a_i * (1 - state_vector[p]) * F + 
                         (1 - a_i) * state_vector[p] * W
                    
                    opp_prob <- p_hat[s, a_opp + 1, opponent_idx]
                    immediate_payoff <- immediate_payoff + opp_prob * u
                }
                
                # Continuation value (using pre-computed multipliers)
                continuation_value <- 0
                for (s_next in 1:game$num_states) {
                    state_next_vector <- as.numeric(game$states[s_next, ])
                    
                    # Expected flow payoff in next state
                    expected_next_payoff <- 0
                    for (a_next_i in game$actions) {
                        for (a_next_opp in game$actions) {
                            actions_next <- numeric(game$num_players)
                            actions_next[p] <- a_next_i
                            actions_next[opponent_idx] <- a_next_opp
                            
                            # Flow payoff in next state
                            u_next <- a_next_i * (pi0 * (1 - a_next_opp) + pi1 * a_next_opp) + 
                                     a_next_i * (1 - state_next_vector[p]) * F + 
                                     (1 - a_next_i) * state_next_vector[p] * W
                            
                            prob_next <- p_hat[s_next, a_next_i + 1, p] * p_hat[s_next, a_next_opp + 1, opponent_idx]
                            expected_next_payoff <- expected_next_payoff + prob_next * u_next
                        }
                    }
                    
                    continuation_value <- continuation_value + cont_mult[s, a_i + 1, p, s_next] * expected_next_payoff
                }
                
                choice_values[s, a_i + 1, p] <- immediate_payoff + continuation_value
            }
        }
    }
    
    # Compute model-implied CCPs
    implied_ccps <- array(0, dim = c(game$num_states, game$num_actions, game$num_players))
    
    for (s in 1:game$num_states) {
        for (p in 1:game$num_players) {
            v_diff <- choice_values[s, 2, p] - choice_values[s, 1, p]
            implied_ccps[s, 2, p] <- pnorm(v_diff)
            implied_ccps[s, 1, p] <- 1 - implied_ccps[s, 2, p]
        }
    }
    
    # Compute objective
    objective_value <- 0
    for (s in 1:game$num_states) {
        for (p in 1:game$num_players) {
            for (a in 1:game$num_actions) {
                diff <- implied_ccps[s, a, p] - p_hat[s, a, p]
                objective_value <- objective_value + diff^2
            }
        }
    }
    
    return(objective_value)
}

# Estimate parameters using optimized Hotz-Miller
estimate_parameters_hotz_miller <- function(simulated_data, initial_params, game, 
                                           beta_fixed = 0.9, H = 100, num_sim_draws = 500) {
    # Step 1: Estimate CCPs (done once)
    p_hat <- get_p_hat(simulated_data, game)
    p_hat_single <- p_hat[1, , , ]
    
    # Laplace smoothing
    for (s in 1:game$num_states) {
        for (p in 1:game$num_players) {
            n_total <- sum(simulated_data$states[1, ] == s)
            if (n_total > 0) {
                for (a in 1:game$num_actions) {
                    p_hat_single[s, a, p] <- (p_hat_single[s, a, p] * n_total + 1) / (n_total + 2)
                }
                p_hat_single[s, , p] <- p_hat_single[s, , p] / sum(p_hat_single[s, , p])
            }
        }
    }
    
    # Step 2: Pre-compute continuation multipliers (done once)
    cont_mult <- compute_continuation_multipliers(game, p_hat_single, H, num_sim_draws)
    
    # Step 3: Optimize (fast - no simulation in objective)
    result <- optim(par = initial_params,
                    fn = hotz_miller_objective_fast,
                    p_hat = p_hat_single,
                    cont_mult = cont_mult,
                    game = game,
                    beta_fixed = beta_fixed,
                    method = "BFGS",
                    control = list(maxit = 500, trace = 0))
    
    return(list(par = result$par, value = result$value, convergence = result$convergence))
}


