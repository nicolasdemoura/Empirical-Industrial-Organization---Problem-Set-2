###############################################################################
# Topic: Forward Simulation of Value Functions à la Hotz et al. (1994)
# Goal: To implement forward simulation method for computing ex-ante value functions
# Keywords: Empirical IO, Dynamic Games, Forward Simulation, Hotz-Miller
# Autor: Nícolas de Moura (<nicolasgoulartdemoura@gmail.com>) and Michel Wachsmann (<michel@wachsmann@gmail.com>)
# Date: 2025-12-08
###############################################################################

# Core function: Simulate one forward path from a starting state
# Returns: vector of discounted deterministic payoffs for each player
simulate_forward_path_hotz <- function(game, ccps, start_state, H) {
    discounted_payoffs <- rep(0, game$num_players)
    current_state <- start_state
    
    # Simulate forward for H periods
    for (t in 1:(H - 1)) {
        # Draw actions for each player according to CCPs
        actions <- numeric(game$num_players)
        for (p in 1:game$num_players) {
            actions[p] <- sample(game$actions, size = 1, prob = ccps[current_state, , p])
        }
        
        # Compute deterministic flow payoffs (no shocks)
        state_vector <- as.numeric(game$states[current_state, ])
        flow_payoffs <- game$payoff_function(actions, state_vector)
        
        # Add discounted flow payoffs
        discounted_payoffs <- discounted_payoffs + (game$beta^t) * flow_payoffs
        
        # Transition to next state (deterministic: s_{t+1} = a_t)
        action_profile_idx <- which(apply(game$action_profiles, 1, function(x) all(x == actions)))
        transition_row <- (current_state - 1) * game$num_action_profiles + action_profile_idx
        next_state_probs <- game$transition_matrix[transition_row, ]
        current_state <- which.max(next_state_probs)
    }
    
    return(discounted_payoffs)
}

# Compute forward-simulated continuation values for all states
# Returns: matrix of dimension num_states × num_players
compute_forward_continuation <- function(game, ccps, H, num_draws = 500) {
    V_forward <- matrix(0, nrow = game$num_states, ncol = game$num_players)
    
    # Use parallel processing if available
    if (requireNamespace("future.apply", quietly = TRUE)) {
        library(future.apply)
        if (inherits(future::plan(), "sequential")) {
            future::plan(future::multisession)
        }
        
        V_list <- future_lapply(1:game$num_states, function(s) {
            draws_payoffs <- matrix(0, nrow = num_draws, ncol = game$num_players)
            for (r in 1:num_draws) {
                draws_payoffs[r, ] <- simulate_forward_path_hotz(game, ccps, s, H)
            }
            return(colMeans(draws_payoffs))
        }, future.seed = TRUE)
        
        V_forward <- do.call(rbind, V_list)
    } else {
        # Sequential fallback
        for (s in 1:game$num_states) {
            draws_payoffs <- matrix(0, nrow = num_draws, ncol = game$num_players)
            for (r in 1:num_draws) {
                draws_payoffs[r, ] <- simulate_forward_path_hotz(game, ccps, s, H)
            }
            V_forward[s, ] <- colMeans(draws_payoffs)
        }
    }
    
    return(V_forward)
}
