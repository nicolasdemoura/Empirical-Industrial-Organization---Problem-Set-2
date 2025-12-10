###############################################################################
# Topic: Estimating incomplete information dynamic games
# Goal: To generate functions to simulate data from a dynamic game à la Pesendorfer and Schmidt-Dengler (2008)
# Keywords: Empirical IO, Dynamic Games, Incomplete Information, ALS, Structural Estimation 
# Autor: Nícolas de Moura (<nicolasgoulartdemoura@gmail.com>) and Michel Wachsmann (<michel@wachsmann@gmail.com>)
# Date: 2025-11-29
###############################################################################

# Function to simulate data from the dynamic game
simulate_game_data <- function(game, draws, num_periods, initial_state, equilibrium) {
    num_players <- game$num_players
    num_states <- game$num_states
    num_actions <- game$num_actions
    
    # Initialize storage for states and actions
    states <- matrix(0, nrow = draws, ncol = num_periods)
    actions <- array(0, dim = c(draws, num_periods, num_players))
    
    # Set initial states
    states[, 1] <- initial_state
    
    for (d in 1:draws) {
        for (t in 1:(num_periods - 1)) {
            current_state <- states[d, t]
            
            # Sample actions based on equilibrium probabilities
            for (p in 1:num_players) {
                action_probs <- equilibrium[current_state, , p]
                actions[d, t, p] <- sample(1:num_actions, size = 1, prob = action_probs)
            }
            
            # Determine next state based on transition matrix
            action_profile <- actions[d, t, ]
            # Find the action profile index
            action_profile_idx <- which(apply(game$action_profiles, 1, function(x) all(x == (action_profile - 1))))
            # Get transition probabilities for current state and action profile
            transition_row_idx <- (current_state - 1) * game$num_action_profiles + action_profile_idx
            next_state_probs <- game$transition_matrix[transition_row_idx, ]
            states[d, t + 1] <- sample(1:num_states, size = 1, prob = next_state_probs)
        }
        
        # Sample actions for the last period
        current_state <- states[d, num_periods]
        for (p in 1:num_players) {
            action_probs <- equilibrium[current_state, , p]
            actions[d, num_periods, p] <- sample(1:num_actions, size = 1, prob = action_probs)
        }
    }
    
    return(list(draws = draws, states = states, actions = actions))
}

