###############################################################################
# Topic: Estimating incomplete information dynamic games
# Goal: To generate functions to solve the model 
# Keywords: Empirical IO, Dynamic Games, Incomplete Information, ALS, Structural Estimation
# Autor: Nícolas de Moura (<nicolasgoulartdemoura@gmail.com>) and Michel Wachsmann (<michel@wachsmann@gmail.com>)
# Date: 2025-11-27
###############################################################################

# Function to set up the game environment
get_game_parameters <- function(pi_0 = 1.2, pi_1 = -1.2, F = 0.2, beta = 0.9) {

    # Set the players 
    num_players <- 2
    players <- 1:num_players

    # Set the actions 
    num_actions <- 2
    actions <- 0:(num_actions - 1)

    # Set the action profiles
    num_action_profiles <- num_actions ^ num_players
    action_profiles <- expand.grid(rep(list(0:(num_actions - 1)), num_players))

    # State profiles
    # Set the states
    num_states <- num_actions ^ num_players
    states <- expand.grid(rep(list(0:(num_actions - 1)), num_players))

    # Set discount factor
    beta <- beta

    # Set parameters for payoffs and transitions
    pi_0 <- pi_0                # Payoff if state is (.,0)
    pi_1 <- pi_1                # Payoff if state is (.,1)
    F <- F                      # Fixed cost
    W <- 0.1                    # Scrap value

    # Payoff function
    payoff_function <- function(actions, state) {
        payoff <- matrix(0, nrow = num_players, ncol = 1)
        payoff[1] <- actions[1] * (pi_0 * (1 - actions[2]) + pi_1 * actions[2]) + actions[1] * (1 - state[1]) * F + (1 - actions[1]) * state[1] * W
        payoff[2] <- actions[2] * (pi_0 * (1 - actions[1]) + pi_1 * actions[1]) + actions[2] * (1 - state[2]) * F + (1 - actions[2]) * state[2] * W
        return(payoff)
    }

    # Transition matrix such that s_{t+1} = a_t with probability 1, independent of current state
    transition_matrix <- matrix(0, nrow = num_action_profiles * num_states, ncol = num_states)
    for(s in 1:num_states) {
        for(a in 1:num_action_profiles) {
            next_state <- as.numeric(action_profiles[a, ])
            next_state_index <- which(apply(states, 1, function(x) all(x == next_state)))
            transition_matrix[(s - 1) * num_action_profiles + a, next_state_index] <- 1
        }
    }
    colnames(transition_matrix) <- paste0("S", sapply(1:num_states, function(x) paste(states[x, ], collapse = "")))
    rownames(transition_matrix) <- paste0("S", sapply(rep(1:num_states, each = num_action_profiles) - 1, function(x) paste(states[x + 1, ], collapse = "")), "_A", sapply(rep(1:num_action_profiles, times = num_states) - 1, function(x) paste(action_profiles[x + 1, ], collapse = "")))

    # Placeholder for game setup
    game <- list(
        num_players = num_players,
        players = players,
        num_actions = num_actions, 
        actions = actions,
        num_action_profiles = num_action_profiles,
        action_profiles = action_profiles,
        num_states = num_states,
        states = states,
        beta = beta,
        payoff_function = payoff_function,
        transition_matrix = transition_matrix        
    )
    return(game)
}


# Set up expected flow payoffs matrix
# Pi_i(s, a) = flow payoff of player i when action profile a is played in state s
get_flow_payoffs <- function(game) {
    Pi <- matrix(0, nrow = game$num_action_profiles * game$num_states, ncol = game$num_players)

    # Calculate flow payoffs
    for(s in 1:game$num_states) {
        for(a in 1:game$num_action_profiles) {
            actions <- as.numeric(game$action_profiles[a, ])
            state <- as.numeric(game$states[s, ])
            payoff <- game$payoff_function(actions, state)
            Pi[(s - 1) * game$num_action_profiles + a, ] <- payoff
        }
    }
    colnames(Pi) <- game$players
    rownames(Pi) <- paste0("S", sapply(rep(1:game$num_states, each = game$num_action_profiles) - 1, function(x) paste(game$states[x + 1, ], collapse = "")), "_A", sapply(rep(1:game$num_action_profiles, times = game$num_states) - 1, function(x) paste(game$action_profiles[x + 1, ], collapse = "")))
    return(Pi)
}

# Set up expected expected shocks matrix - for logit shocks this is Euler's constant
get_expected_shocks <- function(game, beliefs) {
    D <- matrix(0, nrow = game$num_states, ncol = game$num_players)
    for(p in 1:game$num_players) {
        for(s in 1:game$num_states) {
            marginal_belief <- rep(0, game$num_actions)
            for(a in 1:game$num_action_profiles) {
                action_profile <- as.numeric(game$action_profiles[a, ])
                marginal_belief[action_profile[p] + 1] <- marginal_belief[action_profile[p] + 1] + beliefs[s, a, p]
            }
            # For logit shocks: E[eps_i | action i] = -log(P(action i)) - Euler constant
            # But in the normalized logit model, expected shocks integrate to 0
            D[s, p] <- 0
        }
    }
    return(D)
}

# Function to reshape a belief to a num_states x (num_action_profiles x num_states) matrix with off-diagonal blocks of zeros for a player p
reshape_belief_for_player <- function(belief, game, p) {
    reshaped_belief <- matrix(0, nrow = game$num_states, ncol = game$num_action_profiles * game$num_states)
    for(s in 1:game$num_states) {
        for(a in 1:game$num_action_profiles) {
            reshaped_belief[s, (s - 1) * game$num_action_profiles + a] <- belief[s, a, p]
        }
    }
    return(reshaped_belief)
}

# Function to solve for the value functions given parameters and beliefs
# Implements PSD (2008) Equation 4:
# V_i(s; σ, θ) = E[u_i(s,a; θ) + ε_i(a) + β Σ_{s'} V_i(s'; σ, θ) Pr(s'|s,a,σ) | σ, ε]
# Solving: V_i = (I - β Γ(σ))^{-1} [Π_i(σ) + D_i(σ)]
# where Γ(σ) is transition matrix, Π_i(σ) is expected flow payoff, D_i(σ) is expected shock
get_value_functions <- function(game, beliefs) {
    Pi <- get_flow_payoffs(game)
    D <- get_expected_shocks(game, beliefs)
    V <- matrix(0, nrow = game$num_states, ncol = game$num_players)

    for(p in 1:game$num_players) {
        sigma <- reshape_belief_for_player(beliefs, game, p)
        mat_to_invert <- diag(game$num_states) - game$beta * (sigma %*% game$transition_matrix)
        first_term <- solve(mat_to_invert)
        second_term <- (sigma %*% Pi[, p]) + D[, p]
        V[, p] <- first_term %*% second_term
    }
    colnames(V) <- game$players
    rownames(V) <- paste0("S", sapply(1:game$num_states, function(x) paste(game$states[x, ], collapse = "")))
    return(V)
}

# Compute conditional choice probabilities (best response mapping)
# Implements PSD (2008) Equation 5:
# For normal shocks: σ_i(a|s; θ) = Φ(V_i(s,a; θ) - V_i(s,a'; θ))
# This is the best response probability given value functions
get_conditional_action_probabilities <- function(game, beliefs) {
    V <- get_value_functions(game, beliefs)

    P <- array(0, dim = c(game$num_states, game$num_actions, game$num_players))
    for(p in 1:game$num_players) {
        for(s in 1:game$num_states) {
            u <- expected_value_NoS(game, beliefs, V, p)
            
            # Compute probabilities using normal distribution as P(a_i = 0 | s) = Phi(u(s,a_i=0) - u(s,a_i=1))
            for(a_i in 1:game$num_actions) {
                if(a_i == 1) {
                    P[s, a_i, p] <- pnorm(u[s, 1] - u[s, 2])
                } else {
                    P[s, a_i, p] <- 1 - pnorm(u[s, 1] - u[s, 2])
                }
            }
        }
    }
    
    return(P)
}

expected_value_NoS <- function(game, beliefs, V, p) {
    # Compute expected value function for player p choosing action a_i
    # given beliefs about opponents' strategies (net of own shocks)
    # V(s, a_i) = E_{a_{-i}}[u_i(s, a_i, a_{-i}) + β E[V(s')|s, a_i, a_{-i}]]
    u <- matrix(0, nrow = game$num_states, ncol = game$num_actions)
    for(s in 1:game$num_states) {
        for(a_i in 1:game$num_actions) {
            expected_value <- 0
            # Sum over all action profiles where player p takes action a_i
            for(a in 1:game$num_action_profiles) {
                action_profile <- as.numeric(game$action_profiles[a, ])
                # Check if this action profile has player p taking action a_i
                if(action_profile[p] == game$actions[a_i]) {
                    # Flow payoff for this action profile
                    flow_payoff <- game$payoff_function(action_profile, as.numeric(game$states[s, ]))[p]
                    
                    # Expected continuation value
                    expected_next_value <- 0
                    for(s_next in 1:game$num_states) {
                        transition_prob <- game$transition_matrix[(s - 1) * game$num_action_profiles + a, s_next]
                        expected_next_value <- expected_next_value + transition_prob * V[s_next, p]
                    }
                    
                    # Weight by probability of opponents playing their part of this action profile
                    # This is beliefs[s, a, p] which represents the joint probability from player p's perspective
                    prob_this_profile <- beliefs[s, a, p]
                    expected_value <- expected_value + prob_this_profile * (flow_payoff + game$beta * expected_next_value)
                }
            }
            
            # Normalize by the marginal probability of player p taking action a_i
            # to get the expected value conditional on taking action a_i
            marginal_prob_ai <- 0
            for(a in 1:game$num_action_profiles) {
                action_profile <- as.numeric(game$action_profiles[a, ])
                if(action_profile[p] == game$actions[a_i]) {
                    marginal_prob_ai <- marginal_prob_ai + beliefs[s, a, p]
                }
            }
            
            if(marginal_prob_ai > 1e-10) {
                u[s, a_i] <- expected_value / marginal_prob_ai
            } else {
                u[s, a_i] <- 0
            }
        }
    }
    colnames(u) <- game$actions
    rownames(u) <- paste0("S", sapply(1:game$num_states, function(x) paste(game$states[x, ], collapse = "")))
    return(u)
}

# Make beliefs consistent with own action probabilities
# own_beliefs is a num_states x num_actions x num_players matrix
make_consistent_beliefs <- function(game, own_beliefs) {
    beliefs <- array(0, dim = c(game$num_states, game$num_action_profiles, game$num_players))
    for(p in 1:game$num_players) {
        for(s in 1:game$num_states) {
            for(a in 1:game$num_action_profiles) {
                action_profile <- as.numeric(game$action_profiles[a, ])
                prob <- 1
                for(p_other in 1:game$num_players) {
                    prob <- prob * own_beliefs[s, action_profile[p_other] + 1, p_other]
                }
                beliefs[s, a, p] <- prob
            }
        }
    }
    return(beliefs)
}

# Fixed point nleqslv to find conditional action probabilities
# Logistic transforms
logit    <- function(p) {
    # Bound p away from 0 and 1 to avoid infinite values
    p <- pmax(pmin(p, 1 - 1e-10), 1e-10)
    log(p/(1 - p))
}
invlogit <- function(t) 1 / (1 + exp(-t))

# Fixed point solver for conditional action probabilities
solve_conditional_action_probabilities <- function(game, initial_own_beliefs) {

    # Flatten initial beliefs and convert to theta-space
    p0 <- as.vector(initial_own_beliefs)
    t0 <- logit(p0)   # unconstrained parameter

    # Objective defined in theta-space
    objective_function_theta <- function(t_vec) {
        p_vec <- invlogit(t_vec)  # convert back into probabilities ∈ (0,1)

        # reshape into structure: states × actions × players
        own_beliefs <- array(p_vec,
                             dim = c(game$num_states, game$num_actions, game$num_players))

        beliefs <- make_consistent_beliefs(game, own_beliefs)
        P <- get_conditional_action_probabilities(game, beliefs)

        # fixed-point equation: P - own_beliefs
        residual <- as.vector(P - own_beliefs)

        return(residual)  # still on probability scale
    }

    # Run nleqslv in theta-space
    result <- nleqslv(
        t0,
        objective_function_theta,
        method = "Broyden",
        control = list(ftol = 1e-18, xtol = 1e-18, maxit = 500),
        jacobian = FALSE,
        global = "dbldog"
    )

    # Convert solution back to probability space
    p_star <- invlogit(result$x)

    # reshape to array
    consistent_own_beliefs <- array(
        p_star,
        dim = c(game$num_states, game$num_actions, game$num_players)
    )

    return(consistent_own_beliefs)
}
