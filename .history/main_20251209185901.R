###############################################################################
# Topic: Estimating incomplete information dynamic games
# Goal: Monte Carlo comparison of PSD (2008) ALS and HMSS (1994) estimators
# Keywords: Empirical IO, Dynamic Games, Incomplete Information, ALS, Structural Estimation 
# Autor: Nícolas de Moura (<nicolasgoulartdemoura@gmail.com>) and Michel Wachsmann (<michel@wachsmann@gmail.com>)
# Date: 2025-12-08
###############################################################################

###############################################################################
# Setup
###############################################################################

rm(list = ls())
suppressPackageStartupMessages({
    packages <- c("dplyr", "ggplot2", "data.table", "MASS", "nleqslv", "progress", "future.apply", "progressr", "reshape2")
    invisible(lapply(packages, function(pkg) {
        if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
            install.packages(pkg, dependencies = TRUE)
            library(pkg, character.only = TRUE)
        }
    }))
})

set.seed(20251208)

source("model.R")
source("simulation.R")
source("estimation.R")
source("performance.R")
source("hotz_miller.R")

###############################################################################
# Game Setup and Equilibrium
###############################################################################

pi_0 <- 1.2
pi_1 <- -1.2
F <- -0.2
beta <- 0.9
game <- get_game_parameters(pi_0 = pi_0, pi_1 = pi_1, F = F, beta = beta)

initial_beliefs <- array(runif(game$num_states * game$num_actions * game$num_players), 
                         dim = c(game$num_states, game$num_actions, game$num_players))
for (s in 1:game$num_states) {
    for (p in 1:game$num_players) {
        initial_beliefs[s, , p] <- initial_beliefs[s, , p] / sum(initial_beliefs[s, , p])
    }
}

equilibrium <- solve_conditional_action_probabilities(game, initial_beliefs)

###############################################################################
# Monte Carlo Simulation
###############################################################################

num_draws <- 1000
sample_sizes <- c(100, 1000, 10000, 100000)
initial_state <- 1
true_params_3 <- c(pi_0, pi_1, F)
true_params_4 <- c(pi_0, pi_1, F, beta)

for (T in sample_sizes) {
    cat(sprintf("\nT = %d\n\n", T))
    
    simulated_data <- simulate_game_data(game, draws = num_draws, num_periods = T, 
                                         initial_state = initial_state, equilibrium = equilibrium)
    
    ###########################################################################
    # PSD2008 Calibrated Beta (1000 draws, parallelized)
    ###########################################################################
    
    # cat("PSD2008 Calibrated Beta\n")
    
    # # Parallelize over draws
    # library(future.apply)
    # if (inherits(future::plan(), "sequential")) {
    #     future::plan(future::multisession)
    # }
    
    # pb <- progress_bar$new(format = "[:bar] :current/:total", total = num_draws, clear = FALSE, width = 60)
    
    # results_1_list <- future_lapply(1:num_draws, function(d) {
    #     single_draw_data <- list(
    #         draws = 1,
    #         states = matrix(simulated_data$states[d, ], nrow = 1),
    #         actions = array(simulated_data$actions[d, , ], dim = c(1, dim(simulated_data$actions)[2], dim(simulated_data$actions)[3]))
    #     )
        
    #     result <- estimate_parameters_als(single_draw_data, initial_params = c(0, 0, 0), 
    #                                       estimate_beta = FALSE, beta_fixed = beta)
    #     pb$tick()
    #     return(result)
    # }, future.seed = TRUE)
    
    # estimates_1 <- t(sapply(results_1_list, function(x) x$par))
    # convergence_1 <- sapply(results_1_list, function(x) x$convergence)
    # objectives_1 <- sapply(results_1_list, function(x) x$value)
    
    # perf_1 <- measure_als_performance(true_params_3, estimates_1, plot = FALSE)
    # cat("\n")
    # print(perf_1[, c("Parameter", "True_Value", "Mean_Estimate", "Std_Error", "MSE")])
    # cat("\n")
    
    # saveRDS(list(estimates = estimates_1, convergence = convergence_1, objectives = objectives_1, performance = perf_1),
    #         file = sprintf("output/PSD2008_Calibrated_Beta_T%d.rds", T))
    
    # ###########################################################################
    # # PSD2008 Estimated Beta (1000 draws, parallelized)
    # ###########################################################################
    
    # cat("PSD2008 Estimated Beta\n")
    
    # pb <- progress_bar$new(format = "[:bar] :current/:total", total = num_draws, clear = FALSE, width = 60)
    
    # results_2_list <- future_lapply(1:num_draws, function(d) {
    #     single_draw_data <- list(
    #         draws = 1,
    #         states = matrix(simulated_data$states[d, ], nrow = 1),
    #         actions = array(simulated_data$actions[d, , ], dim = c(1, dim(simulated_data$actions)[2], dim(simulated_data$actions)[3]))
    #     )
        
    #     result <- estimate_parameters_als(single_draw_data, initial_params = c(0, 0, 0, 0.5), 
    #                                       estimate_beta = TRUE, beta_fixed = NULL)
    #     pb$tick()
    #     return(result)
    # }, future.seed = TRUE)
    
    # estimates_2 <- t(sapply(results_2_list, function(x) x$par))
    # convergence_2 <- sapply(results_2_list, function(x) x$convergence)
    # objectives_2 <- sapply(results_2_list, function(x) x$value)
    
    # perf_2 <- measure_als_performance(true_params_4, estimates_2, plot = FALSE)
    # cat("\n")
    # print(perf_2[, c("Parameter", "True_Value", "Mean_Estimate", "Std_Error", "MSE")])
    # cat("\n")
    
    # saveRDS(list(estimates = estimates_2, convergence = convergence_2, objectives = objectives_2, performance = perf_2),
    #         file = sprintf("output/PSD2008_Estimated_Beta_T%d.rds", T))
    
    ###########################################################################
    # HMSS1994 - Multiple draws, multiple horizons (parallelize over draws)
    ###########################################################################
    
    horizons <- c(200, 100, 50)
    
    # Set up parallel backend for use within estimation functions
    library(future.apply)
    if (inherits(future::plan(), "sequential")) {
        future::plan(future::multisession)
    }
    
    for (H in horizons) {
        cat(sprintf("HMSS1994 (H = %d)\n", H))
        
        # Initialize storage
        estimates_HM <- matrix(NA, nrow = num_draws, ncol = 3)
        convergence_HM <- rep(NA, num_draws)
        objectives_HM <- rep(NA, num_draws)
        
        pb <- progress_bar$new(format = "[:bar] :current/:total", total = num_draws, clear = FALSE, width = 60)
        
        # Sequential loop over draws (parallelization happens inside estimation)
        for (d in 1:num_draws) {
            # Extract only this draw's data
            single_draw_data <- list(
                draws = 1,
                states = matrix(simulated_data$states[d, ], nrow = 1),
                actions = array(simulated_data$actions[d, , ], dim = c(1, dim(simulated_data$actions)[2], dim(simulated_data$actions)[3]))
            )
            
            # Estimate parameters (parallel computation happens inside)
            result <- estimate_parameters_hotz_miller(single_draw_data, 
                                                      initial_params = c(0, 0, 0),
                                                      game = game, 
                                                      beta_fixed = beta, 
                                                      H = H, 
                                                      num_sim_draws = 50)
            
            # Store results
            estimates_HM[d, ] <- result$par
            convergence_HM[d] <- result$convergence
            objectives_HM[d] <- result$value
            
            pb$tick()
        }
        
        perf_HM <- measure_als_performance(true_params_3, estimates_HM, plot = FALSE)
        cat("\n")
        print(perf_HM[, c("Parameter", "True_Value", "Mean_Estimate", "Std_Error", "MSE")])
        cat("\n")
        
        saveRDS(list(estimates = estimates_HM, convergence = convergence_HM, objectives = objectives_HM, performance = perf_HM),
                file = sprintf("output/HMSS1994_H%d_T%d.rds", H, T))
    }
}

cat("\n=== SIMULATION COMPLETE ===\n")

