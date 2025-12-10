###############################################################################
# Topic: Estimating incomplete information dynamic games
# Goal: To generate functions to measure performance of estimation methods in dynamic games à la Pesendorfer and Schmidt-Dengler (2008)
# Keywords: Empirical IO, Dynamic Games, Incomplete Information, ALS, Structural Estimation 
# Autor: Nícolas de Moura (<nicolasgoulartdemoura@gmail.com>) and Michel Wachsmann (<michel@wachsmann@gmail.com>)
# Date: 2025-11-30
###############################################################################

# Function to measure performance of ALS estimation
# Mean, standard error, MSE
measure_als_performance <- function(true_params, estimated_params, plot = TRUE) {
    # Ensure estimated_params is a matrix (handle single vector case)
    if (is.vector(estimated_params)) {
        estimated_params <- matrix(estimated_params, nrow = 1)
    }
    
    # Calculate errors
    errors <- estimated_params - matrix(true_params, nrow = nrow(estimated_params), ncol = length(true_params), byrow = TRUE)
    
    # Parameter names based on number of parameters
    if (length(true_params) == 3) {
        param_names <- c("pi_0", "pi_1", "F")
    } else if (length(true_params) == 4) {
        param_names <- c("pi_0", "pi_1", "F", "beta")
    } else {
        param_names <- paste0("param_", 1:length(true_params))
    }
    
    # Calculate performance metrics
    performance <- data.frame(
        Parameter = param_names,
        True_Value = true_params,
        Mean_Estimate = colMeans(estimated_params),
        Std_Error = apply(estimated_params, 2, sd),
        MSE = colMeans(errors^2)
    )
    
    # Print performance metrics
    print(performance)
    
    # Plot estimates if required
    if (plot && nrow(estimated_params) > 1) {
        library(ggplot2)
        est_melted <- reshape2::melt(estimated_params)
        colnames(est_melted) <- c("Simulation", "Parameter", "Estimate")
        est_melted$Parameter <- param_names[est_melted$Parameter]
        
        ggplot(est_melted, aes(x = Parameter, y = Estimate)) +
            geom_boxplot() +
            geom_hline(data = performance, aes(yintercept = True_Value), color = "red", linetype = "dashed") +
            labs(title = "Parameter Estimates Distribution", y = "Estimates") +
            theme_minimal()
    }
    
    return(performance)
}