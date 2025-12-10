###############################################################################
# Topic: Markov Chain Visualization for Dynamic Games
# Goal: Create LaTeX/TikZ visualization of equilibrium transition probabilities
# Keywords: Markov Chain, Equilibrium, Visualization
# Date: 2025-12-04
###############################################################################

# Function to generate LaTeX code for Markov chain diagram
generate_markov_chain_tex <- function(game, equilibrium, output_file = "output/markov_chain.tex") {
    
    # Compute transition probabilities induced by equilibrium
    transition_probs <- compute_equilibrium_transitions(game, equilibrium)
    
    # Compute ex-ante value functions at equilibrium
    beliefs_eq <- make_consistent_beliefs(game, equilibrium)
    V <- get_value_functions(game, beliefs_eq)
    
    # Start LaTeX document
    tex_code <- "\\documentclass[border=5pt]{standalone}\n"
    tex_code <- paste0(tex_code, "\\usepackage{tikz}\n")
    tex_code <- paste0(tex_code, "\\usetikzlibrary{positioning,arrows.meta,decorations.markings,bending}\n\n")
    tex_code <- paste0(tex_code, "\\begin{document}\n")
    tex_code <- paste0(tex_code, "\\begin{tikzpicture}[>=Stealth]\n\n")
    
    # Define node positions (4 states in square formation)
    tex_code <- paste0(tex_code, "  % Define states\n")
    tex_code <- paste0(tex_code, sprintf("  \\node[circle, draw, minimum size=1.5cm] (s00) at (0,0) {$(0,0)$};\n"))
    tex_code <- paste0(tex_code, sprintf("  \\node[circle, draw, minimum size=1.5cm] (s01) at (6,0) {$(0,1)$};\n"))
    tex_code <- paste0(tex_code, sprintf("  \\node[circle, draw, minimum size=1.5cm] (s10) at (0,-4) {$(1,0)$};\n"))
    tex_code <- paste0(tex_code, sprintf("  \\node[circle, draw, minimum size=1.5cm] (s11) at (6,-4) {$(1,1)$};\n\n"))
    
    # Add value functions as labels
    tex_code <- paste0(tex_code, "  % Value functions\n")
    state_names <- c("00", "01", "10", "11")
    positions <- list(c(-2, 0), c(8, 0), c(-2, -4), c(8, -4))
    
    for(s in 1:game$num_states) {
        v1 <- sprintf("%.2f", V[s, 1])
        v2 <- sprintf("%.2f", V[s, 2])
        pos <- positions[[s]]
        tex_code <- paste0(tex_code, sprintf("  \\node[align=left] at (%s,%s) {$V_1=%.2f$ \\\\ $V_2=%.2f$};\n", 
                                              pos[1], pos[2], V[s, 1], V[s, 2]))
    }
    tex_code <- paste0(tex_code, "\n")
    
    # Add transitions
    tex_code <- paste0(tex_code, "  % Transitions\n")
    
    # Map state indices to node names
    state_map <- c("s00", "s01", "s10", "s11")
    
    # Define loop positions (outward angles for each state)
    loop_angles <- list(
        c(45, 135),    # s00: top-right, top-left
        c(45, 135),    # s01: top-right, top-left  
        c(-135, -45),  # s10: bottom-left, bottom-right
        c(-135, -45)   # s11: bottom-left, bottom-right
    )
    
    # For each state, add arrows to other states (including self-loops)
    for(from_s in 1:game$num_states) {
        for(to_s in 1:game$num_states) {
            prob <- transition_probs[from_s, to_s]
            if(prob > 0.001) { # Only show transitions with probability > 0.1%
                from_node <- state_map[from_s]
                to_node <- state_map[to_s]
                
                # Scale line width based on probability (0.5pt to 3pt)
                line_width <- 0.5 + (prob * 2.5)
                
                if(from_s == to_s) {
                    # Self-loop (outward)
                    angles <- loop_angles[[from_s]]
                    tex_code <- paste0(tex_code, sprintf("  \\path[->,line width=%.2fpt] (%s) edge[out=%d,in=%d,looseness=8] node[auto,pos=0.5] {%.2f} (%s);\n", 
                                                        line_width, from_node, angles[1], angles[2], prob, from_node))
                } else {
                    # Regular transition
                    bend_dir <- ""
                    label_pos <- "above"
                    # Default pos is 0.5, but 0.33 for diagonal edges
                    pos_val <- 0.5
                    
                    # Horizontal edges (01 <-> 00 and 11 <-> 10)
                    if((from_s == 1 && to_s == 2) || (from_s == 2 && to_s == 1)) {
                        bend_dir <- if(from_s < to_s) "bend left=10" else "bend left=10"
                        label_pos <- if(from_s < to_s) "above" else "below"
                    } else if((from_s == 3 && to_s == 4) || (from_s == 4 && to_s == 3)) {
                        bend_dir <- if(from_s < to_s) "bend left=10" else "bend left=10"
                        label_pos <- if(from_s < to_s) "above" else "below"
                    }
                    # Vertical edges (00 <-> 10 and 01 <-> 11)
                    else if((from_s == 1 && to_s == 3) || (from_s == 3 && to_s == 1)) {
                        bend_dir <- if(from_s < to_s) "bend left=10" else "bend left=10"
                        label_pos <- if(from_s < to_s) "right" else "left"
                    } else if((from_s == 2 && to_s == 4) || (from_s == 4 && to_s == 2)) {
                        bend_dir <- if(from_s < to_s) "bend left=10" else "bend left=10"
                        label_pos <- if(from_s < to_s) "right" else "left"
                    }
                    # Diagonal edges: (0,0) <-> (1,1) and (0,1) <-> (1,0)
                    else {
                        bend_dir <- "bend left=5"
                        label_pos <- "above,sloped"
                        pos_val <- 0.33
                    }
                    
                    tex_code <- paste0(tex_code, sprintf("  \\draw[->,line width=%.2fpt] (%s) edge[%s] node[%s,fill=white,inner sep=2pt,pos=%.2f] {%.2f} (%s);\n", 
                                                        line_width, from_node, bend_dir, label_pos, pos_val, prob, to_node))
                }
            }
        }
    }
    
    tex_code <- paste0(tex_code, "\n\\end{tikzpicture}\n")
    tex_code <- paste0(tex_code, "\\end{document}\n")
    
    # Write to file
    writeLines(tex_code, output_file)
    cat(sprintf("Markov chain diagram saved to: %s\n", output_file))
    
    return(list(transition_probs = transition_probs, values = V))
}

# Compute transition probabilities induced by equilibrium
compute_equilibrium_transitions <- function(game, equilibrium) {
    trans_prob <- matrix(0, nrow = game$num_states, ncol = game$num_states)
    
    for(s in 1:game$num_states) {
        for(a in 1:game$num_action_profiles) {
            # Probability of action profile a in state s under equilibrium
            action_profile <- as.numeric(game$action_profiles[a, ])
            joint_prob <- equilibrium[s, action_profile[1] + 1, 1] * equilibrium[s, action_profile[2] + 1, 2]
            
            # Where does this action profile lead?
            transition_row <- (s - 1) * game$num_action_profiles + a
            for(s_next in 1:game$num_states) {
                trans_prob[s, s_next] <- trans_prob[s, s_next] + 
                    joint_prob * game$transition_matrix[transition_row, s_next]
            }
        }
    }
    
    rownames(trans_prob) <- paste0("S", sapply(1:game$num_states, function(x) paste(game$states[x, ], collapse = "")))
    colnames(trans_prob) <- paste0("S", sapply(1:game$num_states, function(x) paste(game$states[x, ], collapse = "")))
    
    return(trans_prob)
}

# Function to plot state sequence from time series data
# Creates a visualization showing the evolution of states over time
plot_state_sequence <- function(state_series, game = NULL, max_periods = 500, output_file = NULL) {
    library(ggplot2)
    
    # Limit to max_periods if series is longer
    if(length(state_series) > max_periods) {
        state_series <- state_series[1:max_periods]
        cat(sprintf("Note: Showing first %d periods of %d total\n", max_periods, length(state_series)))
    }
    
    # Create state labels
    if(!is.null(game)) {
        state_labels <- sapply(1:game$num_states, function(s) {
            paste0("(", paste(game$states[s, ], collapse = ","), ")")
        })
    } else {
        # Default for 4-state game: (0,0), (0,1), (1,0), (1,1)
        state_labels <- c("(0,0)", "(0,1)", "(1,0)", "(1,1)")
    }
    
    # Create data frame for plotting
    plot_data <- data.frame(
        Period = 1:length(state_series),
        State = factor(state_series, levels = 1:length(state_labels), labels = state_labels)
    )
    
    # Create plot
    p <- ggplot(plot_data, aes(x = Period, y = State)) +
        geom_line(aes(group = 1), color = "steelblue", linewidth = 0.8, alpha = 0.6) +
        geom_point(color = "steelblue", size = 2, alpha = 0.7) +
        theme_minimal(base_size = 12) +
        labs(
            x = "Period (t)",
            y = "State"
        ) +
        theme(
            panel.grid = element_blank(),
            axis.line.y = element_line(color = "black", linewidth = 0.5),
            axis.text.y = element_text(size = 11),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
        )
    
    # Save or display
    if(!is.null(output_file)) {
        ggsave(output_file, p, width = 10, height = 6, dpi = 300, bg = "white")
        cat(sprintf("State sequence plot saved to: %s\n", output_file))
    }
    
    # Also print state frequency table
    state_freq <- table(plot_data$State)
    
    return(p)
}



