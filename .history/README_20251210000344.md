# Empirical Industrial Organization - Problem Set 2

## Estimating Incomplete Information Dynamic Games: A Monte Carlo Comparison

**Authors:** Nícolas de Moura ([nicolasgoulartdemoura@gmail.com](mailto:nicolasgoulartdemoura@gmail.com)) and Michel Wachsmann ([michel@wachsmann@gmail.com](mailto:michel@wachsmann@gmail.com))

**Date:** December 2025

---

## Overview

This project implements and compares two structural estimation methods for dynamic games with incomplete information:

1. **Asymptotic Least Squares (ALS)** - Pesendorfer and Schmidt-Dengler (2008)
2. **Forward Simulation Method** - Hotz, Miller, Sanders, and Smith (1994)

The code conducts extensive Monte Carlo simulations to evaluate the finite-sample performance of these estimators across different sample sizes and simulation horizons.

## Project Structure

```
.
├── main.R                    # Main execution script with Monte Carlo simulations
├── model.R                   # Game setup, equilibrium solver, and core model functions
├── simulation.R              # Data generation from the dynamic game
├── estimation.R              # ALS estimator implementation
├── hotz_miller.R            # Forward simulation method (HMSS 1994)
├── performance.R            # Performance metrics and evaluation functions
├── markov_chain_viz.R       # Visualization utilities for Markov chains
└── output/                  # Simulation results and tables
    ├── *.rds                # Saved R data objects
    └── *.csv                # Summary tables
```
