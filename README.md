
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Hierarchical Matrix-F prior for Graphical Models, or HMFGraph

<!-- badges: start -->
<!-- badges: end -->

HMFGraph implements Bayesian Gaussian Graphical Model using hierarchical
matrix-F prior. A generalized expectation-maximization algorithm (GEM)
and a Gibbs sampler are included.

## Installation

You can install the package using the following commands:

``` r
# install.packages("devtools")
devtools::install_github("AapoKorhonen/HMFGraph")
```

## Example using FDR-control for edge selection

``` r
library(HMFGraph)

n <- 200
p <- 100

set.seed(42)
generated_data <- data_generator(n=n, p = p)

results_HMFGraph_GEM <- HMFGraph_GEM(generated_data$data, alpha = p*5/(p*5+n), beta=0.9)

permutations <- HMFGraph_GEM_permutations(generated_data$data, results_HMFGraph_GEM, number_of_permutations = 100, parallel = F)

results_FDR <- HMFGraph_GEM_FDR_control(results_HMFGraph_GEM, permutations, target_FDR = 0.2)

library(qgraph)

qgraph(results_FDR$adjacency_matrix)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

## Example using optimal CI for edge selection

``` r
library(HMFGraph)

n <- 200
p <- 100

set.seed(42)
generated_data <- data_generator(n=n, p = p)

results_HMFGraph_GEM <- HMFGraph_GEM(generated_data$data, alpha =  p*5/(p*5+n), beta=0.9)

permutations <- HMFGraph_GEM_permutations(generated_data$data, results_HMFGraph_GEM, number_of_permutations = 100, parallel = F)

results_optimal_CI <- HMFGraph_GEM_optimal_CI(results_HMFGraph_GEM, permutations, expected_connections = 100)

library(qgraph)

qgraph(results_optimal_CI$adjacency_matrix)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

## Example using CI for edge selection

``` r
library(HMFGraph)

n <- 200
p <- 100

set.seed(42)
generated_data <- data_generator(n=n, p = p)

results_HMFGraph_GEM <- HMFGraph_GEM(generated_data$data, alpha = p*5/(p*5+n), beta=0.9)

results_CI <- HMFGraph_GEM_CI(results_HMFGraph_GEM, CI = 0.8)

library(qgraph)

qgraph(results_CI$adjacency_matrix)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

## Example using Gibbs sampler with CI for edge selection

``` r
library(HMFGraph)

n <- 200
p <- 100

set.seed(42)
generated_data <- data_generator(n=n, p = p)

results_HMFGraph_gibbs <- HMFGraph_gibbs_sampler(generated_data$data, alpha = p*5/(p*5+n), beta=0.9, iters = 5000, burn_in = 1000)

results_gibbs_CI <- HMFGraph_gibbs_CI(results_HMFGraph_gibbs, CI = 0.8)

library(qgraph)

qgraph(results_gibbs_CI$adjacency_matrix)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />
