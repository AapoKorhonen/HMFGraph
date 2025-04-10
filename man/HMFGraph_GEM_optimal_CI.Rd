% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HMFGraph_GEM_optimal_CI.R
\name{HMFGraph_GEM_optimal_CI}
\alias{HMFGraph_GEM_optimal_CI}
\title{Selects the edges for the graph based on an optimal credible interval.}
\usage{
HMFGraph_GEM_optimal_CI(
  HMFGraph_GEM_RESULTS,
  permutations,
  expected_connections = NULL
)
}
\arguments{
\item{HMFGraph_GEM_RESULTS}{Results from the function HMFGraph_GEM}

\item{permutations}{Results from the function HMFGraph_GEM_permutation}

\item{expected_connections}{An expected number of connections in the real network. The default value is p, or the number of variables.}
}
\value{
Returns the adjacency matrix, the map estimate, the variance matrix, lower and upper credible interval point matrices.
}
\description{
Selects the edges for the graph based on an optimal credible interval.
}
\examples{
library(HMFGraph)
n <- 200
p <- 100
set.seed(42)
generated_data <- data_generator(n=n, p = p)
results_HMFGraph_GEM <- HMFGraph_GEM(generated_data$data, alpha = p * 5 / ( p * 5+n), beta=0.9)
permutations <- HMFGraph_GEM_permutations(generated_data$data, results_HMFGraph_GEM, number_of_permutations = 100, parallel = F)
results_optimal_CI <- HMFGraph_GEM_optimal_CI(results_HMFGraph_GEM, permutations, expected_connections = p)
library(qgraph)
qgraph(results_optimal_CI$adjacency_matrix)
}
