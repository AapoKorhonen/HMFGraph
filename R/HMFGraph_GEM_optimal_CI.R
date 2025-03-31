
#' Selects the edges for the graph based on an optimal credible interval.
#'
#' @param HMFGraph_GEM_RESULTS  Results from the function HMFGraph_GEM
#' @param permutations   Results from the function HMFGraph_GEM_permutation
#' @param expected_connections  An expected number of connections in the real network. The default value is p, or the number of variables.
#' @param MCC If True then MCC will be used in the place of F1.
#'
#' @return  Returns the adjacency matrix, the map estimate, the variance matrix, lower and upper credible interval point matrices.
#' @export
#'
#' @examples library(HMFGraph)
#' @examples n <- 200
#' @examples p <- 100
#' @examples set.seed(42)
#' @examples generated_data <- data_generator(n=n, p = p)
#' @examples results_HMFGraph_GEM <- HMFGraph_GEM(generated_data$data, alpha = p * 5 / ( p * 5+n), beta=0.9)
#' @examples permutations <- HMFGraph_GEM_permutations(generated_data$data, results_HMFGraph_GEM, number_of_permutations = 100, parallel = F)
#' @examples results_optimal_CI <- HMFGraph_GEM_optimal_CI(results_HMFGraph_GEM, permutations, expected_connections = p)
#' @examples library(qgraph)
#' @examples qgraph(results_optimal_CI$adjacency_matrix)
HMFGraph_GEM_optimal_CI <- function(HMFGraph_GEM_RESULTS, permutations, expected_connections= NULL, MCC=F){
  
  #param MCC # If TRUE, then  the optimal credible interval will be determined by an approximative MCC. Other wise an approximative F1 -value.
  
  if(is.null( expected_connections )){
    expected_connections <- HMFGraph_GEM_RESULTS$p
  }
  
  TP <- permutations$qp_connections - permutations$qp_connections_permutation 
  
  FP <- permutations$qp_connections_permutation
  
  FN <- expected_connections - TP
  
  
  F1 <- (2*TP)/(2*TP+FP+FN)
  
  N <- (HMFGraph_GEM_RESULTS$p^2/2 - HMFGraph_GEM_RESULTS$p)/2 - expected_connections
  
  TN <- N - permutations$qp_connections_permutation

  F1[is.nan(F1)] <- 0
  F1[is.na(F1)] <- 0
  
  if(MCC){
    
    MCC <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    
    MCC[is.na(MCC)] <- 0
    F1 <- MCC
    
  }
  
  
  quantile_point <- permutations$quantile_points[F1==max(F1[!is.nan(F1)])][1]

  MAP_estimate <- HMFGraph_GEM_RESULTS$omega 
  
  variance_matrix <- HMFGraph_GEM_RESULTS$varmat 
  
  
  lower_CI <- MAP_estimate - quantile_point*sqrt(variance_matrix)
  
  upper_CI <- MAP_estimate + quantile_point*sqrt(variance_matrix)
  
  multiplication <- lower_CI*upper_CI
  
  adjacency_matrix <- lower_CI*upper_CI
  
  
  adjacency_matrix[multiplication <= 0] <- 0
  
  adjacency_matrix[multiplication > 0] <- 1
  
  diag(adjacency_matrix) <- 0
  
  adjacency_matrix[is.na(adjacency_matrix)] <- 0
  
  adjacency_matrix[is.nan(adjacency_matrix)] <- 0
  
  adjacency_matrix[is.null(adjacency_matrix)] <- 0
  
  return(list(adjacency_matrix = adjacency_matrix, MAP_estimate =  MAP_estimate,
              variance_matrix =  variance_matrix, lower_CI = lower_CI,
              upper_CI = upper_CI, F1=max(F1[!is.nan(F1)]), quantile_point =quantile_point ))  
}
