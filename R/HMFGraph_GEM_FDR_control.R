
#' Edge selection for GEM-Map estimate based on a target FDR.
#'
#' @param HMFGraph_GEM_RESULTS  Results from the function HMFGraph_GEM
#' @param permutation  Results from the function HMFGraph_GEM_permutation
#' @param target_FDR  The target FDR
#'
#' @return  Returns the adjacency matrix, the map estimate, the variance matrix, lower and upper credible interval point matrices.
#' @export
#'
#' @examples library(HMFGraph)
#' @examples 
#' @examples n <- 200
#' @examples p <- 100
#' @examples
#' @examples set.seed(42)
#' @examples generated_data <- data_generator(n=n, p = p)
#' @examples results_HMFGraph_GEM <- HMFGraph_GEM(generated_data$data, alpha = p * 5 / ( p * 5+n), beta=0.9)
#' @examples
#' @examples permutations <- HMFGraph_GEM_permutations(generated_data$data, results_HMFGraph_GEM, number_of_permutations = 100, parallel = F)
#' @examples
#' @examples results_FDR <- HMFGraph_GEM_FDR_control(results_HMFGraph_GEM, permutations, target_FDR = 0.2) library(qgraph)
#' @examples
#' @examples library(qgraph)
#' @examples qgraph(results_FDR$adjacency_matrix)

HMFGraph_GEM_FDR_control <- function(HMFGraph_GEM_RESULTS,permutation , target_FDR = 0.2 ){
  
  
  quantile_point <- permutation$quantile_points[permutation$qp_connections_permutation / permutation$qp_connections <= target_FDR][1]
  
  if(is.null(quantile_point)){
    quantile_point <- max(permutation$quantile_points)
  }
  
  print(quantile_point)
  
  MAP_estimate <- HMFGraph_GEM_RESULTS$omega 
  
  variance_matrix <- HMFGraph_GEM_RESULTS$varmat 
  
  
  lower_CI <- MAP_estimate - quantile_point*sqrt(variance_matrix)
  
  upper_CI <- MAP_estimate + quantile_point*sqrt(variance_matrix)
  
  multiplication <- lower_CI*upper_CI
  
  adjacency_matrix <- lower_CI*upper_CI
  
  
  adjacency_matrix[multiplication <= 0] <- 0
  
  adjacency_matrix[multiplication > 0] <- 1
  
  diag(adjacency_matrix) <- 0
  
  return(list(adjacency_matrix = adjacency_matrix, MAP_estimate =  MAP_estimate,
              variance_matrix =  variance_matrix, lower_CI = lower_CI,
              upper_CI = upper_CI) ) 
  
  
}
