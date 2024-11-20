
#' Selects the edges for the graph based on an optimal credibility interval.
#'
#' @param HMFGraph_GEM_RESULTS  Results from the function HMFGraph_GEM
#' @param permutations   Results from the function HMFGraph_GEM_permutation
#' @param expected_connections  An expected number of connections in the real network. The default value is p, or the number of variables.
#'
#' @return  Returns the adjacency matrix, the map estimate, the variance matrix, lower and upper credibility interval point matrices.
#' @export
#'
#' @examples
HMFGraph_GEM_optimal_CI <- function(HMFGraph_GEM_RESULTS, permutations, expected_connections= NULL){
  
  #param MCC # If TRUE, then  the optimal credibility interval will be determined by an approximative MCC. Other wise an approximative F1 -value.
  
  if(is.null( expected_connections )){
    expected_connections <- HMFGraph_GEM_RESULTS$p
  }
  
  TP <- permutations$qp_connections - permutations$qp_connections_permutation 
  
  FP <- permutations$qp_connections_permutation
  
  FN <- expected_connections -TP
  
  
  F1 <- (2*TP)/(2*TP+FP+FN)
  
  N <- (HMFGraph_GEM_RESULTS$p^2/2 - HMFGraph_GEM_RESULTS$p)/2 - expected_connections
  
  TN <- N - permutations$qp_connections_permutation
  
  #MCC <- ( (TP*TN)-(FP*FN) )/ (sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)  )   )

  F1[is.nan(F1)] <- 0
  F1[is.na(F1)] <- 0
  
  # MCC[is.nan(MCC)] <- 0
  # MCC[is.na(MCC)] <- 0
  
  quantile_point <- permutations$quantile_points[F1==max(F1[!is.nan(F1)])][1]
  print(quantile_point)
  
  # if(MCC){
  #   quantile_point <- permutations$quantile_points[MCC==max(MCC[!is.nan(MCC)])][1]
  # }
  # 
  
  
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
