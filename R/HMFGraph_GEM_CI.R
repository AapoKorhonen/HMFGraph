#' Edge selection for GEM-MAP estimate based on a credible interval.
#'
#' @param HMFGraph_GEM_RESULTS  Results from the function HMFGraph_GEM
#' @param CI  Credible interval. The default value is 0.9, or 90% credible interval 
#'
#' @return  Returns the adjacency matrix, the map estimate, the variance matrix, lower and upper credible interval point matrices.
#' @export
#'
#' @examples
HMFGraph_GEM_CI <- function(HMFGraph_GEM_RESULTS, CI=0.9){
  
  a = (1-CI)/2
  
  MAP_estimate <- HMFGraph_GEM_RESULTS$omega 
  
  variance_matrix <- HMFGraph_GEM_RESULTS$varmat 
  
  # Quantile point z calculated using qnorm function. Each off diagonal element's 
  # posterior distributions is estimated to be normal distributions
  
  z <- qnorm(1-a)
  
  lower_CI <- MAP_estimate - z*sqrt(variance_matrix)
  
  upper_CI <- MAP_estimate + z*sqrt(variance_matrix)

  multiplication <- lower_CI*upper_CI
  
  adjacency_matrix <- lower_CI*upper_CI
  
  
  adjacency_matrix[multiplication <= 0] <- 0
  
  adjacency_matrix[multiplication > 0] <- 1
  
  diag(adjacency_matrix) <- 0
  
  return(list(adjacency_matrix = adjacency_matrix, MAP_estimate =  MAP_estimate,
              variance_matrix =  variance_matrix, lower_CI = lower_CI,
              upper_CI = upper_CI) ) 
  
}

