#' Edge selection based on credible intervals
#'
#' @param posterior_samples Posterior samples from the HMFGraph_gibbs_sampler -function
#' @param CI The credible interval to be used for selection edges. On defaut the value is 0.9, or 90 % credible interval
#'
#' @return Returns an adjacency matrix, the posterior mean and median, and lower and upper credible interval points
#' @export
#'
#' @examples 
#' 

HMFGraph_gibbs_CI <- function(posterior_samples, CI=0.9){
  
  a = (1-CI)/2
  
  gibbs_mean <- apply(posterior_samples$omega, c(1, 2),mean) 
  
  gibbs_median <- apply(posterior_samples$omega, c(1, 2),median) 
  
  
  lower_CI <- apply(posterior_samples$omega, c(1, 2),quantile,probs = a) 
  
  upper_CI <- apply(posterior_samples$omega, c(1, 2),quantile,probs = 1-a) 
  
  
  multiplication <- lower_CI*upper_CI
  
  adjacency_matrix <- lower_CI*upper_CI
  
  
  adjacency_matrix[multiplication <= 0] <- 0
  
  adjacency_matrix[multiplication > 0] <- 1
  
  
  return(list(adjacency_matrix = adjacency_matrix, median=  gibbs_median  , mean=  gibbs_mean, lower_CI = lower_CI, upper_CI = upper_CI) ) 
  
}

