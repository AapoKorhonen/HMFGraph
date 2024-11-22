#' Edge selection based on credible intervals
#'
#' @param posterior_samples Posterior samples from the HMFGraph_gibbs_sampler -function
#' @param CI The credible interval to be used for selection edges. On defaut the value is 0.9, or 90 % credible interval
#'
#' @return Returns an adjacency matrix, the posterior mean and median, and lower and upper credible interval points
#' @export
#'
#' @examples library(HMFGraph)
#' @examples n <- 200
#' @examples p <- 100
#' @examples set.seed(42)
#' @examples generated_data <- data_generator(n=n, p = p)
#' @examples results_HMFGraph_gibbs <- HMFGraph_gibbs_sampler(generated_data$data, alpha = p*5/(p*5+n), beta=0.9, iters = 5000, burn_in = 1000)
#' @examples results_gibbs_CI <- HMFGraph_gibbs_CI(results_HMFGraph_gibbs, CI = 0.8)
#' @examples library(qgraph)
#' @examples qgraph(results_gibbs_CI$adjacency_matrix)
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

