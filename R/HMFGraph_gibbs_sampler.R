
#' Gibbs sampler for Hierarchical matrix F using Rcpp
#'
#' @param data The data matrix. Variables on columns and samples on rows. n x p
#' @param p The number of variables. The value will be derived from the data matrix dimensions if not specified.
#' @param n The number of samples. The value will be derived from the data matrix dimensions if not specified.
#' @param alpha The value for hyperparameter alpha
#' @param beta The value for hyperparameter beta
#' @param iters The number of iterarions
#' @param burn_in Burn in length for Gibbs sampling procedure. 
#' @param epsilon1 A shape parameter for Inverse_gamma prior. If a flat prior on logarithmic scale is used, this values should be set to 0.
#' @param epsilon2 A scale parameter for Inverse_gamma prior. If a flat prior on logarithmic scale is used, this values should be set to 0.
#' @param B The initial value for target matrix B. If fixed_B is TRUE, then this value will not change for the whole sampling process.
#' @param fixed_B Boolean value. If TRUE, then no prior for target matrix B is used.
#' @param print_t If FALSE, no progress bar will be printed.
#' 
#' @return Posterior samples of Omega. Also returns nu, delta, alpha and beta values, and the last samples of B and Phi matrices.
#' @export
#'
#' @examples library(HMFGraph)
#' @examples n <- 200
#' @examples p <- 100
#' @examples set.seed(42)
#' @examples generated_data <- data_generator(n=n, p = p)
#' @examples results_HMFGraph_gibbs <- HMFGraph_gibbs_sampler(generated_data$data, alpha = p*5/(p*5+n), beta=0.9, iters = 5000, burn_in = 1000)
HMFGraph_gibbs_sampler <- function(data, p = 0,  n = 0, 
                              alpha = 0.9,  beta = 0.9, iters = 5000, burn_in = 1000, epsilon1 = 0.001, epsilon2 = 0.001, B = NULL,fixed_B = F, print_t=T,
                              a_lim = 0.0001,  b_lim =10^5, omega0=NULL){
  
  
  # If n is not specified, the number of samples are derived from the dimension of the data matrix
  if(n == 0){
    n = dim(data)[1]
  }
  
  if(p == 0){
    p = dim(data)[2]
  }
  
  if(is.null(B)){
    B = diag(p)
  }
  
  if(is.null(omega0)){
    omega0 = diag(p)
  }
  
  nu <- alpha_to_nu(p = p, n = n , alpha = alpha)
  delta <- beta_to_delta(p = p,n = n, nu = nu , beta = beta)
  
  
  posterior_samples <-  gibbs_algorithm_cpp(iters+burn_in, cov(data),B= B, p =p, n = n, 
                         delta = delta,  nu = nu, epsilon1 = epsilon1 , epsilon2 = epsilon2, fixed_B, print_t=print_t,
                         a_lim, b_lim, Omega0=omega0)
  
  # Removing the burn in
  posterior_samples$omega <- posterior_samples$omega[,, (burn_in+1):(iters+burn_in) ]  
  
  posterior_samples <- c(posterior_samples, beta = beta)
  
  posterior_samples <- c(posterior_samples, alpha = alpha)
  
  return(posterior_samples)
}
