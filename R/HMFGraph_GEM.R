#' GEM-algorithm for solving MAP-estimate for Bayesian GGM using hierarchical matrix-f prior
#'
#' @param data The data matrix. Variables on columns and samples on rows. n x p
#' @param p The number of variables. The value will be derived from the data matrix dimensions if not specified.
#' @param n The number of samples. The value will be derived from the data matrix dimensions if not specified.
#' @param alpha The value for hyperparameter alpha
#' @param beta The value for hyperparameter beta
#' @param max_iters The number of maximum iterarions
#' @param stop_criterion Stopping criterion for GEM algorithm
#' @param epsilon1 A shape parameter for Inverse_gamma prior. If a flat prior on logarithmic scale is used, this values should be set to 0.
#' @param epsilon2 A scale parameter for Inverse_gamma prior. If a flat prior on logarithmic scale is used, this values should be set to 0.
#' @param B The initial value for target matrix B. If fixed_B is TRUE, then this value will not change for the whole sampling process.
#' @param fixed_B Boolean value. If TRUE, then no prior for target matrix B is used.
#' @param inter Interval of iterations when info about algorithm is printed out.
#' @param print_t If FALSE, then no text is printed during the function running process.
#' 
#' @return Returns the map estimate for precision matrix. Also returns all values that are used for the algorithm.
#' @export
#'
#' @examples library(HMFGraph)
#' @examples n <- 200
#' @examples p <- 100
#' @examples set.seed(42)
#' @examples generated_data <- data_generator(n=n, p = p)
HMFGraph_GEM <- function(data, p = 0,  n = 0, 
                    alpha = -1,  beta = 0.9, max_iters = 5000,
                    stop_criterion = 10^(-6)  ,epsilon1 = 0, epsilon2 = 0,
                    B = NULL,fixed_B = F, inter=500,  print_t=T){
  
  


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
  
  if(alpha == -1){
    
    fn <- function(x){
      tulos <- logmarginal(x,n,p,eigen(  diag(p), only.values = T)$values, eigen(  solve(diag(p)) %*% (n*cov(scale(data))), only.values = T)$values)
      return(Re(tulos))
    }
    
    nu <- optimize( fn,interval = c(alpha_to_nu(p,n,0.001),  alpha_to_nu(p,n,0.999)), maximum = TRUE )$maximum
    alpha <- nu_to_alpha(p,n,nu)
  }
  
  else{
    
    nu <- alpha_to_nu(p = p, n = n , alpha = alpha)
  
  }
  
  delta <- beta_to_delta(p = p,n = n, nu = nu , beta = beta)
  
  
  HMFGraph_GEM_MAP <-  HMFGraph_Gem_algorithm_cpp(iters = max_iters, S = cov(data) ,
                        B= B,  p = p, n = n, stop_criterion = stop_criterion,
                        delta = delta, nu = nu,  inter = inter, 
                        epsilon1 = epsilon1, epsilon2 = epsilon2,
                        fixed_B = fixed_B, print_t = print_t)  
  
  
  return(list(omega = HMFGraph_GEM_MAP$omega, B_i = HMFGraph_GEM_MAP$B_i, 
            phi = HMFGraph_GEM_MAP$phi , varmat= HMFGraph_GEM_MAP$varmat , 
            n = n, p = p , alpha = alpha, beta = beta,
           max_iters = max_iters, stop_criterion = stop_criterion, 
           epsilon1 = epsilon1,epsilon2 = epsilon2, B = B, 
           fixed_B = fixed_B, inter=inter,  print_t=print_t))
}

