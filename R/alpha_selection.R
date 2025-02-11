#' Selecting suitable alpha value based on a condition number constraint method.
#'
#' @param data The data 
#'
#' @return
#' @export
#'
#' @examples

alpha_binary_search <- function(data, p = NULL,  n = NULL, 
                                alpha = -1,  beta = 0.9, max_iters = 10000,
                                stop_criterion = 10^(-6)  ,epsilon1 = 0, epsilon2 = 0,
                                B = diag(p),fixed_B = F, inter=500,  print_t=T,
                                kappa_max=NULL, omega_0 = diag(p),
                                max_steps=50, threshold=0.05, lower_alpha = p/(p+n),
                                print_binary_search=F){
  

  if(is.null(n)){
    n = dim(data)[1]
  }
  
  if(is.null(p)){
    p = dim(data)[2]
  }
  
  if(is.null(kappa_max)){
    
    est <- nlshrink::linshrink_cov(data)
    kappa_max <- max(eigen(est)$values)/min(eigen(est)$values)
    
  }
  
  if(is.null(lower_alpha)){
    lower_alpha = p/(n+p)
  }
  
  lower_alpha <- lower_alpha
  upper_alpha <- 0.999
  
  alpha_p <- (lower_alpha+upper_alpha)/2
  
  omega_0 = diag(p)
  
  if(print_binary_search){
    cat("Initial alpha: " , alpha_p, "\n")
  }
  
  
  eigen_max <- 0
  eigen_min <- 0
  
  alpha <- 0.999
  
  if(print_binary_search){
    cat("Kappa max: " , kappa_max, "\n")
  }
  
  for(i in c(1:max_steps)){
    
    nu <- alpha_to_nu(p = p, n = n , alpha = alpha_p)
    delta <- max(beta_to_delta(p = p,n = n, nu = nu , beta = beta),1)
    
    HMFGraph_GEM_MAP <-  HMFGraph_Gem_algorithm_cpp(iters = max_iters, S = cov(data) ,
                                                    B = B,  p = p, n = n, stop_criterion = stop_criterion,
                                                    delta = delta, nu = nu,  inter = inter, 
                                                    epsilon1 = epsilon1, epsilon2 = epsilon2,
                                                    fixed_B = fixed_B, print_t = F,
                                                    omega_0 = omega_0)  
    
    
    
    B <- HMFGraph_GEM_MAP$B_i
    omega_0 <- HMFGraph_GEM_MAP$omega
    eigen_min <-min(eigen(HMFGraph_GEM_MAP$omega)$values ) 
    
    eigen_max <-max(eigen(HMFGraph_GEM_MAP$omega)$values ) 
    
    err <- (eigen_max/eigen_min - kappa_max) / kappa_max
    
    if(print_binary_search){
      cat("Condition number: " ,eigen_max/eigen_min, "\n")
      cat("Relative difference: " ,abs(err), "\n")
      
    }
    
    
    if(err < 0){
      if(abs(err) < threshold){
        alpha <- alpha_p
        break
      }
    }
    if(err > 0){
      lower_alpha <- alpha_p
      alpha_p <- (lower_alpha+upper_alpha)/2
    }
    else{
      upper_alpha <- alpha_p
      alpha_p <- (lower_alpha+upper_alpha)/2
    }
    if(alpha_p == 0.999){
      alpha <- alpha_p
      break
    }
    
    if(print_binary_search){
      cat("New alpha: " ,alpha_p, "\n")
    }
    
  }
  
  if(print_binary_search){
    cat("Number of iterarions needed for binary search: ", i, "\n")
    cat("Selected alpha: " ,alpha, "\n")
  }
  return(alpha)
  
}