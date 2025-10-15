
#' HMFGraph, with GEM algorithm, permutations and an optimal CI.
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
#' @param B The initial value for the target matrix B. If fixed_B is TRUE, then this value will not change for the whole sampling process.
#' @param omega_0 The initial value for Omega. On default it is an identity matrix. A correctly selected initial value can speed up the algorithm, e.i. warm start.
#' @param fixed_B Boolean value. If TRUE, then no prior for target matrix B is used.
#' @param inter Interval of iterations when info about algorithm is printed out.
#' @param print_t If FALSE, then no text is printed during the function running process.
#' @param number_of_permutations The number of permutations to be run
#' @param parallel If TRUE, then parallel computing is used.
#' @param seed A seed for permutation function.
#' @param n_cores The number of cores for the permutation calculations
#' @param median_p The median of permutation is used, otherwise uses mean.
#' @param omega_0 The initial value for Omega for the algorithm.
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
#' @examples results <- HMFGraph(generated_data$data)
#' @examples library(qgraph)
#' @examples qgraph(results$adjacency_matrix)
#' 
HMFGraph <- function(data, p = NULL,  n = NULL, 
                         alpha = -1,  beta = 0.9, max_iters = 10000,
                         stop_criterion = 10^(-6)  ,epsilon1 = 0.001, epsilon2 = 0.001,
                         B = NULL,fixed_B = F, inter=500,  print_t=T,
                         kappa_max=NULL, omega_0 = NULL,
                         max_steps=50, threshold=0.05, alpha_selection=F, 
                         lower_alpha = p/(p+n), print_binary_search = F,print_t_alpha=F, 
                         number_of_permutations = 50, parallel = TRUE, seed = FALSE, 
                         n_cores = 0, median_p=T, 
                         expected_connections= NULL, MCC=F){
  
  
  # If n is not specified, the number of samples are derived from the dimension of the data matrix
  if(is.null(n)){
    n = dim(data)[1]
  }
  
  if(is.null(p)){
    p = dim(data)[2]
  }
  
  if(is.null(B)){
    B = diag(p)
  }
  if(is.null(omega_0)){
    omega_0 = diag(p)
  }
  
  
  if(alpha == -1 || alpha_selection){
    
    cat("Selecting an optimal alpha value", "\n")
    
    alpha <- alpha_binary_search(data= data, p = p,  n = n, 
                                 alpha = alpha,  beta = beta, max_iters = max_iters,
                                 stop_criterion = stop_criterion  ,epsilon1 = epsilon1, 
                                 epsilon2 = epsilon2, B = B,fixed_B = fixed_B, inter=inter, print_t=print_t_alpha,
                                 omega_0 = omega_0, max_steps=max_steps, threshold=threshold,
                                 lower_alpha = lower_alpha, print_binary_search=print_binary_search)
    
    cat("The optimal alpha: " ,alpha, "\n")
    
  }
  
  nu <- alpha_to_nu(p = p, n = n , alpha = alpha)
  
  delta <- beta_to_delta(p = p,n = n, nu = nu , beta = beta)
  
  
  HMFGraph_GEM_MAP <-  HMFGraph_Gem_algorithm_cpp(iters = max_iters, S = cov(data) ,
                                                  B= B,  p = p, n = n, stop_criterion = stop_criterion,
                                                  delta = delta, nu = nu,  inter = inter, 
                                                  epsilon1 = epsilon1, epsilon2 = epsilon2,
                                                  fixed_B = fixed_B, print_t = print_t,
                                                  omega_0 = omega_0)  
  

  if(is.null(omega_0)){
    omega_0 = diag(p)
  }
  
  quantile_points <- seq(0.00001,7, length.out=5000)
  
  if(seed==F){
    seed <- sample(.Random.seed,1)
  }
  
  if(parallel){
    
    if(n_cores <= 0){
      n_cores <- parallel::detectCores() - 1
    }
    
    n_clusters <- number_of_permutations
    
    cat("Starting permutations", "\n")
    
    cl <- parallel::makeCluster(min(n_cores, n_clusters) , type = "SOCK")
    
    doSNOW::registerDoSNOW(cl)
    
    pb <- progress::progress_bar$new(format = "Permutations :percent [:bar] :elapsed | eta: :eta",
                                     total = number_of_permutations+1, width = 80)
    
    progress <- function() pb$tick()
    
    pb$tick()
    
    opts <- list(progress = progress)
    
    `%dopar%` <- foreach::`%dopar%`
    
    permutation_results <- foreach::foreach(i = 1:n_clusters, .combine = 'rbind', .options.snow = opts)  %dopar%  {
      
      data_mixed <- shuffle_matrix(data, seed_number = seed + i)
      
      HMFGraph_GEM_results_permutation <- HMFGraph_GEM(data_mixed, n = n, p = p, 
                                                       alpha = alpha, beta = beta,
                                                       max_iters = max_iters,
                                                       stop_criterion = stop_criterion, 
                                                       epsilon1 = epsilon1,
                                                       epsilon2 = epsilon2, B = B, 
                                                       fixed_B = fixed_B, inter=inter,
                                                       print_t=F, omega_0 = omega_0 )
      
      
      
      qp_connections <- rajat_Cpp(HMFGraph_GEM_results_permutation$omega,sqrt(HMFGraph_GEM_results_permutation$varmat),p,quantile_points)
      
      qp_connections
    }
    
    
    parallel::stopCluster(cl)
    
    if(median_p){
      qp_connections_permutation<-  apply(permutation_results,2,median)
      
    }
    else{
      qp_connections_permutation<-  apply(permutation_results,2,mean)
      
    }
    
    
  }
  else{
    
    pb <- progress::progress_bar$new(format = "Permutations :percent [:bar] :elapsed | eta: :eta",
                                     total = number_of_permutations, width = 80)
    
    cat("Starting permutations", "\n")
    
    qp_connections_permutation <- matrix(0,nrow= number_of_permutations ,ncol=length(quantile_points))
    
    for(i in 1:number_of_permutations){
      
      pb$tick()
      
      data_mixed <- shuffle_matrix(data, seed_number = seed + i)
      
      HMFGraph_GEM_results_permutation <- HMFGraph_GEM(data_mixed, n = n, p = p, 
                                                       alpha = alpha, beta = beta,
                                                       max_iters = max_iters,
                                                       stop_criterion = stop_criterion, 
                                                       epsilon1 = epsilon1,
                                                       epsilon2 = epsilon2, B = B, 
                                                       fixed_B = fixed_B, inter=inter,
                                                       print_t=F, omega_0 = omega_0 )
      
      
      qp_connections_permutation[i,] <- rajat_Cpp(HMFGraph_GEM_results_permutation$omega,sqrt(HMFGraph_GEM_results_permutation$varmat),p,quantile_points)
      
      
    }
    if(median_p){
      qp_connections_permutation<-  apply(qp_connections_permutation,2,median)
      
    }
    else{
      qp_connections_permutation<-  apply(qp_connections_permutation,2,mean)
      
    }
    
  }
  
  
  cat("Permutations finished", "\n")
  
  qp_connections <- rajat_Cpp(HMFGraph_GEM_MAP$omega,sqrt(HMFGraph_GEM_MAP$varmat),p,quantile_points)
  
  
  
  if(is.null( expected_connections )){
    expected_connections <- p
  }
  TP <- qp_connections - qp_connections_permutation 
  
  FP <- qp_connections_permutation
  
  FN <- expected_connections - TP
  
  F1 <- (2*TP)/(2*TP+FP+FN)
  
  N <- (p^2/2 - p)/2 - expected_connections
  
  TN <- N - qp_connections_permutation
  
  F1[is.nan(F1)] <- 0
  F1[is.na(F1)] <- 0
  
  if(MCC){
    
    MCC <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    
    MCC[is.na(MCC)] <- 0
    F1 <- MCC
    
  }
  
  
  quantile_point <- quantile_points[F1==max(F1[!is.nan(F1)])][1]
  
  MAP_estimate <- HMFGraph_GEM_MAP$omega 
  
  variance_matrix <- HMFGraph_GEM_MAP$varmat 
  
  
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
              upper_CI = upper_CI, F1=max(F1[!is.nan(F1)]), quantile_point =quantile_point, alpha=alpha,beta=beta))  
}
