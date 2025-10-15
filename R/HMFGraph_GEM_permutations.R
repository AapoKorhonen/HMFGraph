
#' HMFGraph GEM algorithm with permutated data.
#'
#' @param data The original data set 
#' @param HMFGraph_GEM_RESULTS Results from the function HMFGraph_GEM
#' @param number_of_permutations The number of permutations to be run
#' @param parallel If TRUE, then parallel computing is used.
#' @param seed A seed for permutation function.
#' @param n_cores The number of cores for the permutation calculations
#' @param median_p The median of permutation is used, otherwise uses mean.
#' @param omega_0 The initial value for Omega for the algorithm.
#'
#' @return Returns two vectors and vector on quantile points. The first vector consists the number of estimated connections for each quantile point with the original data and the second consist median number of connections for permutated data.
#' @export
#'
#' @examples library(HMFGraph)
#' @examples n <- 200
#' @examples p <- 100
#' @examples set.seed(42)
#' @examples generated_data <- data_generator(n=n, p = p)
#' @examples results_HMFGraph_GEM <- HMFGraph_GEM(generated_data$data, alpha = p * 5 / ( p * 5+n), beta=0.9)
#' @examples permutations <- HMFGraph_GEM_permutations(generated_data$data, results_HMFGraph_GEM, number_of_permutations = 100, parallel = F)
HMFGraph_GEM_permutations <- function(data, HMFGraph_GEM_RESULTS, number_of_permutations = 50, parallel = TRUE, seed = FALSE, n_cores = 0, median_p=T, omega_0 = NULL){
  
  
  
  if(is.null(omega_0)){
    omega_0 = diag(HMFGraph_GEM_RESULTS$p)
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
    
    # pb <- utils::txtProgressBar(min = 0, max = number_of_permutations, style = 3)
    # 
    # progress <- function(n) utils::setTxtProgressBar(pb, n)
    
    pb <- progress::progress_bar$new(format = "Permutations :percent [:bar] :elapsed | eta: :eta",
                                     total = number_of_permutations+1, width = 80)
    
    progress <- function() pb$tick()
    
    pb$tick()
    
    opts <- list(progress = progress)
    
    `%dopar%` <- foreach::`%dopar%`
    
    permutation_results <- foreach::foreach(i = 1:n_clusters, .combine = 'rbind', .options.snow = opts)  %dopar%  {
      
      data_mixed <- shuffle_matrix(data, seed_number = seed + i)
      
      HMFGraph_GEM_results_permutation <- HMFGraph_GEM(data_mixed, n = HMFGraph_GEM_RESULTS$n, p = HMFGraph_GEM_RESULTS$p, 
              alpha = HMFGraph_GEM_RESULTS$alpha, beta = HMFGraph_GEM_RESULTS$beta,
              max_iters = HMFGraph_GEM_RESULTS$max_iters,
              stop_criterion = HMFGraph_GEM_RESULTS$stop_criterion, 
              epsilon1 = HMFGraph_GEM_RESULTS$epsilon1,
              epsilon2 = HMFGraph_GEM_RESULTS$epsilon2, B = HMFGraph_GEM_RESULTS$B, 
              fixed_B = HMFGraph_GEM_RESULTS$fixed_B, inter=HMFGraph_GEM_RESULTS$inter,
              print_t=F, omega_0 = omega_0 )
      
      
      
      qp_connections <- rajat_Cpp(HMFGraph_GEM_results_permutation$omega,sqrt(HMFGraph_GEM_results_permutation$varmat),HMFGraph_GEM_RESULTS$p,quantile_points)
      
      qp_connections
    }
    
    
    
    #base::close(pb)
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
      
      HMFGraph_GEM_results_permutation <- HMFGraph_GEM(data_mixed, n = HMFGraph_GEM_RESULTS$n, p = HMFGraph_GEM_RESULTS$p, 
              alpha = HMFGraph_GEM_RESULTS$alpha, beta = HMFGraph_GEM_RESULTS$beta,
              max_iters = HMFGraph_GEM_RESULTS$max_iters,
              stop_criterion = HMFGraph_GEM_RESULTS$stop_criterion, 
              epsilon1 = HMFGraph_GEM_RESULTS$epsilon1,
              epsilon2 = HMFGraph_GEM_RESULTS$epsilon2, B = HMFGraph_GEM_RESULTS$B, 
              fixed_B = HMFGraph_GEM_RESULTS$fixed_B, inter=HMFGraph_GEM_RESULTS$inter,
              print_t=F, omega_0 = omega_0)
      
      qp_connections_permutation[i,] <- rajat_Cpp(HMFGraph_GEM_results_permutation$omega,sqrt(HMFGraph_GEM_results_permutation$varmat),HMFGraph_GEM_RESULTS$p,quantile_points)
      
      
    }
    if(median_p){
      qp_connections_permutation<-  apply(qp_connections_permutation,2,median)
      
    }
    else{
      qp_connections_permutation<-  apply(qp_connections_permutation,2,mean)
      
    }
    
  }
  
  
  cat("Permutations finished", "\n")
  
  qp_connections <- rajat_Cpp(HMFGraph_GEM_RESULTS$omega,sqrt(HMFGraph_GEM_RESULTS$varmat),HMFGraph_GEM_RESULTS$p,quantile_points)
  
  return(list(qp_connections = qp_connections, qp_connections_permutation = qp_connections_permutation, quantile_points=quantile_points))
}
