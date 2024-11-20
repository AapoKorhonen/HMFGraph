#' Generates data with a scale free network structure.
#'
#' @param n The number of samples
#' @param p The number of variables
#' @param d The number of connected variables. Because the function produces scale-free network with maximum of p-1 connections, this values should be between 1 and p.
#'
#' @return Returns the generated data with precision, covariance and adjacency matrices.
#' @export
#'
#' @examples
#' 
data_generator <- function(n, p, d=0){

  
  if(d==0){
    d = p
  }
  
  # Firstly the graph structure is generated
  
  presicion_matrix <- diag(p)
  
  connected_nodes <- sample(c(1:p))[1:d]
  
  round(d/2)
  
  connections_S <- rep(0, d)
  
  for(i in 1:d){
    if (0.5 < runif(1)){
      connections_S[i] <- runif(1, 0.2, 0.9)
    }
    else{
      connections_S[i] <- runif(1, -0.9, -0.2)
      
      
    }
  }
  
  connections_S <- sample(connections_S)
  
  connections_sums <- rep(0, d)
  
  for(i in 1:d){
    
    if(i == 1){
      presicion_matrix[ connected_nodes[1],connected_nodes[2]] <- connections_S[i]
      presicion_matrix[ connected_nodes[2],connected_nodes[1]] <- connections_S[i]
      
      connections_sums[1] = 1 + connections_sums[1]
      connections_sums[2] = 1 + connections_sums[2]
    }
    
    else{
      node <- sample(connected_nodes[-i], size = 1, replace = F, prob = connections_sums[-i]/i)
      presicion_matrix[ connected_nodes[i],node] <- connections_S[i]
      presicion_matrix[ node,connected_nodes[i]] <- connections_S[i]
      
      connections_sums[i] = 1 + connections_sums[i]
      connections_sums[node == connected_nodes] = 1 + connections_sums[node == connected_nodes]
    }
  }
  
  omega <- presicion_matrix
  
  omega_old <- omega
  
  diag(omega) <- 0
  
  diag(omega) = abs(min(eigen(omega)$values)) + 0.2
  
  sigma = cov2cor(solve(omega))
  
  omega = solve(sigma)
  
  # Using Rcpp for sampling. This is much faster than mvrnorm in R if p >> 100.
  
  x = mvrnorm_cpp(n, rep(0, p), sigma)
  
  adjacency_matrix <- omega_old
  
  adjacency_matrix[abs(omega_old) > 0] <- 1
  
  diag(adjacency_matrix) <- 0
  
  return(list(data = x, precision_matrix = omega, covariance_matrix = sigma, adjacency_matrix=adjacency_matrix))
}

