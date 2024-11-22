

#' Changes alpha parameter value to nu
#'
#' @param p The number of variables
#' @param n  The number of samples
#' @param alpha Value for alpha parameter
#'
#' @return returns nu
alpha_to_nu  <- function(p,n, alpha){
  
  nu <- (alpha*n-alpha*p-alpha+p+1)/(1-alpha)
  
  return(nu)
  
} 


#' Changes nu parameter value to alpha
#'
#' @param p The number of variables
#' @param n  The number of samples 
#' @param nu Value for nu parameter
#'
#' @return returns alpha
nu_to_alpha  <- function(p,n, nu){
  
  alpha <- (nu-p-1)/(nu+n-p-1)
  
  return(alpha)
  
} 


#' Changes beta parameter value to delta
#'
#' @param p The number of variables
#' @param n The number of samples 
#' @param nu Value for nu parameter
#' @param beta Value for beta parameter
#'
#' @return returns delta
beta_to_delta  <- function(p,n, nu,beta){
  
  delta <- (beta*nu-p+1-2*beta)/(1-beta)
  
  return(delta)
}


#' Changes delta  parameter value to beta
#'
#' @param p The number of variables
#' @param n The number of samples 
#' @param nu Value for nu parameter
#' @param delta Value for delta parameter
#'
#' @return returns beta
delta_to_beta  <- function(p,n, nu,delta){
  
  beta <- (delta+p-1)/(delta+nu-2)
  
  return(beta)
  
} 


# 
# 
# logmarginal_determinant <- function(x,n,p,D,S){
#   
#   p6 <- (x/2)*determinant(D*(x-p-1),logarithm = T)$modulus[1]
#   
#   p5 <- ((x+n)/2)*determinant( D*(x-p-1) + S ,logarithm = T)$modulus[1]
#   
#   p1 <- ((n*p)/2)*log(pi)
#   
#   p2 <- (CholWishart::lmvgamma((x+n)/2, p)) 
#   
#   p3 <- (CholWishart::lmvgamma((x)/2, p))
#   
#   
#   pp <- -p1 +p2 -p3 - p5 + p6 
#   
#   return(pp)
# }


#' Calculates the log-marginal likelihood 
#'
#' @param x Values for the parameter nu
#' @param n The number of samples 
#' @param p The number of variables
#' @param eigen_D Eigenvalues of matrix D
#' @param eigen_inv_DS Eiganvalues of inverse of matrix product of S and D
#'
#' @return Returns the value of the log-marginal likelihood
logmarginal <- function(x,n,p,eigen_D,eigen_inv_DS){
  
  logdetD <- (sum(log(eigen_D)))
  
  p1 <- ((n*p)/2)*log(pi)
  
  p2 <- (CholWishart::lmvgamma((x+n)/2, p)) 
  
  p3 <- (CholWishart::lmvgamma((x)/2, p))
  
  p4 <- (x/2)*p*log(x-p-1)
  
  p5 <- ((x+n)/2)*sum(log(x-p-1+eigen_inv_DS))
  
  p6 <- (n)/2*logdetD
  
  pp <- -p1 +p2 -p3 + p4 - p5 - p6 
  
  
  return(pp)
}


