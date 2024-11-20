

#' Title
#'
#' @param p 
#' @param n 
#' @param alpha 
#'
#' @return
#'
#' @examples
alpha_to_nu  <- function(p,n, alpha){
  
  nu <- (alpha*n-alpha*p-alpha+p+1)/(1-alpha)
  
  return(nu)
  
} 


#' Title
#'
#' @param p 
#' @param n 
#' @param nu 
#'
#' @return
#'
#' @examples
nu_to_alpha  <- function(p,n, nu){
  
  alpha <- (nu-p-1)/(nu+n-p-1)
  
  return(alpha)
  
} 


#' Title
#'
#' @param p 
#' @param n 
#' @param nu 
#' @param beta 
#'
#' @return
#'
#' @examples
beta_to_delta  <- function(p,n, nu,beta){
  
  delta <- (beta*nu-p+1-2*beta)/(1-beta)
  
  return(delta)
}


#' Title
#'
#' @param p 
#' @param n 
#' @param nu 
#' @param delta 
#'
#' @return
#'
#' @examples
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


#' Title
#'
#' @param x 
#' @param n 
#' @param p 
#' @param eigen_D 
#' @param eigen_inv_DS 
#'
#' @return
#'
#' @examples
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


