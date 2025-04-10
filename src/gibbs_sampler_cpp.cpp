#include <RcppArmadillo.h>
using namespace Rcpp;
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// Gibbs-algoritmi Rcpp:llä
// [[Rcpp::export]]
List gibbs_algorithm_cpp(int iters, const arma::mat& S, arma::mat B, int p, int n, 
                     double delta, double nu, double epsilon1, double epsilon2,  bool fixed_B = false, bool print_t = true,
                     double a_lim = 0.0001, double b_lim =10000 , arma::mat Omega0 = NULL) {
  
  
  
  Progress progress_bar(iters, print_t);
  
  arma::cube Omega(p, p, iters + 1, arma::fill::zeros);
  arma::mat Phi;
  arma::mat B_i;
  
  Omega.slice(0).eye();  // diag(p)
  
  B_i = B;
  Phi = B;
  arma::mat Omega10 = Omega0;
  
  double shape = 1; 
  double rate = 1; 
  
  arma::mat L1;
  arma::mat L2;
  
  for (int i = 0; i < iters; ++i) {
    
    progress_bar.increment();
    
    if (!fixed_B) {
      for (int ii = 0; ii < p; ++ii) {
        
        //shape
        shape = (delta + p - 1)*0.5 + epsilon1;
        
        //rate
        rate = (delta +p- 1)*(Phi(ii,ii) *0.5 ) +  epsilon2 ;
        
        // Sampling B-matrix. Using shape and scale parametrization 
        B_i(ii, ii) = R::rgamma( shape ,  1 / rate  );  
      }
    }
    
    L1 = arma::inv(  (delta + p - 1)*B_i + (nu - p - 1)*Omega.slice(i));
    
    // Sampling the Phi matrix
    Phi = arma::wishrnd ( L1, nu + delta + p - 1);
    
    //
    L2 = arma::inv( (nu - p - 1)*Phi + n*S);
    
    // Sampling the precision matrix
    Omega10 = arma::wishrnd(L2,n + nu);
    
    Omega.slice(i + 1) = Omega10;
    
  }
  
  return List::create(Named("phi") = Phi,
                      Named("omega") = Omega,
                      Named("B") = B_i,
                      Named("delta") = delta,
                      Named("nu") = nu);
}
