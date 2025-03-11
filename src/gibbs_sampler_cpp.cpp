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
  //arma::mat B_i;
  arma::cube  B_i(p, p, iters + 1, arma::fill::zeros);
  
  Omega.slice(0).eye();  // diag(p)
  B_i.slice(0) = B;
  
  //B_i = B;
  Phi = B;
  arma::mat Omega10 = Omega0;
  
  for (int i = 0; i < iters; ++i) {
    
    progress_bar.increment();
    
    arma::mat L1 = arma::inv(  (delta + p - 1)*B_i.slice(i) + (nu - p - 1)*Omega.slice(i));
    
    // Sampling the Phi matrix
    arma::mat Phi = arma::wishrnd ( L1, nu + delta + p - 1);
    
    //
    arma::mat L2 = arma::inv( (nu - p - 1)*Phi + n*S);
    
    // Sampling the precision matrix
    arma::mat Omega10 = arma::wishrnd(L2,n + nu);
    
    Omega.slice(i + 1) = Omega10;
    
    
    if (!fixed_B) {
      for (int ii = 0; ii < p; ++ii) {
        
        //shape
        double shape = (delta + p - 1)*0.5 + epsilon1;
        
        //rate
        double rate = (delta +p- 1)*(Phi(ii,ii) *0.5 ) +  epsilon2 ;
        
        // Sampling B-matrix. Using shape and scale parametrization 
        //B_i(ii, ii) = R::rgamma( shape ,  1 / rate  );  
        double B_d_ii = R::rgamma( shape ,  1 / rate  );
        arma::mat B_ii = B_i.slice(i + 1);
        B_ii(ii,ii) = B_d_ii;
        B_i.slice(i + 1) = B_ii;
      }
    }
    else{
      B_i.slice(i + 1) = B_i.slice(i);
    }
  }
  
  return List::create(Named("phi") = Phi,
                      Named("omega") = Omega,
                      Named("B") = B_i,
                      Named("delta") = delta,
                      Named("nu") = nu);
}
