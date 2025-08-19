#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <list>
#include <tuple>
#include <boost/math/tools/minima.hpp>
#include <boost/math/special_functions/beta.hpp>

using boost::math::tools::brent_find_minima;
using boost::math::ibetac;

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List Matrix_F_Gem_algorithm_cpp(int iters, arma::mat S, const arma::mat B, int p, int n
                                  , double stop_criterion, double delta, double nu, int inter, 
                                  bool print_t, const arma::mat omega_0
) {
  
  arma::mat B_i = B;
  
  arma::mat phi_new = B;
  
  arma::mat omega_new = omega_0;
  
  arma::mat omega_old = arma::eye(p, p);
  
  arma::mat varmat = arma::eye(p, p);
  
  double l1_k = (nu + delta +p - 1);
  
  double l2_k = (n + nu - p - 1);
  
  arma::mat L2 = arma::eye(p, p);
  
  arma::mat L1 = arma::eye(p, p);
  
  double ero;
  double suhde;
  
  for(int i = 0; i < iters; i++) {
    
    
    L1 = B_i + omega_new;
    
    phi_new = l1_k * arma::inv(L1);
    
    
    L2 = arma::inv(phi_new + n*S);
    omega_new =  l2_k * L2 ;
    
    
    ero = norm(omega_new - omega_old, "fro");
    suhde = ero / norm(omega_old, "fro");
    omega_old = omega_new;
    if(i % inter == 0) {
      if(print_t){
        Rcout << "Iteration: " << i << ", Relative Difference: " << suhde << std::endl;
      }
    }
    
    if(suhde < stop_criterion) {
      if(i > 2) {
        if(print_t){
          Rcout << "Convergence reached at iteration: " << i << ", Relative Difference: " << suhde << std::endl;
        }
        break;
      }
    }
  }
  
  
  return List::create(
    Named("omega") = omega_new,
    Named("B_i") = B_i,
    Named("phi") = phi_new
  );
}
