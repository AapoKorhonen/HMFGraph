#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <list>
#include <tuple>
#include <boost/math/tools/minima.hpp>
#include <boost/math/special_functions/beta.hpp>

using boost::math::tools::brent_find_minima;
using boost::math::ibetac;

double alphaToDelta(double alpha, int n, int p){
  return (alpha*n+(1-alpha)*p+(1-alpha))/(1-alpha);
}

double deltaToAlpha(double delta, int n, int p){
  return (delta-p-1)/(n+delta-p-1);
}


double BetatoDelta(double beta, double nu,  int n, int p){
  return (beta*nu-p+1-2*beta)/(1-beta);
}


double DeltatoBeta(double delta, double nu,  int n, int p){
  return (delta+p-1)/(delta+nu-2);
}


using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List HMFGraph_Gem_algorithm_cpp(int iters, arma::mat S, const arma::mat B, int p, int n
                                  , double stop_criterion, double delta, double nu, int inter, double epsilon1
                                  , double epsilon2, bool fixed_B,bool print_t, const arma::mat omega_0
                                  ) {
  
  arma::mat B_i = B;
  
  double alpha = deltaToAlpha( nu,  n,  p);
  
  double beta = DeltatoBeta(delta, nu, n,  p);
  
  arma::mat phi_new = B;
  
  arma::mat omega_new = omega_0;
  
  arma::mat omega_old = arma::eye(p, p);
  
  arma::mat varmat = arma::eye(p, p);
  
  double l1_k = (nu + delta +p - 1 - p - 1);
  
  double l2_k = (n + nu - p - 1);
  
  arma::mat L2 = arma::eye(p, p);
  
  arma::mat L1 = arma::eye(p, p);
  
  double shape;
  double rate;
  double ero;
  double suhde;
  
  for(int i = 0; i < iters; i++) {
     
     
     if(!fixed_B){
       
       for(int ii = 0; ii < p; ii++) {
         
         shape = (delta + p - 1) * 0.5 + epsilon1 ;
         
         rate = ((delta +p- 1))* ( (phi_new(ii, ii)) * 0.5 ) + epsilon2;
         
         B_i(ii, ii) =  (shape - 1) / (rate)  ;
       } 
       
       
     }
     
    L1 = B_i*((delta +p- 1)) + omega_new*(nu-p- 1);
    
    phi_new = l1_k * arma::inv(L1);
      

    L2 = arma::inv(phi_new*(nu-p- 1) + n*S);
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

  for(int i = 0; i < (p-1); i++){
    for(int ii = i+1; ii < p; ii++){
      varmat(i,ii) =  (n + nu )*(   pow((L2(i,ii) ),2) + L2(i,i)*L2(ii,ii));
      varmat(ii,i) = varmat(i,ii);
    }
  }

  
  return List::create(
    Named("omega") = omega_new,
    Named("B_i") = B_i,
    Named("phi") = phi_new,
    Named("varmat") = varmat,
    Named("alpha") = alpha,
    Named("beta") = beta
  );
}
