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

double lpvarGamma(const double x, const int p) {
  double ans = (p * (p - 1) * 0.25) * log(3.14159);
  for(int j = 1; j < (p + 1); ++j){
    ans += std::lgamma(x - ((j - 1.0) * 0.5));
  }
  return ans;
}


double logML(const double nu, const int p, const int n, arma::mat B, arma::mat omega, double gamma1, double gamma2){

  
  std::complex<double> juts = (gamma1-1)*log(nu)-gamma2*nu-nu*p/2*log(2)- arma::trace(  ((nu-p-1)*B)*omega  )/2+ ((nu-p-1)/2)*arma::log_det(omega)+(nu/2)*arma::log_det((nu-p-1)*B)-1*lpvarGamma(nu/2,p);
  
  double kats = juts.real() ;
  
  return(kats);
}


double get_nu(const int p, const int n, arma::mat B, arma::mat Omega, double gamma1, double gamma2){
  
  const double lowerVal = alphaToDelta(0.001, n, p);
  const double upperVal = alphaToDelta(0.999, n, p);
  
  const auto obj = [p, n, B, Omega,  gamma1,  gamma2](double x) { return -logML(x, p, n, B, Omega,  gamma1,  gamma2); };
  
  
  boost::uintmax_t it = 1000;
  
  const auto result = brent_find_minima(obj, lowerVal, upperVal, 1000, it);
  auto nu = 0.0, value = 0.0;
  std::tie(nu, value) = result;
  
  return(nu);
}

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List hwo_Gem_algorithm_cpp(int iters, arma::mat S, const arma::mat B, int p, int n
                                              , double stop_criterion, double delta, double nu, int inter, double epsilon1
                                              , double epsilon2, bool fixed_B,bool print_t, const arma::mat omega_0
                                              , double gamma1, double gamma2
) {
  
  arma::mat B_i = B;
  
  arma::mat omega_new = omega_0;
  
  arma::mat omega_old = arma::eye(p, p);
  
  arma::mat varmat = arma::eye(p, p);
  
  //double l2_k = (n + nu - p - 1);
  
  arma::mat L2 = arma::eye(p, p);
  
  
  double shape;
  double rate;
  double ero;
  double suhde;
  
  for(int i = 0; i < iters; i++) {
    
    
    if(!fixed_B){
      
      for(int ii = 0; ii < p; ii++) {
        
        shape = (nu-p- 1) * 0.5 + epsilon1 ;
        
        rate = ((nu-p- 1))* ( (omega_new(ii, ii)) * 0.5 ) + epsilon2;
        
        B_i(ii, ii) =  (shape) / (rate)  ;
      } 
      
      
    }
    
    
    L2 = arma::inv(B_i*(nu-p- 1) + n*S);
    omega_new =  (n + nu - p - 1) * L2 ;
    
    nu = get_nu(p,n,B_i,omega_new,gamma1, gamma2);
    
    
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
    Named("varmat") = varmat,
    Named("nu") = nu
  );
}
