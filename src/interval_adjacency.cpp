#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix interval_adjacency(NumericMatrix lower, NumericMatrix upper, int k) {
  
  NumericMatrix adjcency(k);
  
  for (int i = 0; i < k; i++) {
    
    for (int ii = 0; ii < i; ii++) {
      
      double val1 = lower(i, ii);
      
      double val2 = upper(i, ii);
      
      if ((val1 - val2) * (val1 + val2) > 0) {
        
        adjcency(i,ii) = 1;
        
        adjcency(ii,i) = 1;
        
      }
      
    }
    
  }
  
  
  
  return adjcency;
}
