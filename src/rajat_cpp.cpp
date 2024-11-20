#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rajat_Cpp(NumericMatrix mat1, NumericMatrix mat2, int k, NumericVector lvt) {
  NumericVector yht(lvt.length());
  
  for (int j = 0; j < lvt.length(); j++) {
    double lvt_j = lvt[j];
    for (int i = 0; i < k; i++) {
      for (int ii = 0; ii < i; ii++) {
        double val1 = mat1(i, ii);
        double val2 = mat2(i, ii) * lvt_j;
        if ((val1 - val2) * (val1 + val2) > 0) {
          yht[j] += 1;
        }
      }
    }
  }
  return yht;
}
