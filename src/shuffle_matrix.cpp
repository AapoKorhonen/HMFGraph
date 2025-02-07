#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <ctime>

// [[Rcpp::export]]
Rcpp::NumericMatrix shuffle_matrix(const Rcpp::NumericMatrix data, int seed_number) {
  
  
  Rcpp::NumericMatrix data100=Rcpp::clone(data);
  int n = data100.nrow();
  int p = data100.ncol();
  
  std::srand(seed_number); // Seed for random number generator
  
  for (int i = 0; i < n; ++i) {
    // Shuffle rows within the i-th row
    std::vector<double> row(p);
    for (int j = 0; j < p; ++j) {
      row[j] = data100(i, j);
    }
    std::random_shuffle(row.begin(), row.end());
    for (int j = 0; j < p; ++j) {
      data100(i, j) = row[j];
    }
  }
  
  return data100;
}


Rcpp::NumericMatrix shuffle_matrix_rows(const Rcpp::NumericMatrix data, int seed_number) {
  
  Rcpp::NumericMatrix data100=Rcpp::clone(data);
  int n = data100.nrow();
  int p = data100.ncol();
  
  std::srand(seed_number); // Seed for random number generator
  
  for (int i = 0; i < n; ++i) {
    // Shuffle rows within the i-th row
    std::vector<double> row(p);
    for (int j = 0; j < p; ++j) {
      row[j] = data100(i, j);
    }
    std::random_shuffle(row.begin(), row.end());
    for (int j = 0; j < p; ++j) {
      data100(i, j) = row[j];
    }
  }
  
  return data100;
}


Rcpp::NumericMatrix shuffle_matrix_cols(const Rcpp::NumericMatrix data, int seed_number) {
  
  Rcpp::NumericMatrix data100=Rcpp::clone(data);
  int n = data100.nrow();
  int p = data100.ncol();
  
  std::srand(seed_number); // Seed for random number generator
  
  for (int i = 0; i < n; ++i) {
    // Shuffle rows within the i-th row
    std::vector<double> col(n);
    for (int j = 0; j < p; ++j) {
      col[j] = data100(j, i);
    }
    std::random_shuffle(col.begin(), col.end());
    for (int j = 0; j < p; ++j) {
      data100(j, i) = col[j];
    }
  }
  
  return data100;
}