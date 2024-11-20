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



// 
// 
// #include <Rcpp.h>
// #include <algorithm>
// #include <vector>
// #include <cstdlib>
// #include <ctime>
// 
// // [[Rcpp::export]]
// Rcpp::NumericMatrix shuffle_matrix(const Rcpp::NumericMatrix data, int seed_number) {
//   
//   Rcpp::NumericMatrix data100=Rcpp::clone(data);
//   int n = data100.nrow();
//   int p = data100.ncol();
//   
//   std::srand(seed_number); // Seed for random number generator
//   
//   // for (int i = 0; i < p; ++i) {
//   //   // Shuffle columns in the i-th column
//   //   std::vector<double> column(n);
//   //   for (int j = 0; j < n; ++j) {
//   //     column[j] = data100(j, i);
//   //   }
//   //   std::random_shuffle(column.begin(), column.end());
//   //   for (int j = 0; j < n; ++j) {
//   //     data100(j, i) = column[j];
//   //   }
//   // }
//   
//   for (int i = 0; i < n; ++i) {
//     // Shuffle rows within the i-th row
//     std::vector<double> row(p);
//     for (int j = 0; j < p; ++j) {
//       row[j] = data100(i, j);
//     }
//     std::random_shuffle(row.begin(), row.end());
//     for (int j = 0; j < p; ++j) {
//       data100(i, j) = row[j];
//     }
//   }
//   
//   return data100;
// }
