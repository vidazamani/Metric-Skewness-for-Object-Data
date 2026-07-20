#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix distance_matrix_mv_cpp(const NumericMatrix& X) {
  
  int n = X.nrow();
  int d = X.ncol();
  
  NumericMatrix D(n, n);
  
  for(int i = 0; i < n - 1; i++) {
    for(int j = i + 1; j < n; j++) {
      
      double d2 = 0.0;
      
      for(int k = 0; k < d; k++) {
        double diff = X(i, k) - X(j, k);
        d2 += diff * diff;
      }
      
      D(i, j) = d2;
      D(j, i) = d2;
    }
  }
  
  return D;
}

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix distance_matrix_to_flip_cpp(const NumericMatrix& X) {
  
  int n = X.nrow();
  int d = X.ncol();
  
  NumericMatrix Dflip(n, n);
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      
      double d2 = 0.0;
      
      for(int k = 0; k < d; k++) {
        double sumv = X(i, k) + X(j, k);
        d2 += sumv * sumv;
      }
      
      Dflip(i, j) = d2;
    }
  }
  
  return Dflip;
}