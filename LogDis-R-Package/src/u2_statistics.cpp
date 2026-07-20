// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double u2_statistic_rcpp(const NumericMatrix& D,
                         const NumericMatrix& G) {
  
  int n = D.nrow();
  
  if (n < 3)
    stop("U2 statistic requires at least n >= 3 observations.");
  
  if (D.ncol() != n || G.nrow() != n || G.ncol() != n)
    stop("D and G must be square matrices of the same size.");
  
  double total = 0.0;
  
  for (int i = 0; i < n - 2; ++i) {
    for (int j = i + 1; j < n - 1; ++j) {
      
      double Dij = D(i, j) - G(i, j);
      
      for (int k = j + 1; k < n; ++k) {
        
        double Dik = D(i, k) - G(i, k);
        double Djk = D(j, k) - G(j, k);
        
        total +=
          Dij * Dik +   // r(i, j, k)
          Djk * Dij +   // r(j, k, i)
          Dik * Djk;    // r(k, i, j)
      }
    }
  }
  
  // 3 * choose(n, 3) = n * (n-1) * (n-2) / 2
  double denom = (double)n * (n - 1) * (n - 2) / 2.0;
  
  return total / denom;
}
