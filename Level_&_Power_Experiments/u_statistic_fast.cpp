#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double u_statistic_rcpp(NumericVector x) {
  int n = x.size();
  double sum_val = 0.0;
  
  for (int i = 0; i < n - 2; i++) {
    for (int j = i + 1; j < n - 1; j++) {
      for (int k = j + 1; k < n; k++) {
        double x1 = x[i];
        double x2 = x[j];
        double x3 = x[k];
        double kernel = (16.0/3.0) * x1 * x2 * x3 * (x1 + x2 + x3);
        sum_val += kernel;
      }
    }
  }
  
  double denom = (double)n * (n - 1) * (n - 2) / 6.0; // choose(n,3)
  return sum_val / denom;
}