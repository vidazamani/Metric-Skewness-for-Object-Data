#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Compute matrix logarithm of an SPD matrix via eigen-decomposition
arma::mat logm_spd_cpp(const arma::mat& X) {
  arma::vec eigval;
  arma::mat eigvec;
  if (!arma::eig_sym(eigval, eigvec, X)) {
    stop("Eigen decomposition failed.");
  }
  if (arma::any(eigval <= 0)) {
    stop("Matrix not SPD (non-positive eigenvalues).");
  }
  arma::mat L = eigvec * arma::diagmat(arma::log(eigval)) * eigvec.t();
  return L;
}

// [[Rcpp::export]]
NumericMatrix distance_logeuclid_cpp(const arma::cube& mats) {
  // const int p = mats.n_rows;
  const int n = mats.n_slices;
  
  // Precompute logs
  std::vector<arma::mat> logs(n);
  for (int i = 0; i < n; ++i) {
    logs[i] = logm_spd_cpp(mats.slice(i));
  }
  
  NumericMatrix D(n, n);
  
  for (int i = 0; i < n - 1; ++i) {
    const arma::mat& Li = logs[i];
    for (int j = i + 1; j < n; ++j) {
      const arma::mat& Lj = logs[j];
      arma::mat diff = Li - Lj;
      double d2 = arma::accu(diff % diff); // Frobenius norm squared
      D(i, j) = d2;
      D(j, i) = d2;
    }
  }
  return D;
}

// [[Rcpp::export]]
NumericMatrix distance_to_inverse_logeuclid_cpp(const arma::cube& mats) {
 // const int p = mats.n_rows;
  const int n = mats.n_slices;
  
  // Precompute logs
  std::vector<arma::mat> logs(n);
  for (int i = 0; i < n; ++i) {
    logs[i] = logm_spd_cpp(mats.slice(i));
  }
  
  NumericMatrix G(n, n);
  
  // Use log(X^{-1}) = -log(X)
  for (int i = 0; i < n; ++i) {
    const arma::mat& Li = logs[i];
    for (int j = 0; j < n; ++j) {
      const arma::mat& Lj = logs[j];
      arma::mat diff = Li + Lj; // Li - (-Lj)
      double d2 = arma::accu(diff % diff);
      G(i, j) = d2;
    }
  }
  return G;
}
