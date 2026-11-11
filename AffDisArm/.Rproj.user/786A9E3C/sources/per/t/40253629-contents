
#include "RcppArmadillo.h"
#include <vector>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// ------------------------------------------------------------
// Compute inverse square root of an SPD matrix
// ------------------------------------------------------------

arma::mat inv_sqrt_spd_cpp(const arma::mat& X) {
  
  arma::vec eigval;
  arma::mat eigvec;
  
  if (!arma::eig_sym(eigval, eigvec, X)) {
    stop("Eigen decomposition failed.");
  }
  
  if (arma::any(eigval <= 0)) {
    stop("Matrix is not SPD (non-positive eigenvalues).");
  }
  
  return eigvec *
    arma::diagmat(1.0 / arma::sqrt(eigval)) *
    eigvec.t();
}


// ------------------------------------------------------------
// Compute ordinary affine-invariant Riemannian distance
//
// d_AI(A,B) =
// sqrt(sum(log(lambda_k)^2))
//
// where lambda_k are eigenvalues of
// A^(-1/2) B A^(-1/2)
// ------------------------------------------------------------

double distance_ai_cpp(
    const arma::mat& A_inv_sqrt,
    const arma::mat& B
) {
  
  // Compute A^(-1/2) B A^(-1/2)
  arma::mat C = A_inv_sqrt * B * A_inv_sqrt;
  
  // Correct small numerical asymmetry
  C = 0.5 * (C + C.t());
  
  // Eigenvalues of transformed matrix
  arma::vec eigval;
  
  if (!arma::eig_sym(eigval, C)) {
    stop("Eigen decomposition failed.");
  }
  
  if (arma::any(eigval <= 0)) {
    stop("Transformed matrix is not SPD.");
  }
  
  // Ordinary AIRM distance
  double d2 = arma::accu(
    arma::square(arma::log(eigval))
  );
  
  return std::sqrt(d2);
}


// ------------------------------------------------------------
// Pairwise ordinary affine-invariant distances
// ------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix affinv_metric_cpp(
    const arma::cube& mats
) {
  
  const int n = mats.n_slices;
  
  // Precompute inverse square roots
  std::vector<arma::mat> inv_sqrts(n);
  
  for (int i = 0; i < n; ++i) {
    inv_sqrts[i] = inv_sqrt_spd_cpp(mats.slice(i));
  }
  
  NumericMatrix D(n, n);
  
  // Pairwise distances
  for (int i = 0; i < n - 1; ++i) {
    
    for (int j = i + 1; j < n; ++j) {
      
      double d = distance_ai_cpp(
        inv_sqrts[i],
                 mats.slice(j)
      );
      
      D(i, j) = d;
      D(j, i) = d;
    }
  }
  
  // Diagonal is zero
  for (int i = 0; i < n; ++i) {
    D(i, i) = 0.0;
  }
  
  return D;
}


// ------------------------------------------------------------
// Pairwise ordinary AIRM distances between A_i and inverse(A_j)
// ------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix affinv_metric_inverse_cpp(
    const arma::cube& mats
) {
  
  const int n = mats.n_slices;
  
  // Precompute inverse square roots and inverses
  std::vector<arma::mat> inv_sqrts(n);
  std::vector<arma::mat> inverses(n);
  
  for (int i = 0; i < n; ++i) {
    
    const arma::mat& A = mats.slice(i);
    
    inv_sqrts[i] = inv_sqrt_spd_cpp(A);
    
    if (!arma::inv_sympd(inverses[i], A)) {
      stop("Matrix inversion failed.");
    }
  }
  
  NumericMatrix G(n, n);
  
  // Distances between A_i and inverse(A_j)
  for (int i = 0; i < n; ++i) {
    
    for (int j = 0; j < n; ++j) {
      
      G(i, j) = distance_ai_cpp(
        inv_sqrts[i],
                 inverses[j]
      );
    }
  }
  
  return G;
}
