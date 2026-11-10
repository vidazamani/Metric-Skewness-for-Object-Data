
#include "RcppArmadillo.h"
#include <vector>

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
// Compute squared affine-invariant distance
// d_AI^2(A,B) = ||log(A^(-1/2) B A^(-1/2))||_F^2
// ------------------------------------------------------------

double distance_ai_squared_cpp(
    const arma::mat& A_inv_sqrt,
    const arma::mat& B
) {
  
  // Compute A^(-1/2) B A^(-1/2)
  arma::mat C = A_inv_sqrt * B * A_inv_sqrt;
  
  // Correct small numerical asymmetry
  C = 0.5 * (C + C.t());
  
  // Eigenvalues of the transformed matrix
  arma::vec eigval;
  
  if (!arma::eig_sym(eigval, C)) {
    stop("Eigen decomposition failed.");
  }
  
  if (arma::any(eigval <= 0)) {
    stop("Transformed matrix is not SPD.");
  }
  
  // Squared affine-invariant distance
  return arma::accu(arma::square(arma::log(eigval)));
}


// ------------------------------------------------------------
// Pairwise squared affine-invariant distances
// ------------------------------------------------------------

//' Compute pairwise squared affine-invariant distances
 //'
 //' @param mats A three-dimensional array containing SPD matrices.
 //' @return A symmetric matrix of pairwise squared AIRM distances.
 //' @export
 // [[Rcpp::export]]
 NumericMatrix affinv_metric_cpp(const arma::cube& mats) {
   
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
       
       double d2 = distance_ai_squared_cpp(
         inv_sqrts[i],
                  mats.slice(j)
       );
       
       D(i, j) = d2;
       D(j, i) = d2;
     }
   }
   
   return D;
 }


// ------------------------------------------------------------
// Pairwise squared distances between A_i and A_j^(-1)
// ------------------------------------------------------------

//' Compute pairwise squared affine-invariant distances
 //' between original and inverted SPD matrices
 //'
 //' @param mats A three-dimensional array containing SPD matrices.
 //' @return A matrix of squared distances between A_i and A_j^(-1).
 //' @export
 // [[Rcpp::export]]
 NumericMatrix affinv_metric_inverse_cpp(const arma::cube& mats) {
   
   const int n = mats.n_slices;
   
   // Precompute inverse square roots and matrix inverses
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
       
       G(i, j) = distance_ai_squared_cpp(
         inv_sqrts[i],
                  inverses[j]
       );
     }
   }
   
   return G;
 }
