
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;


// ------------------------------------------------------------
// Symmetric eigenvalue decomposition using Jacobi rotations
// ------------------------------------------------------------

void jacobi_eigen(
    std::vector<double>& A,
    std::vector<double>& eigenvalues,
    std::vector<double>& eigenvectors,
    int n
) {

  eigenvectors.assign(n * n, 0.0);

  for (int i = 0; i < n; ++i) {
    eigenvectors[i * n + i] = 1.0;
  }

  const int max_iter = 100 * n * n;
  const double tol = 1e-12;

  for (int iter = 0; iter < max_iter; ++iter) {

    int p = 0;
    int q = 1;
    double max_val = 0.0;

    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {

        double val = std::abs(A[i * n + j]);

        if (val > max_val) {
          max_val = val;
          p = i;
          q = j;
        }
      }
    }

    if (max_val < tol) {
      break;
    }

    double app = A[p * n + p];
    double aqq = A[q * n + q];
    double apq = A[p * n + q];

    double tau = (aqq - app) / (2.0 * apq);

    double t = (tau >= 0.0 ? 1.0 : -1.0) /
      (std::abs(tau) + std::sqrt(1.0 + tau * tau));

    double c = 1.0 / std::sqrt(1.0 + t * t);
    double s = t * c;

    // Update diagonal elements
    A[p * n + p] = app - t * apq;
    A[q * n + q] = aqq + t * apq;

    A[p * n + q] = 0.0;
    A[q * n + p] = 0.0;

    // Update remaining matrix elements
    for (int k = 0; k < n; ++k) {

      if (k != p && k != q) {

        double akp = A[k * n + p];
        double akq = A[k * n + q];

        A[k * n + p] = c * akp - s * akq;
        A[p * n + k] = A[k * n + p];

        A[k * n + q] = s * akp + c * akq;
        A[q * n + k] = A[k * n + q];
      }
    }

    // Update eigenvectors
    for (int k = 0; k < n; ++k) {

      double vkp = eigenvectors[k * n + p];
      double vkq = eigenvectors[k * n + q];

      eigenvectors[k * n + p] = c * vkp - s * vkq;
      eigenvectors[k * n + q] = s * vkp + c * vkq;
    }
  }

  eigenvalues.resize(n);

  for (int i = 0; i < n; ++i) {
    eigenvalues[i] = A[i * n + i];
  }
}


// ------------------------------------------------------------
// Extract a symmetric matrix from an R matrix
// ------------------------------------------------------------

std::vector<double> get_symmetric_matrix(
    const NumericMatrix& X
) {

  int n = X.nrow();

  std::vector<double> A(n * n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {

      A[i * n + j] =
        0.5 * (X(i, j) + X(j, i));
    }
  }

  return A;
}


// ------------------------------------------------------------
// Compute squared AIRM distance
// ------------------------------------------------------------

double distance_ai_squared_cpp(
    const NumericMatrix& A,
    const NumericMatrix& B
) {

  int n = A.nrow();

  // Eigen-decomposition of A
  std::vector<double> eigA;
  std::vector<double> vecA;

  std::vector<double> symA = get_symmetric_matrix(A);

  jacobi_eigen(symA, eigA, vecA, n);

  // Check positive definiteness
  for (int i = 0; i < n; ++i) {
    if (eigA[i] <= 0.0) {
      stop("Matrix A is not SPD.");
    }
  }

  // Construct A^(-1/2)
  std::vector<double> inv_sqrt_A(n * n, 0.0);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {

      double value = 0.0;

      for (int k = 0; k < n; ++k) {

        value +=
          vecA[i * n + k] *
          vecA[j * n + k] /
            std::sqrt(eigA[k]);
      }

      inv_sqrt_A[i * n + j] = value;
    }
  }

  // Compute C = A^(-1/2) B A^(-1/2)
  std::vector<double> C(n * n, 0.0);

  // First multiplication: temp = A^(-1/2) B
  std::vector<double> temp(n * n, 0.0);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {

        temp[i * n + j] +=
          inv_sqrt_A[i * n + k] * B(k, j);
      }
    }
  }

  // Second multiplication: C = temp A^(-1/2)
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {

        C[i * n + j] +=
          temp[i * n + k] *
          inv_sqrt_A[k * n + j];
      }
    }
  }

  // Enforce numerical symmetry
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {

      double value =
        0.5 * (C[i * n + j] + C[j * n + i]);

      C[i * n + j] = value;
      C[j * n + i] = value;
    }
  }

  // Eigen-decomposition of C
  std::vector<double> eigC;
  std::vector<double> vecC;

  jacobi_eigen(C, eigC, vecC, n);

  // Calculate squared AIRM distance
  double distance = 0.0;

  for (int i = 0; i < n; ++i) {

    if (eigC[i] <= 0.0) {
      stop("Transformed matrix is not SPD.");
    }

    distance += std::pow(std::log(eigC[i]), 2.0);
  }

  return distance;
}


// ------------------------------------------------------------
// Pairwise squared AIRM distances
// ------------------------------------------------------------

//' @export
 // [[Rcpp::export]]
 NumericMatrix affinv_metric_cpp(NumericVector mats) {

   IntegerVector dims = mats.attr("dim");

   if (dims.size() != 3) {
     stop("Input must be a 3D array of matrices.");
   }

   int p = dims[0];
   int q = dims[1];
   int n = dims[2];

   if (p != q) {
     stop("Matrices must be square.");
   }

   NumericMatrix D(n, n);

   for (int i = 0; i < n; ++i) {

     NumericMatrix A(p, p);

     for (int r = 0; r < p; ++r) {
       for (int c = 0; c < p; ++c) {
         A(r, c) = mats[r + p * c + p * p * i];
       }
     }

     for (int j = i + 1; j < n; ++j) {

       NumericMatrix B(p, p);

       for (int r = 0; r < p; ++r) {
         for (int c = 0; c < p; ++c) {
           B(r, c) = mats[r + p * c + p * p * j];
         }
       }

       double d2 = distance_ai_squared_cpp(A, B);

       D(i, j) = d2;
       D(j, i) = d2;
     }
   }

   return D;
 }


// ------------------------------------------------------------
// Pairwise squared AIRM distances between A_i and A_j^(-1)
// ------------------------------------------------------------

//' @export
 // [[Rcpp::export]]
 NumericMatrix affinv_metric_inverse_cpp(NumericVector mats) {

   IntegerVector dims = mats.attr("dim");

   if (dims.size() != 3) {
     stop("Input must be a 3D array of matrices.");
   }

   int p = dims[0];
   int q = dims[1];
   int n = dims[2];

   if (p != q) {
     stop("Matrices must be square.");
   }

   NumericMatrix G(n, n);

   // Construct all matrix inverses
   std::vector<NumericMatrix> inverses(n);

   for (int i = 0; i < n; ++i) {

     NumericMatrix A(p, p);

     for (int r = 0; r < p; ++r) {
       for (int c = 0; c < p; ++c) {
         A(r, c) = mats[r + p * c + p * p * i];
       }
     }

     // Gauss-Jordan inversion
     NumericMatrix augmented(p, 2 * p);

     for (int r = 0; r < p; ++r) {
       for (int c = 0; c < p; ++c) {
         augmented(r, c) = A(r, c);
         augmented(r, c + p) = (r == c ? 1.0 : 0.0);
       }
     }

     for (int k = 0; k < p; ++k) {

       int pivot = k;

       for (int r = k + 1; r < p; ++r) {
         if (std::abs(augmented(r, k)) >
               std::abs(augmented(pivot, k))) {
           pivot = r;
         }
       }

       if (std::abs(augmented(pivot, k)) < 1e-14) {
         stop("Matrix inversion failed.");
       }

       if (pivot != k) {
         for (int c = 0; c < 2 * p; ++c) {
           std::swap(augmented(k, c), augmented(pivot, c));
         }
       }

       double divisor = augmented(k, k);

       for (int c = 0; c < 2 * p; ++c) {
         augmented(k, c) /= divisor;
       }

       for (int r = 0; r < p; ++r) {

         if (r != k) {

           double factor = augmented(r, k);

           for (int c = 0; c < 2 * p; ++c) {
             augmented(r, c) -=
               factor * augmented(k, c);
           }
         }
       }
     }

     NumericMatrix invA(p, p);

     for (int r = 0; r < p; ++r) {
       for (int c = 0; c < p; ++c) {
         invA(r, c) = augmented(r, c + p);
       }
     }

     inverses[i] = invA;
   }

   // Calculate distances between A_i and inverse(A_j)
   for (int i = 0; i < n; ++i) {

     NumericMatrix A(p, p);

     for (int r = 0; r < p; ++r) {
       for (int c = 0; c < p; ++c) {
         A(r, c) = mats[r + p * c + p * p * i];
       }
     }

     for (int j = 0; j < n; ++j) {

       G(i, j) = distance_ai_squared_cpp(
         A,
         inverses[j]
       );
     }
   }

   return G;
 }
