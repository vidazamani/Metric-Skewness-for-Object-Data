# Load necessary libraries
library(CovTools) 
library(shapes)
library(ICtest)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)
library(purrr)
library(Rcpp)
library(CompQuadForm)
library(RcppHungarian)

#####
remotes::install_github('vidazamani/Metric-Skewness-for-Object-Data/LogDis-R-Package@V3')
library(LogDis)
# ## OR 
# sourceCpp("/Users/vizama/Documents/Papers/2nd paper/Codes/u2_statistic_rcpp.cpp")
# sourceCpp("/Users/vizama/Documents/Papers/2nd paper/Codes/spd_distances.cpp")
# ## OR
# devtools::install('/Users/vizama/Documents/Papers/2nd paper/Codes/LogDis')
# library(LogDis)
#####



### This is a way to generate data (corr and Cov)

generate_matrices <- function(sample_size,dim,mu,sig) {
  
  replicate(sample_size, {
    
    U <- rorth(dim)
    D <- diag(exp(rnorm(dim,mu,sig)))
    S <- U%*%D%*%t(U)
    
    ## if you want to generate corr matrices, then run the code lines below
    # diagS <- diag(S)
    # diag(1/sqrt(diagS))%*%S%*%diag(1/sqrt(diagS))
    
  })
  
}




# Step 2: Compute average squared distance using distcov
# average_squared_distances <- function(matrices) {
#   n <- dim(matrices)[3]
#   sapply(1:n, function(i) {
#     mean(sapply(1:n,
#                 function(j) if (i != j) distcov(matrices[, , i],
#                                                 matrices[, , j],
#                                                 method ='LogEuclidean')^2 else 0))
#   })
# }




# OR
# rowMeans(distance_logeuclid_cpp(matrices))



# Step 3: Invert all generated matrices
invert_matrices <- function(matrices) {
  n <- dim(matrices)[3] # Number of matrices
  dim <- dim(matrices)[1] # Dimension of each matrix
  inverted_matrices <- array(0, dim = c(dim, dim, n)) # Create an empty 3D array
  for (i in 1:n) {
    inverted_matrices[, , i] <- solve(matrices[, , i]) # Invert each matrix
  }
  return(inverted_matrices)
}

# Step 4: Compute average squared distance for inverted matrices
# average_squared_distances_inverted <- function(matrices) {
#   iverted_matrices <- invert_matrices(matrices)
#   n <- dim(matrices)[3]
#   sapply(1:n, function(i) {
#     mean(sapply(1:n, function(j) distcov(matrices[, , i], iverted_matrices[, , j],
#                                          method ='LogEuclidean')^2))
#   })
# }

# OR
# rowMeans(distance_to_inverse_logeuclid_cpp(matrices))


# --------------------------------------------------------
# Function to compute skewness
# --------------------------------------------------------

metric_skew_fun <- function(matrices) {
  
  # Step 2: Compute average squared distance for original matrices
  avg_dist_original <- rowMeans(distance_logeuclid_cpp(matrices))
  # avg_dist_original <- average_squared_distances(matrices)
  
  
  # Step 3: Compute average squared distance for inverted matrices
  avg_dist_inverted <- rowMeans(distance_to_inverse_logeuclid_cpp(matrices))
  # avg_dist_inverted <- average_squared_distances_inverted(matrices)
  
  
  # Step 4: Compute the final average squared distance
  final_avg_distance <- mean((avg_dist_original - avg_dist_inverted)^2)/mean(avg_dist_original^2)
  
  # Return results
  return(Metric_skewness = final_avg_distance)
}


############ Permutation test ##################################

Perm_test <- function(matrices, iter, regularize){
  
  
  
  # Compute the observed Metric Skewness
  T0 <- metric_skew_fun(matrices)
  
  n <- dim(matrices)[3]
  Tperm <- numeric(iter)
  
  for(b in 1:iter) {
    # random flips: TRUE = invert, FALSE = keep
    flips <- sample(c(TRUE, FALSE), n, replace = TRUE)
    perm_sample <- array(0, dim = dim(matrices))
    for(i in 1:n) {
      if(flips[i]) {
        A <- matrices[,,i]
        if(regularize > 0) A <- A + regularize * diag(dim(generate_matrices(sample_size, dim, mu, sig))[1])
        perm_sample[,,i] <- solve(A)
      } else {
        perm_sample[,,i] <- matrices[,,i]
      }
    }
    Tperm[b] <- metric_skew_fun(perm_sample)
  }
  
  
  
  
  pval <- mean(abs(Tperm) >= abs(T0))
  
  result <-list(statistic = T0, p.value = pval, Tperm = Tperm)
  names(result) <- c("observed statistic","p_value", "permutation stats")
  result
}

matrices <- generate_matrices(5,3,0,1)
Perm_test(matrices,20,0)




############## Asymptotic test for Metric Skewness ###########################


##### Distances

## Using distcove function

# distance_matrix <- function(mats) {
#   n <- dim(mats)[3]
#   D <- matrix(0, n, n)
#   
#   for (i in 1:(n - 1)) {
#     for (j in (i + 1):n) {
#       d2 <- distcov(
#         mats[, , i],
#         mats[, , j],
#         method = "LogEuclidean"
#       )^2
#       
#       D[i, j] <- d2
#       D[j, i] <- d2
#     }
#   }
#   
#   D
# }
# 
# distance_matrix_to_inverse <- function(mats) {
#   inv_mats <- invert_matrices(mats)
#   n <- dim(mats)[3]
#   
#   Dinv <- matrix(0, n, n)
#   
#   for (i in 1:n) {
#     Xi <- mats[, , i]
#     for (j in i:n) {
#       d2 <- distcov(Xi, inv_mats[, , j],
#                     method = "LogEuclidean")^2
#       
#       Dinv[i, j] <- d2
#       Dinv[j, i] <- d2
#     }
#   }
#   
#   Dinv
# }





#### OR
## computing LogEuc Distance manually


# distance_matrix_logeuclid <- function(mats) {
#   # mats: p x p x n array of SPD matrices
#   n <- dim(mats)[3]
#   logs <- lapply(seq_len(n), function(i) logm_spd(mats[,,i]))
#   D <- matrix(0, n, n)
#   for (i in 1:(n-1)) {
#     Li <- logs[[i]]
#     for (j in (i+1):n) {
#       diff <- Li - logs[[j]]
#       d2 <- sum(diff * diff)  # Frobenius norm squared
#       D[i, j] <- d2
#       D[j, i] <- d2
#     }
#   }
#   D
# }
# 
# distance_matrix_to_inverse_logeuclid <- function(mats) {
#   # g(X)=X^{-1}; log(X^{-1}) = -log(X)
#   n <- dim(mats)[3]
#   logs <- lapply(seq_len(n), function(i) logm_spd(mats[,,i]))
#   Dinv <- matrix(0, n, n)
#   for (i in 1:n) {
#     Li <- logs[[i]]
#     for (j in 1:n) {
#       diff <- Li + logs[[j]]  # Li - (-Lj)
#       Dinv[i, j] <- sum(diff * diff)
#     }
#   }
#   Dinv
# }
# 
# D1 <- distance_matrix(mats)
# D2 <- distance_matrix_to_inverse(mats)
# D3 <- distance_logeuclid_cpp(mats)
# D4 <- distance_to_inverse_logeuclid_cpp(mats)
# round(max(abs(D2 - D4)),5)
# max(abs(D1 - D3))


###### U statistics


u1_statistic <- function(D) {
  n <- nrow(D)
  total <- 0
  
  for (i in 1:(n - 2)) {
    for (j in (i + 1):(n - 1)) {
      Dij <- D[i, j]
      
      for (k in (j + 1):n) {
        Dik <- D[i, k]
        Djk <- D[j, k]
        
        total <- total +
          Dij * Dik +
          Djk * Dij +
          Dik * Djk
      }
    }
  }
  
  total / (3 * choose(n, 3))
}



u2_statistic <- function(D, G) {
  n <- nrow(D)
  total <- 0
  
  if (n < 3) {
    stop("U2 statistic requires at least n >= 3 observations.")
  }
  
  stopifnot(
    is.matrix(D), is.matrix(G),
    nrow(D) == ncol(D),
    nrow(G) == ncol(G),
    nrow(D) == nrow(G)
  )
  
  total <- 0
  
  
  for (i in 1:(n - 2)) {
    for (j in (i + 1):(n - 1)) {
      Dij <- D[i, j] - G[i, j]
      
      for (k in (j + 1):n) {
        Dik <- D[i, k] - G[i, k]
        Djk <- D[j, k] - G[j, k]
        
        total <- total +
          Dij * Dik +   # r(i, j, k)
          Djk * Dij +   # r(j, k, i)
          Dik * Djk     # r(k, i, j)
      }
    }
  }
  
  total / (3 * choose(n, 3))
}

# #### Examples 
# 
# set.seed(1)
# 
# n   <- 10
# dim <- 5
# mu  <- 1
# sig <- 0.5
# 
# mats <- generate_matrices(n, dim, mu, sig)
# 
# 
# 
# D <- distance_logeuclid_cpp(mats)
# G <- distance_to_inverse_logeuclid_cpp(mats)
# 
# 
# U2_cpp <- u2_statistic_rcpp(D, G)
# U2_r   <- u2_statistic(D, G)
# 
# all.equal(U2_cpp, U2_r)





################################################



vec0 <- function(S) {
  # S: p x p symmetric matrix
  p <- nrow(S)
  out <- numeric(p * (p + 1) / 2)
  idx <- 1L
  for (j in 1:p) {
    for (i in 1:j) {
      if (i == j) {
        out[idx] <- S[i, j]
      } else {
        out[idx] <- sqrt(2) * S[i, j]
      }
      idx <- idx + 1L
    }
  }
  out
}




##### Sigma hat

logm_spd <- function(X) {
  # X: p x p SPD matrix
  ee <- eigen(X, symmetric = TRUE)
  Q  <- ee$vectors
  d  <- ee$values
  if (any(d <= 0)) stop("Matrix not SPD (non-positive eigenvalues).")
  Q %*% diag(log(d), length(d)) %*% t(Q)
}


cov_vec0_log <- function(mats) {
  # mats: p x p x n
  n <- dim(mats)[3]
  V <- lapply(seq_len(n), function(i) vec0(logm_spd(mats[,,i])))
  V <- do.call(rbind, V)  # n x q, q = p(p+1)/2
  cov(V)                  # sample covariance (q x q)
}


imhof_cdf <- function(x, weights) {
  res <- CompQuadForm::imhof(q = 1, 
                             lambda = weights/x,
                             epsabs = 1e-15, 
                             epsrel = 1e-15)
  
  res$Qq
}



#################################################################

#### Asymptotic Test

Asymp_metric_skewness_spd <- function(mats) {
  # mats: p x p x n SPD array
  n <- dim(mats)[3]
  if (n < 3) return(list(statistic = NA_real_, p.value = NA_real_))
  
  # 1) Log–Euclidean distances (it is actullay Frobenius norm squared)
  D <- distance_logeuclid_cpp(mats)
  G <- distance_to_inverse_logeuclid_cpp(mats)
  
  
  # 2) U2 from Rcpp
  U2 <- u2_statistic_rcpp(D, G)
  
  # 3) Sigma from vec0(log(X))
  Sigma_hat <- cov_vec0_log(mats)
  
  
  
  # 4) Eigenvalues and trace(Sigma^2)
  eigvals <- eigen(Sigma_hat, symmetric = TRUE, only.values = TRUE)$values
  theta_sq <- eigvals^2
  trSigma2 <- sum(theta_sq)  # = tr(Sigma_hat^2)
  
  # 5) Test statistic (correct parentheses)
  Tn <- (n * U2) / 16 + trSigma2
  
  # # 6) Upper-tail p-value under weighted chi-square limit
  cdf_val <- imhof_cdf(Tn, theta_sq)
  pval <- cdf_val
  
  
  # # 6) CDF via Imhof
  # cdf_val <- imhof_cdf(Tn, theta_sq)
  # 
  # # two-sided p-value
  # pval <- 2 * min(cdf_val, 1 - cdf_val)
  # 
  # # keep within [0,1] numerically
  # pval <- min(max(pval, 0), 1)
  
  list(
    statistic = Tn,
    p.value = pval,
    U2 = U2,
    trSigma2 = trSigma2,
    eigenvalues = eigvals
  )
}


Asymp_metric_skewness_spd(mats)


#################### Wasserstein-2 bootstrap-permutation test for skewness #########
# Assumes that n is even
#
# D = squared distance matrix
# G = squared distance matrix from original to flipped
# R = number of bootstrap reps
# B = number of premutation reps
wasserstein_test <- function(mats, R, B){
  pval_vec <- rep(0, R)
  
  n <- dim(mats)[3]
  
  # Distance matrices
  D <- distance_logeuclid_cpp(mats)
  G <- distance_to_inverse_logeuclid_cpp(mats)
  
  
  # Bootstrap replicates
  for(i in 1:R){
    b_indices <- sample(1:n)
    
    D0 <- D[b_indices, b_indices]
    G0 <- G[b_indices, b_indices]
    
    G_corner <- G0[1:(n/2), (n/2 + 1):n]
    
    cost0 <- HungarianSolver(G_corner)$cost
    
    Dhn <- rbind(cbind(D0[1:(n/2), 1:(n/2)], G0[1:(n/2), (n/2 + 1):n]),
                 cbind(t(G0[1:(n/2), (n/2 + 1):n]), D0[(n/2 + 1):n, (n/2 + 1):n]))
    
    perm_costs <- rep(0, B)
    
    # Permutation replicates
    for(j in 1:B){
      p_indices <- sample(1:n)
      Dhn0 <- Dhn[p_indices, p_indices]
      Dhn0_corner <- Dhn0[1:(n/2), (n/2 + 1):n]
      
      perm_costs[j] <- HungarianSolver(Dhn0_corner)$cost
    }
    
    pval_vec[i] <- mean(cost0 < perm_costs)
  }
  
  mean(pval_vec)
}



# Wasserstein-2 permutation test for skewness
# Single-rep version suggested by Janne
#
# Assumes that n is even
# B = number of permutation reps

wasserstein_test_2 <- function(mats, B){
  
  n <- dim(mats)[3]
  
  
  # Distance matrices
  D <- distance_logeuclid_cpp(mats)
  G <- distance_to_inverse_logeuclid_cpp(mats)
  
  G_corner <- G[1:(n/2), (n/2 + 1):n]
  
  cost0 <- HungarianSolver(G_corner)$cost
  
  Dhn <- rbind(cbind(D[1:(n/2), 1:(n/2)], G[1:(n/2), (n/2 + 1):n]),
               cbind(t(G[1:(n/2), (n/2 + 1):n]), D[(n/2 + 1):n, (n/2 + 1):n]))
  
  perm_costs <- rep(0, B)
  
  # Permutation replicates
  for(j in 1:B){
    p_indices <- sample(1:n)
    Dhn0 <- Dhn[p_indices, p_indices]
    Dhn0_corner <- Dhn0[1:(n/2), (n/2 + 1):n]
    
    perm_costs[j] <- HungarianSolver(Dhn0_corner)$cost
  }
  
  mean(cost0 <= perm_costs)
}




################################## Power Evaluation ########################

### Mac

power_fixed_n <- function(
    n,
    dim, 
    mu, 
    sig,
    nrep,
    R,
    B,
    alpha
) {
  
  
  
  ncores = detectCores() - 1
  
  
  power <- matrix(NA, nrow = length(mu), ncol = 4)
  colnames(power) <- c("Metric_perm", "Metric_asymp","Wasserstein1", "Wasserstein2")
  
  for (k in seq_along(mu)) {
    
    mu_skew <- mu[k]
    
    
    pvals_list <- mclapply(seq_len(nrep), function(r) {
      
      
      mats <- generate_matrices(n, dim, mu_skew, sig)
      
      
      c(
        Metric_perm  = Perm_test(mats, B,0)$p_value,
        Metric_asymp = Asymp_metric_skewness_spd(mats)$p.value,
        Wasserstein1  = wasserstein_test(mats, R, B),
        Wasserstein2  = wasserstein_test_2(mats, B)
      )
      
    }, mc.cores = ncores)
    
    
    
    # convert list → matrix
    pvals <- do.call(cbind, pvals_list)
    
    power[k, ] <- rowMeans(pvals < alpha, na.rm = TRUE)
    
    cat("n =", n,
        "mu =", paste(mu_skew, collapse = ","),
        "power =", round(power[k, ], 3), "\n")
  }
  
  data.frame(
    mu,
    power
  )
}

#######


sig_grid <- c(0.05, 0.1, 0.3)
n_values <- c(20, 50, 100, 200)
mu <- seq(0,0.06,0.01)
dim <- 3
nrep <- 100
R <- 50
B <- 50
alpha <- 0.05

run_power_all_sigma <- function(sig_grid, n_values){
  
  results_list <- list()
  
  for(sig in sig_grid){
    
    for(n in n_values){
      
      cat("Running for sigma =", sig, "n =", n, "\n")
      
      res <- power_fixed_n(n , dim, mu, sig, nrep, R, B, alpha)
      
      df <- res |>
        pivot_longer(
          cols = -mu,
          names_to = "Test",
          values_to = "Power"
        ) |>
        mutate(
          n = n,
          sigma = sig
        )
      
      results_list[[paste(n, sig, sep = "_")]] <- df
    }
  }
  
  bind_rows(results_list)
}



df_power <- run_power_all_sigma(sig_grid, n_values)

df_power$Test <- factor(df_power$Test,
                        levels = c("Metric_perm", "Metric_asymp", "Wasserstein1", "Wasserstein2"),
                        labels = c("Permutation", "Asymptotic", "Wasserstein1", "Wasserstein2"))

df_power$n <- factor(df_power$n)
df_power$sigma <- factor(df_power$sigma)



p <- ggplot(
  df_power,
  aes(x = mu,
      y = Power,
      color = Test,
      group = Test)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  facet_grid(sigma ~ n,
             labeller = labeller(
               n = function(x) paste0("n = ", x),
               sigma = function(x) paste0("sigma = " , x)
             )) +
  scale_color_manual(values = c("#0072B2", "#D55E00","#E69F00","#56B4E9")) +
  scale_y_continuous(
    limits = c(-0.03, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = expression(mu),
    y = "Power",
    color = "Test"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth = 0.4)
  )


ggsave("/Users/vizama/Documents/Papers/2nd paper/Simulation results/pics/Power Evaluation Cov Data.pdf",
       plot = p, width = 12, height = 7)





############# Level evaluation - Simulation across sample sizes to estimate level

## Mac
level_test_parallel <- function(sample_sizes, 
                                dim, 
                                mu, 
                                sig,
                                nrep,
                                R,
                                B, 
                                alpha) {
  
  
  
  # for reproducibility across forked processes
  RNGkind("L'Ecuyer-CMRG")
  
  ncores = detectCores() - 1
  
  
  results_hhat         <- numeric(length(sample_sizes))
  results_hhat_asym    <- numeric(length(sample_sizes))
  results_Wasserstein1  <- numeric(length(sample_sizes))
  results_Wasserstein2  <- numeric(length(sample_sizes))
  
  
  for (i in seq_along(sample_sizes)) {
    n <- sample_sizes[i]
    
    
    pvals_list <- mclapply(seq_len(nrep), function(r) {
      
      
      mats <- generate_matrices(n, dim, mu, sig)
      
    
      
      c(
        Metric_perm  = Perm_test(mats, B,0)$p_value,
        Metric_asymp = Asymp_metric_skewness_spd(mats)$p.value,
        Wasserstein1  = wasserstein_test(mats, R, B),
        Wasserstein2  = wasserstein_test_2(mats, B)
      )
      
    }, mc.cores = ncores)
    
    
    # convert list → matrix
    pvals <- do.call(cbind, pvals_list)
    
    
    
    results_hhat[i]      <- mean(pvals[1, ] < alpha, na.rm = TRUE)
    results_hhat_asym[i] <- mean(pvals[2, ] < alpha, na.rm = TRUE)
    results_Wasserstein1[i] <- mean(pvals[3, ] < alpha, na.rm = TRUE)
    results_Wasserstein2[i] <- mean(pvals[4, ] < alpha, na.rm = TRUE)
    
    cat(
      "Done n =", n,
      "| metric perm:", round(results_hhat[i], 3),
      "| metric asym:", round(results_hhat_asym[i], 3),
      "| Wasserstein1:", round(results_Wasserstein1[i], 3),
      "| Wasserstein2:", round(results_Wasserstein2[i], 3), "\n")
    
  }
  
  
  
  list(sample_sizes = sample_sizes,
       metric_perm = results_hhat,
       metric_asym = results_hhat_asym,
       Wasserstein1 = results_Wasserstein1,
       Wasserstein2 = results_Wasserstein2)
}






set.seed(1111)

dim  <- 3
mu   <- 0
sig  <- 0.1
alpha <- 0.05
sample_sizes <- seq(10,300,20)

nrep <- 100
R <- 50
B <- 50

level_cov <- level_test_parallel(
  sample_sizes,
  dim,
  mu, 
  sig,
  nrep,
  R,
  B,
  alpha
)






# -----------------------------
# Plot results
# -----------------------------
df_spd <- tibble(
  n = level_cov$sample_sizes,
  Metric_perm = level_cov$metric_perm,
  Metric_asym = level_cov$metric_asym,
  Wasserstein1 = level_cov$Wasserstein1,
  Wasserstein2 = level_cov$Wasserstein2
) |>
  pivot_longer(-n, names_to = "Statistic", values_to = "Rejection") 





df_spd$Statistic <- factor(df_spd$Statistic)

p2 <- ggplot(df_spd,
       aes(x = n,
           y = Rejection,
           color = Statistic,
           shape = Statistic,
           group = Statistic)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.05,
             linetype = "dashed",
             color = "black",
             linewidth = 0.5) +
  scale_color_manual(values = c("#009E73", "#CC79A7","#E69F00","#56B4E9")) +
  scale_shape_manual(values = c(18, 15, 8, 16)) +
  coord_cartesian(ylim = c(0, 0.15)) +
  labs(
    x = "Sample size",
    y = expression("Rejection proportion  (" * alpha == 0.05 * ")"),
    color = "Statistic",
    shape = "Statistic"
  ) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 15),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth = 0.5)
  )


ggsave("/Users/vizama/Documents/Papers/2nd paper/Simulation results/pics/matrix/Level Evaluation Cov Data 2.pdf",
       plot = p2, width = 10, height = 7)
