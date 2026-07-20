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
library(TreeDist)

#####
#remotes::install_github('vidazamani/Metric-Skewness-for-Object-Data/LogDis-R-Package')
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
    
    
    cost0 <- LAPJV(G_corner)[[1]]
    
    
    Dhn <- rbind(cbind(D0[1:(n/2), 1:(n/2)], G0[1:(n/2), (n/2 + 1):n]),
                 cbind(t(G0[1:(n/2), (n/2 + 1):n]), D0[(n/2 + 1):n, (n/2 + 1):n]))
    
    perm_costs <- rep(0, B)
    
    # Permutation replicates
    for(j in 1:B){
      p_indices <- sample(1:n)
      Dhn0 <- Dhn[p_indices, p_indices]
      Dhn0_corner <- Dhn0[1:(n/2), (n/2 + 1):n]
      
      perm_costs[j] <- LAPJV(Dhn0_corner)[[1]]
    }
    
    pval_vec[i] <- mean(cost0 <= perm_costs)
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
  
  cost0 <- LAPJV(G_corner)[[1]]
  
  Dhn <- rbind(cbind(D[1:(n/2), 1:(n/2)], G[1:(n/2), (n/2 + 1):n]),
               cbind(t(G[1:(n/2), (n/2 + 1):n]), D[(n/2 + 1):n, (n/2 + 1):n]))
  
  perm_costs <- rep(0, B)
  
  # Permutation replicates
  for(j in 1:B){
    p_indices <- sample(1:n)
    Dhn0 <- Dhn[p_indices, p_indices]
    Dhn0_corner <- Dhn0[1:(n/2), (n/2 + 1):n]
    
    perm_costs[j] <- LAPJV(Dhn0_corner)[[1]]
  }
  
  mean(cost0 <= perm_costs)
}









  

  
  

  
# ─────────────────────────────────────────────────────────────
# Timing Study
# ─────────────────────────────────────────────────────────────
  
  # Fixed parameters
  dim       <- 3       # matrix dimension
  mu        <- 0       # mean 
  sig       <- 1       # Sigma
  n_rep     <- 50      # number of repetitions per sample size
  iter_perm <- 200     # iterations for Perm_test
  R_wass    <- 50      # bootstrap replicates for wasserstein_test
  B_wass    <- 200     # permutation replicates for wasserstein tests
  reg       <- 0       # regularization parameter for Perm_test
  
  sample_sizes <- c(20, 50, 100, 200)
  
# Storage: one row per sample size, columns = mean & SD for each test
  results <- data.frame(
    sample_size          = sample_sizes,
    perm_mean            = NA_real_,
    perm_sd              = NA_real_,
    asymp_mean           = NA_real_,
    asymp_sd             = NA_real_,
    wass_mean            = NA_real_,
    wass_sd              = NA_real_,
    wass2_mean           = NA_real_,
    wass2_sd             = NA_real_
  )
  
  # ─────────────────────────────────────────────────────────────
  # Main timing loop
  # ─────────────────────────────────────────────────────────────
  
  for (s in seq_along(sample_sizes)) {
    
    n <- sample_sizes[s]
    cat(sprintf("\n── Sample size: %d ──\n", n))
    
    # Pre-generate all replicate datasets to ensure fair comparison
    datasets <- lapply(1:n_rep, function(x) generate_matrices(n, dim, mu, sig))
    
    # ── 1. Perm_test ──────────────────────────────────────────
    cat("  Running Perm_test...\n")
    t_perm <- numeric(n_rep)
    for (r in 1:n_rep) {
      t_perm[r] <- system.time(
        Perm_test(datasets[[r]], iter = iter_perm, regularize = reg)
      )["elapsed"]
    }
    
    # ── 2. Asymp_metric_skewness_spd ──────────────────────────
    cat("  Running Asymp_metric_skewness_spd...\n")
    t_asymp <- numeric(n_rep)
    for (r in 1:n_rep) {
      t_asymp[r] <- system.time(
        Asymp_metric_skewness_spd(datasets[[r]])
      )["elapsed"]
    }
    
    # ── 3. wasserstein_test ───────────────────────────────────
    cat("  Running wasserstein_test...\n")
    t_wass <- numeric(n_rep)
    for (r in 1:n_rep) {
      t_wass[r] <- system.time(
        wasserstein_test(datasets[[r]], R = R_wass, B = B_wass)
      )["elapsed"]
    }
    
    # ── 4. wasserstein_test_2 ─────────────────────────────────
    cat("  Running wasserstein_test_2...\n")
    t_wass2 <- numeric(n_rep)
    for (r in 1:n_rep) {
      t_wass2[r] <- system.time(
        wasserstein_test_2(datasets[[r]], B = B_wass)
      )["elapsed"]
    }
    
    # ── Store results ─────────────────────────────────────────
    results[s, "perm_mean"]  <- mean(t_perm);   results[s, "perm_sd"]  <- sd(t_perm)
    results[s, "asymp_mean"] <- mean(t_asymp);  results[s, "asymp_sd"] <- sd(t_asymp)
    results[s, "wass_mean"]  <- mean(t_wass);   results[s, "wass_sd"]  <- sd(t_wass)
    results[s, "wass2_mean"] <- mean(t_wass2);  results[s, "wass2_sd"] <- sd(t_wass2)
    
  }
  
# ─────────────────────────────────────────────────────────────
# Print summary table
# ─────────────────────────────────────────────────────────────
  
  cat("\n\n══════════════════════════════════════════════════════════════\n")
  cat("Timing Study Results (seconds) — averaged over", n_rep, "replications\n")
  cat("══════════════════════════════════════════════════════════════\n\n")
  
  print(
    data.frame(
      n             = results$sample_size,
      Perm          = sprintf("%.4f (%.4f)", results$perm_mean,  results$perm_sd),
      Asymp         = sprintf("%.4f (%.4f)", results$asymp_mean, results$asymp_sd),
      Wasserstein   = sprintf("%.4f (%.4f)", results$wass_mean,  results$wass_sd),
      Wasserstein_2 = sprintf("%.4f (%.4f)", results$wass2_mean, results$wass2_sd)
    ),
    row.names = FALSE
  )
  
  
library(ggplot2)
library(tidyr)
library(dplyr)
  
# ─────────────────────────────────────────────────────────────
# Reshape results into long format for ggplot2
# ─────────────────────────────────────────────────────────────
  
results_long <- results %>%
    # Means
    select(sample_size, perm_mean, asymp_mean, wass_mean, wass2_mean) %>%
    pivot_longer(
      cols      = -sample_size,
      names_to  = "test",
      values_to = "mean_time"
    ) %>%
    # SDs
    left_join(
      results %>%
        select(sample_size, perm_sd, asymp_sd, wass_sd, wass2_sd) %>%
        pivot_longer(
          cols      = -sample_size,
          names_to  = "test_sd",
          values_to = "sd_time"
        ) %>%
        mutate(test = gsub("_sd", "_mean", test_sd)) %>%
        select(-test_sd),
      by = c("sample_size", "test")
    ) %>%
    # Clean up test labels
    mutate(
      test = recode(test,
                    "perm_mean"  = "Metric Perm",
                    "asymp_mean" = "Metric Asymp",
                    "wass_mean"  = "Wass (Boot+Perm)",
                    "wass2_mean" = "Wass (Perm)"
      ),
      sample_size = factor(sample_size)
    )



results_long <- results_long %>%
  group_by(test) %>%
  mutate(time_norm = mean_time / max(mean_time))
  
# ─────────────────────────────────────────────────────────────
# Plot 1: Line plot with error bands (mean ± SD) 
# ─────────────────────────────────────────────────────────────
  
p1 <- ggplot(results_long, 
               aes(x    = as.numeric(as.character(sample_size)),
                   y    = mean_time,
                   col  = test,
                   fill = test)) +
    geom_ribbon(aes(ymin = mean_time - sd_time,
                    ymax = mean_time + sd_time),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_x_continuous(breaks = c(20, 50, 100, 200)) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette  = "Set1") +
    labs(
      title    = "",
      x        = "Sample Size (n)",
      y        = "Time (seconds)",
      col      = "Test",
      fill     = "Test"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position  = c(0.3,0.8),
      legend.background = element_rect(fill="white"),,
      legend.key.size = unit(0.35, "cm"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.box.background = element_rect(color="black"),
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
# ─────────────────────────────────────────────────────────────
# Plot 2: Faceted bar plot with error bars (mean ± SD)
# ─────────────────────────────────────────────────────────────
  
p2 <- ggplot(results_long,
               aes(x    = sample_size,
                   y    = mean_time,
                   fill = test)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_errorbar(aes(ymin = mean_time - sd_time,
                      ymax = mean_time + sd_time),
                  position = position_dodge(0.7),
                  width    = 0.25,
                  linewidth = 0.6) +
    facet_wrap(~ test, scales = "free_y") +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title    = "",
      x        = "Sample Size (n)",
      y        = "Time (seconds)",
      fill     = "Test"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position  = "none",
      plot.title       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(face = "bold")
    ) 
  
# ─────────────────────────────────────────────────────────────
# Display plots
# ─────────────────────────────────────────────────────────────
  
  print(p1)
  print(p2)
  
  
library(patchwork)
  
p <- (p1 | p2) 


ggsave("/Users/vizama/Documents/Papers/2nd paper/Simulation results/pics/matrix/timing.png", plot = p, width = 8, height = 5, dpi = 300)
ggsave("timing_lineplot.png", plot = p1, width = 8, height = 5, dpi = 300)
ggsave("timing_barplot.png",  plot = p2, width = 9, height = 6, dpi = 300)
