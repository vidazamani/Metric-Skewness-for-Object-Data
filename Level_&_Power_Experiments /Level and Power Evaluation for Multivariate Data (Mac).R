# Load necessary libraries

library(sn)
library(MASS)
library(parallel)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mvnormalTest)
library(Rcpp)
library(CompQuadForm)

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

# --------------------------------------------------------
# Function to compute h-hat(p) (Metric skewness)
# --------------------------------------------------------
metric_skew <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  M <- ncol(X)
  
  ## ---------- Numerator ----------
  # term_j = -(4/n) * sum_i sum_m x_{j m} x_{i m}
  # note: sum_i x_{i m} can be precomputed
  col_sum <- colSums(X)                     # length M
  term_j <- -(4 / n) * (X %*% col_sum)      # n x 1 vector
  
  numerator <- mean(term_j^2)
  
  ## ---------- Denominator ----------
  denom_j <- numeric(n)
  
  for (j in 1:n) {
    diff_sq <- sweep(X, 2, X[j, ], "-")^2   # (x_jm - x_im)^2
    denom_j[j] <- (1 / n) * sum(diff_sq)
  }
  
  denominator <- mean(denom_j^2)
  
  ## ---------- h-hat(p) ----------
  numerator / denominator
}



# --------------------------------------------------------
# Mardia's multivariate skewness
# --------------------------------------------------------
mardia_skewness <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  Xc <- scale(X, center = TRUE, scale = FALSE)
  S  <- cov(Xc)
  S_inv <- solve(S)
  
  # Mahalanobis inner products
  A <- Xc %*% S_inv %*% t(Xc)
  
  # Mardia's skewness
  sum(A^3) / n^2
}


############################################################

# Example: n observations, M variables
set.seed(1)
X <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3)

metric_skew(X)




# Example: multivariate normal 

n <- 200
p <- 3


gen_mvn <- function(n, p) {
  mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
}


metric_skew(gen_mvn(n,p))

# Example: skewed data

##### Azzalini’s skew-normal 

n <- 200              # sample size
p <- 3                # dimension
## With the parameter alpha, you can control the skewness and 
## Direction of skewness controlled by sign of alpha

alpha <- c(1, 1, 1)   # skewness (shape) parameter 
alpha <- rep(0, p)    # zero skewness




gen_azzalini <- function(n, p, alpha) {
  xi <- rep(0, p)  # location (mean)
  Omega <- diag(p) # covariance matrix (must be positive definite)
  
  rmsn(n , xi , Omega , alpha)
}



metric_skew(gen_azzalini(n,p,alpha))


# --------------------------------------------------------
# Sahu–Dey–Branco multivariate skew-normal generator, X ~ SN(xi, Sigma, Lambda)
# --------------------------------------------------------
rmsn_sdb <- function(n, xi, Lambda, Sigma) {
  xi <- as.vector(xi)
  M  <- length(xi)
  q  <- ncol(Lambda)
  
  if (!all(dim(Lambda) == c(M, q)))
    stop("Lambda must be M x q")
  
  if (!all(dim(Sigma) == c(M, M)))
    stop("Sigma must be M x M")
  
  # latent variables
  U <- matrix(rnorm(n * q), n, q)
  U_abs <- abs(U)
  
  # Gaussian noise
  Eps <- MASS::mvrnorm(n, mu = rep(0, M), Sigma = Sigma)
  
  # generate X
  X <- matrix(rep(xi, each = n), n, M) +
    U_abs %*% t(Lambda) +
    Eps
  
  X
}


n <- 200              # sample size
p <- 3                # dimension

### with the Lambda you can control the skewness
Lambda <- matrix(
  c( 2,  0,
     0,  1,
     -1,  1),
  nrow = p, byrow = TRUE
)

Lambda <- matrix(0, nrow = p, ncol = 2) # zero skewness

gen_sdb_sym <- function(n, p, Lambda) {
  xi <- rep(0, p)
  Sigma <- diag(p)
  rmsn_sdb(n, xi, Lambda, Sigma) 
}


metric_skew(gen_sdb_sym(n,p,Lambda))





# --------------------------------------------------------
# Permutation test for h-hat(p) (Metric Skewness) = symmetry
# --------------------------------------------------------
perm_test_metric <- function(X, B) {
  
  X <- as.matrix(X)
  n <- nrow(X)
  
  # observed statistic
  t_obs <- metric_skew(X)
  
  # permutation distribution (row-wise sign flipping)
  t_perm <- replicate(B, {
    signs <- sample(c(-1, 1), n, replace = TRUE)
    X_perm <- X * signs
    metric_skew(X_perm)
  })
  
  # two-sided p-value
  mean(abs(t_perm) >= abs(t_obs))
}




# --------------------------------------------------------
# permutation test for Mardia's skewness
# --------------------------------------------------------

perm_test_mardia <- function(X, B) {
  t_obs <- mardia_skewness(X)
  n <- nrow(X)
  
  t_perm <- replicate(B, {
    signs <- sample(c(-1, 1), n, replace = TRUE)
    mardia_skewness(X * signs)
  })
  
  mean(t_perm >= t_obs)
}

# --------------------------------------------------------
# Asymptotic test for Mardia multivariate skewness
# --------------------------------------------------------


asymp_pvalue_mardia <- function(X) {
  
  as.numeric(as.vector(mardia(X)$mv.test$`p-value`)[1])
  
}





############## Asymptotic test for Metric Skewness ###########################

distance_matrix_mv <- function(X) {
  # X: n x d matrix
  n <- nrow(X)
  D <- matrix(0, n, n)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      d2 <- sum((X[i, ] - X[j, ])^2)
      D[i, j] <- d2
      D[j, i] <- d2
    }
  }
  
  D
}


distance_matrix_to_flip <- function(X) {
  n <- nrow(X)
  Dflip <- matrix(0, n, n)
  
  for (i in 1:n) {
    Xi <- X[i, ]
    for (j in 1:n) {
      Dflip[i, j] <- sum((Xi + X[j, ])^2)
    }
  }
  
  Dflip
}


#### Kernels 

s1_kernel <- function(j, k, l, D) {
  (D[j, k] * D[j, l] +
     D[k, l] * D[k, j] +
     D[l, j] * D[l, k]) / 3
}

s2_kernel <- function(j, k, l, D, G) {
  r_jkl <- (D[j, k] - G[j, k]) * (D[j, l] - G[j, l])
  r_klj <- (D[k, l] - G[k, l]) * (D[k, j] - G[k, j])
  r_ljk <- (D[l, j] - G[l, j]) * (D[l, k] - G[l, k])
  
  (r_jkl + r_klj + r_ljk) / 3
}

############ s_hat_j #######################

s_hat_j <- function(j, D, G, p) {
  n <- nrow(D)
  if (n < 3) return(NA_real_)
  
  idx <- setdiff(seq_len(n), j)
  denom <- choose(n - 1, 2)
  total <- 0
  
  for (k in idx) {
    
    for (l in setdiff(seq(k, n), j)) {   # <-- L >= K
      
      if (p == 1) {
        total <- total + s1_kernel(j, k, l, D)
      } else if (p == 2) {
        total <- total + s2_kernel(j, k, l, D, G)
      }
    }
  }
  
  total / denom
}

################################################

s_hat_vector <- function(D, G, p) {
  n <- nrow(D)
  
  vapply(
    seq_len(n),
    function(j) s_hat_j(j, D, G, p),
    numeric(1)
  )
}

##### U statistics 

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


#### example

set.seed(1)
n <- 100
d <- 3

X <- matrix(rnorm(n * d), n, d)


D <- distance_matrix_mv(X)
G <- distance_matrix_to_flip(X)


U2_cpp <- u2_statistic_rcpp(D, G)
U2_r   <- u2_statistic(D, G)

all.equal(U2_cpp, U2_r)



#### Asymptotic Test

imhof_cdf <- function(x, weights) {
  res <- CompQuadForm::imhof(
    q = x,
    lambda = weights,
    epsabs = 1e-10,
    epsrel = 1e-10
  )
  
  # Defensive check (works across versions)
  if (!is.null(res$status) && res$status != 0) {
    warning("Imhof algorithm returned non-zero status.")
  }
  
  if (!is.finite(res$Qq)) {
    warning("Imhof returned non-finite CDF value.")
    return(NA_real_)
  }
  
  res$Qq
}


Asymp_metric_test <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  
  if (n < 3)
    return(list(statistic = NA, p.value = NA))
  
  # 1) Distance matrices
  D <- distance_matrix_mv(X)
  G <- distance_matrix_to_flip(X)
  
  # 2) U2 statistic (C++ version)
  U2 <- u2_statistic_rcpp(D, G)
  
  # 3) Sample covariance
  Sigma_hat <- cov(X)
  
  # 4) Eigenvalues and trace(Sigma^2)
  eigvals <- eigen(Sigma_hat, symmetric = TRUE, only.values = TRUE)$values
  theta_sq <- eigvals^2
  trSigma2 <- sum(theta_sq)
  
  # 5) Test statistic
  
  Tn <- ((n * U2) / 16) + trSigma2
  
  
  # 6) CDF via Imhof
  cdf_val <- imhof_cdf(Tn, theta_sq)

  pval <-  cdf_val

  
  ## OR two-sided p-value
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

X <- gen_azzalini(n = 100, p = 3, alpha = c(1, 0, 0))
Asymp_metric_test(X)





################################## Power Evaluation ########################

power_fixed_n <- function(
    n,
    alpha_grid,
    nrep,
    B,
    alpha
) {
  
  
  p <- length(alpha_grid[[1]])
  xi <- rep(0, p)
  Omega <- diag(p)
  
  power <- matrix(NA, nrow = length(alpha_grid), ncol = 4)
  colnames(power) <- c("Metric_perm", "Metric_asymp","Mardia_perm", "Mardia_asymp")
  
  for (k in seq_along(alpha_grid)) {
    
    alpha_skew <- alpha_grid[[k]]
    
    ncores = max(1, detectCores() - 1)
    
    pvals_list <- mclapply(seq_len(nrep), function(r) {
      
      X <- gen_azzalini(n, p, alpha_skew)
      
      
      
      c(
        Metric_perm  = perm_test_metric(X, B),
        Metric_asymp = Asymp_metric_test(X)$p.value,
        Mardia_perm  = perm_test_mardia(X, B),
        Mardia_asymp = asymp_pvalue_mardia(X)
      )
      
    }, mc.cores = ncores)
    
    
    # convert list → matrix
    pvals <- do.call(cbind, pvals_list)
    
    power[k, ] <- rowMeans(pvals < alpha, na.rm = TRUE)
    
    cat("n =", n,
        "alpha =", paste(alpha_skew, collapse = ","),
        "power =", round(power[k, ], 3), "\n")
  }
  
  data.frame(
    alpha_norm = sapply(alpha_grid, function(a) sqrt(sum(a^2))),
    power
  )
}




alpha_grid <- list(
  c(0, 0, 0),
  c(1, 0, 0),
  c(2, 0, 0),
  c(3, 0, 0),
  c(4, 0, 0)
)

start <- Sys.time()

res_n20  <- power_fixed_n(n = 20, alpha_grid = alpha_grid, nrep = 1000 , B = 500, 0.05)

res_n200 <- power_fixed_n(n = 200, alpha_grid = alpha_grid, nrep = 1000 , B = 500, 0.05)


end <- Sys.time()
running_time <- end - start


df_n20 <- res_n20 |>
  pivot_longer(
    cols = -alpha_norm,
    names_to = "Test",
    values_to = "Power"
  ) |>
  mutate(n = 20)

df_n200 <- res_n200 |>
  pivot_longer(
    cols = -alpha_norm,
    names_to = "Test",
    values_to = "Power"
  ) |>
  mutate(n = 200)

df_power <- bind_rows(df_n20, df_n200)



ggplot(
  df_power,
  aes(x = alpha_norm, y = Power, color = Test)
) +
  geom_line(linewidth = 1) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    color = "black",
    linewidth = 0.7
  )+
  geom_point(size = 2) +
  facet_wrap(~ n, labeller = label_both) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = expression("norm of "* alpha *" (Skewness magnitude)"),
    y = "RP",
    title = "Power Comparison under Azzalini Skew-Normal Data"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )







###################################################################################
# --------------------------------------------------------
# Parallel level test for different data generator
# --------------------------------------------------------
level_test_parallel <- function(gen_fun,
                                      sample_sizes,
                                      nrep,
                                      B,
                                      alpha,
                                      p,
                                      alpha_skew) {
  
  # for reproducibility across forked processes
  RNGkind("L'Ecuyer-CMRG")
  
  
  ncores = detectCores() - 1
  
  results_hhat         <- numeric(length(sample_sizes))
  results_hhat_asym    <- numeric(length(sample_sizes))
  results_mardia       <- numeric(length(sample_sizes))
  results_mardia_asym  <- numeric(length(sample_sizes))
  
  for (i in seq_along(sample_sizes)) {
    
    n <- sample_sizes[i]
    
    # ---- PARALLEL MONTE CARLO ----
    pvals_list <- mclapply(seq_len(nrep), function(r) {
      
      X <- gen_fun(n, p, alpha_skew)
      
      
      
      c(
        metric_perm  = perm_test_metric(X, B),
        metric_asym  = Asymp_metric_test(X)$p.value,
        mardia_perm  = perm_test_mardia(X, B),
        mardia_asym  = asymp_pvalue_mardia(X)
      )
      
    }, mc.cores = ncores)
    
    # convert list → matrix
    pvals <- do.call(cbind, pvals_list)
    
  
    
    results_hhat[i]      <- mean(pvals[1, ] < alpha, na.rm = TRUE)
    results_hhat_asym[i] <- mean(pvals[2, ] < alpha, na.rm = TRUE)
    results_mardia[i]    <- mean(pvals[3, ] < alpha, na.rm = TRUE)
    results_mardia_asym[i] <- mean(pvals[4, ] < alpha, na.rm = TRUE)
    
    cat(
      "Done n =", n,
      "| metric perm:", round(results_hhat[i], 3),
      "| metric asym:", round(results_hhat_asym[i], 3),
      "| mardia perm:", round(results_mardia[i], 3),
      "| mardia asym:", round(results_mardia_asym[i], 3), "\n"
    )
  }
  
  list(
    sample_sizes = sample_sizes,
    metric_perm = results_hhat,
    metric_asym = results_hhat_asym,
    mardia_perm = results_mardia,
    mardia_asym = results_mardia_asym
  )
}







set.seed(1)

## Parameters 
sample_sizes <- seq(20,300,20)
nrep = 10000
B = 500
alpha = 0.05
p = 3


# res <- level_test_parallel(
#   gen_fun = gen_mvn,
#   sample_sizes = sample_sizes,
#   nrep = 2,
#   B = 500,
#   alpha = 0.05,
#   p = 3
# )

start <- Sys.time()


res_sdb <- level_test_parallel(
  gen_fun = gen_sdb_sym,
  sample_sizes = sample_sizes,
  nrep,
  B,
  alpha ,
  p,
  alpha_skew = matrix(0, nrow = p, ncol = 2)
)


res_az <- level_test_parallel(
  gen_fun = gen_azzalini,
  sample_sizes = sample_sizes,
  nrep,
  B,
  alpha,
  p,
  alpha_skew = c(rep(0,p))
)


end <- Sys.time()

running_time <- end - start


#### Visualization

df_az <- tibble(
  n = res_az$sample_sizes,
  Metric_perm = res_az$metric_perm,
  Metric_asym = res_az$metric_asym,
  Mardia_perm = res_az$mardia_perm,
  Mardai_Asym = res_az$mardia_asym
) |>
  pivot_longer(-n, names_to = "Statistic", values_to = "Rejection") |>
  mutate(Dataset = "Azzalini")

df_sdb <- tibble(
  n = res_sdb$sample_sizes,
  Metric_perm = res_sdb$metric_perm,
  Metric_asym = res_sdb$metric_asym,
  Mardia_perm = res_sdb$mardia_perm,
  Mardai_Asym = res_sdb$mardia_asym
) |>
  pivot_longer(-n, names_to = "Statistic", values_to = "Rejection") |>
  mutate(Dataset = "SDB")

df_all <- bind_rows(df_az, df_sdb)



p_az <- ggplot(df_az,
               aes(x = n, y = Rejection,
                   color = Statistic, shape = Statistic)) +
  geom_line(linewidth = 0.9, linetype = 'solid') +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0.05,
             linetype = "dashed", color = "red") +
  scale_color_manual(values = c("black", "blue",'green','purple')) +
  scale_shape_manual(values = c(19, 17,18,15)) +
  labs(
    title = "Level Evaluation (Azzalini Data)",
    x = "Sample Size",
    y = expression("Proportion of Rejection (p < " * alpha * ")")
  ) + 
  ylim(0.02, 0.08) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )


p_sdb <- ggplot(df_sdb,
                aes(x = n, y = Rejection,
                    color = Statistic, shape = Statistic)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0.05,
             linetype = "dashed", color = "red") +
  scale_color_manual(values = c("black", "blue",'green','purple')) +
  scale_shape_manual(values = c(19, 17,18,15)) +
  labs(
    title = "Level Evaluation (SDB Data)",
    x = "Sample Size",
    y = expression("Proportion of Rejection (p < " * alpha * ")")
  ) +  
  ylim(0, 0.08) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )





(p_az | p_sdb) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.3, "cm"),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )








