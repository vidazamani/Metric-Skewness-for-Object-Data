# --------------------------------------------------------
# Parameters
# --------------------------------------------------------
alpha <- 0.05
nrep  <- 100                    # Repetitions for each sample size
sample_sizes <- seq(20, 300, 30)  # Example sample sizes

# --------------------------------------------------------
# Function to compute skewness
# --------------------------------------------------------
skew_fun <- function(x) {
  m <- mean(x)
  n <- length(x)
  sum((x - m)^3) / n / (sd(x)^3)
}

# --------------------------------------------------------
# Permutation test for skewness = 0
# --------------------------------------------------------
perm_test_skew <- function(x, B = 1000) {
  
  # observed skewness
  t_obs <- skew_fun(x)
  
  # permutation distribution: sign flipping preserves symmetry
  t_perm <- replicate(B, {
    x_perm <- x * sample(c(-1, 1), length(x), replace = TRUE)
    skew_fun(x_perm)
  })
  
  # two-sided p-value
  pval <- mean(abs(t_perm) >= abs(t_obs))
  return(pval)
}

# --------------------------------------------------------
# Simulation to estimate level
# --------------------------------------------------------
results <- numeric(length(sample_sizes))

for(i in seq_along(sample_sizes)) {
  
  n <- sample_sizes[i]
  
  pvals <- replicate(nrep, {
    x <- rnorm(n, 0, 1)       # generate symmetric data
    perm_test_skew(x)         # compute p-value
  })
  
  results[i] <- mean(pvals < alpha)
  cat("Done n =", n, " → rejection proportion:", results[i], "\n")
}

# --------------------------------------------------------
# Plot: sample size vs proportion of rejections
# --------------------------------------------------------
plot(sample_sizes, results, type = "b", pch = 19,
     xlab = "Sample Size",
     ylab = "Proportion Rejected (p < alpha)",
     main = "Level Evaluation of Permutation Skewness Test")
abline(h = alpha, col = "red", lty = 2)   # ideal level line
