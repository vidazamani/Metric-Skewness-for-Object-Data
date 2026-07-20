library(LogDis)
library(mvdistmetricskew)
library(RcppHungarian)


generate_rademacher <- function(n, k, p = 0.5){
  matrix(sample(c(-1,1), n*k, replace = TRUE, prob = c(1-p, p)),
         nrow = n, ncol = k)
}


set.seed(1)

n <- 100
k <- 5

X <- generate_rademacher(n, k)

wasserstein_test <- function(X, R, B){
  pval_vec <- rep(0, R)
  
  
  n <- nrow(X)
  
  # Distance matrices
  D <- distance_matrix_mv_cpp(X)
  G <- distance_matrix_to_flip_cpp(X)
  
  
  
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


set.seed(1)

n <- 100
k <- 5

X <- generate_rademacher(n, k)

wasserstein_test(X, R = 50, B = 100)



X_asym <- generate_rademacher(n, k, p = 0.7)

wasserstein_test(X_asym, R = 50, B = 100)



estimate_power <- function(n, k, p, R, B, n_sim = 100, alpha = 0.05){
  
  rejections <- rep(0, n_sim)
  
  for(i in 1:n_sim){
    
    X <- generate_rademacher(n, k, p)
    
    pval <- wasserstein_test(X, R, B)
    
    rejections[i] <- (pval < alpha)
  }
  
  mean(rejections)
}


n <- 100
k <- 5

p_values <- seq(0.5, 0.7, by = 0.02)

power_vals <- sapply(p_values, function(p){
  estimate_power(n, k, p, R = 30, B = 50, n_sim = 1000)
})


plot(p_values, power_vals,
     type = "b",
     pch = 19,
     xlab = "Rademacher parameter p",
     ylab = "Power",
     ylim = c(0,1),
     main = "Power of Wasserstein Symmetry Test")
abline(h = 0.05)

