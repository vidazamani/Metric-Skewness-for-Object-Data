# Testing the asymptotic distribution
# for Euclidean data in the case p = 1

# We use data from N(0, 1), which means that the covariance matrix \Sigma
# in Theorem 3 is equal to the variance of this distribution = 1.
# This means that tr(\Sigma^2) = 1 and \theta_1 = 1 as well

# Thus Theorem 3 says that n U2/16 + 1 converges to the distribution
# \chi^2_1 and we next verify this with a simulation

library(Rcpp)
sourceCpp("/Users/jomivi/Library/Mobile Documents/com~apple~CloudDocs/Work/Supervision/Vida/Paper_2/u_statistic_fast.cpp")

n <- 200
iter <- 1000

res <- replicate(iter, {
  X <- rnorm(n)
  n*u_statistic_rcpp(X)/16 + 1
})

hist(res, breaks = 30, freq = FALSE)
curve(dchisq(x, df = 1), add = TRUE)

