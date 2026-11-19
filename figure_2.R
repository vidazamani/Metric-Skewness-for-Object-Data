library(sn)

set.seed(2026)

## ---------------------------------------------------------------
## Settings
## ---------------------------------------------------------------
n     <- 10000                 # "medium" sample size
a_seq <- c(1, 4, 20)         # slant (skewness) left -> right
s_seq <- c(1, 5, 10)         # elongation in the y-direction

## ---------------------------------------------------------------
## Zero-mean SN_2 skewed and elongated along the y-axis
##   Omega = diag(1, s^2), alpha = (0, a)
##   E(X)  = xi + omega * delta * sqrt(2/pi) = 0
## ---------------------------------------------------------------
make_sn <- function(a, s) {
  delta <- a / sqrt(1 + a^2)
  list(xi    = c(0, -s * delta * sqrt(2 / pi)),
       Omega = diag(c(1, s^2)),
       alpha = c(0, a))
}

## ---------------------------------------------------------------
## Estimator of h(P)_a with g(x) = -x and Euclidean distance
##   m_a(X_i; P)   ~ mean_{j != i} ||X_i - X_j||^a
##   m_a(X_i; P_g) ~ mean_{j != i} ||X_i + X_j||^a
## ---------------------------------------------------------------
h_est <- function(X, a) {
  sq  <- rowSums(X^2)
  G   <- X %*% t(X)
  D   <- sqrt(pmax(outer(sq, sq, "+") - 2 * G, 0))   # d(X_i,  X_j)
  Dg  <- sqrt(pmax(outer(sq, sq, "+") + 2 * G, 0))   # d(X_i, -X_j)
  Da  <- D^a;  Dga <- Dg^a
  diag(Da) <- NA; diag(Dga) <- NA                    # leave-one-out
  m   <- rowMeans(Da,  na.rm = TRUE)
  mg  <- rowMeans(Dga, na.rm = TRUE)
  mean((m - mg)^2) / mean(m^2)
}

## ---------------------------------------------------------------
## Plot (common axes for all panels)
## ---------------------------------------------------------------
xlim <- c(-4, 4)
ylim <- c(-10, 25)
gx   <- seq(xlim[1], xlim[2], length.out = 200)
gy   <- seq(ylim[1], ylim[2], length.out = 500)
grid <- as.matrix(expand.grid(gx, gy))


pdf("YOUR FOLDER HERE/fig_example_skewnesses.pdf", width = 10, height = 4)

op <- par(mfrow = c(1, 3), mar = c(4, 4, 4.5, 1))

for (k in seq_along(a_seq)) {
  p <- make_sn(a_seq[k], s_seq[k])
  
  z <- matrix(dmsn(grid, xi = p$xi, Omega = p$Omega, alpha = p$alpha),
              nrow = length(gx))
  
  X  <- rmsn(n, xi = p$xi, Omega = p$Omega, alpha = p$alpha)
  h1 <- h_est(X, 1)
  h2 <- h_est(X, 2)
  h3 <- h_est(X, 3)
  
  contour(gx, gy, z, nlevels = 10, drawlabels = FALSE, col = "steelblue",
          xlim = xlim, ylim = ylim,
          xlab = expression(x[1]), ylab = expression(x[2]))
  # points(X, pch = 16, cex = 0.35, col = adjustcolor("grey30", 0.35))
  abline(h = 0, v = 0, lty = 3, col = "grey60")
  points(0, 0, pch = 3, col = "red", lwd = 2)          # true mean (0,0)
  
  title(main = bquote(hat(h)[1](P) == .(sprintf("%.4f", h1)) * "," ~~
                      hat(h)[2](P) == .(sprintf("%.4f", h2)) * "," ~~
                        hat(h)[3](P) == .(sprintf("%.4f", h3))),
        line = 2.2)
  # mtext(bquote(a == .(a_seq[k]) * "," ~~ s == .(s_seq[k])),
  #       side = 3, line = 0.6, cex = 0.8)
}

par(op)

dev.off()

