
## Libraries 
library(dplyr)
Sys.setenv(RGL_USE_NULL=TRUE)
library(CovTools)
library(ICtest)
library(ggplot2)
library(fda)
library(igraph)

#####
remotes::install_github('vidazamani/Metric-Skewness-for-Object-Data/LogDis-R-Package@V3')
library(LogDis)


# Step 2: Compute average squared distance using distcov
average_squared_distances <- function(matrices) {
  n <- dim(matrices)[3]
  sapply(1:n, function(i) {
    mean(sapply(1:n,
                function(j) if (i != j) distcov(matrices[, , i],
                                                matrices[, , j],
                                                method ='LogEuclidean')^2 else 0))
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

# Step 4: Compute average squared distance for inverted matrices
average_squared_distances_inverted <- function(matrices) {
  iverted_matrices <- invert_matrices(matrices)
  n <- dim(matrices)[3]
  sapply(1:n, function(i) {
    mean(sapply(1:n, function(j) distcov(matrices[, , i], iverted_matrices[, , j],
                                         method ='LogEuclidean')^2))
  })
}



# --------------------------------------------------------
# Function to compute skewness
# --------------------------------------------------------

metric_skew_fun <- function(matrices) {
  
  # Step 2: Compute average squared distance for original matrices
  avg_dist_original <- average_squared_distances(matrices)
  
  
  # Step 3: Compute average squared distance for inverted matrices
  avg_dist_inverted <- average_squared_distances_inverted(matrices)
  
  # Step 4: Compute the final average squared distance
  final_avg_distance <- mean((avg_dist_original - avg_dist_inverted)^2)/mean(avg_dist_original^2)
  
  # Return results
  return(Metric_skewness = final_avg_distance)
}


############ Permutation test ##################################

Perm_test <- function(matrices, iter, seed, regularize){
  
  
  set.seed(seed)
  
  
  
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
  res <- CompQuadForm::imhof(q = x, 
                             lambda = weights,
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
  
  # 1) Log–Euclidean distances
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


################ Import data ######################
### There are 102 patients all together 

setwd('/Users/vizama/Documents/Papers/2nd paper/Dataset/PET data')


# vector of organ names
organNames <- c('spleen','kidney_right','kidney_left','liver','pancreas',
                'LUL','LLL','RUL','RML','RLL',
                'colon','urinary_bladder','heart','aorta','brain')

# get all csv files
all_files <- list.files()

# helper function to read and add organ Names
read_with_organs <- function(file) {
  df <- read.csv(file)
  df <- df[, -1]                          # remove first column
  df <- cbind(organ = organNames, df)     # add organ column
  return(df)
}



# split by category
rest_files   <- grep("rest", all_files, value = TRUE)
stress_files <- grep("stress", all_files, value = TRUE)



######################################################################
#### In this part, a data frame per each patient is created 


rest_data   <- lapply(rest_files, read_with_organs)
stress_data <- lapply(stress_files, read_with_organs)



get_pair_covariances <- function(data_list, organA, organB) {
  
  n_patients <- length(data_list)
  S_array <- array(0, dim = c(2, 2, n_patients))
  
  for (i in seq_len(n_patients)) {
    
    patient_df <- data_list[[i]]
    
    # Extract numeric time columns
    numeric_mat <- as.matrix(
      patient_df[, sapply(patient_df, is.numeric)]
    )
    
    # Set row names to organ names
    rownames(numeric_mat) <- patient_df$organ
    
    # Check organs exist
    if (!all(c(organA, organB) %in% rownames(numeric_mat))) {
      stop("Organ name not found in patient ", i)
    }
    
    Xi <- t(numeric_mat[c(organA, organB), ])
    
    S_array[, , i] <- cov(Xi)
  }
  
  S_array
}


get_pair_covariances(rest_data,'heart' ,'spleen')

get_organ_names <- function(data_list) {
  patient_df <- data_list[[1]]
  patient_df$organ
}



pairwise_metric_tests <- function(data_list,
                                  organs = NULL,
                                  n_organs = NULL,
                                  test_type = c("perm", "asymp"),
                                  B ,
                                  alpha = 0.05,
                                  seed = 1818) {
  
  test_type <- match.arg(test_type)
  
  all_organs <- get_organ_names(data_list)
  
  # ----- Organ selection logic -----
  
  if (!is.null(n_organs)) {
    if (!is.null(seed)) set.seed(seed)
    organs <- sample(all_organs, n_organs)
  }
  
  if (is.null(organs)) {
    organs <- all_organs
  }
  
  if (!all(organs %in% all_organs)) {
    stop("Some selected organs not found in dataset.")
  }
  
  # ----------------------------------
  
  pairs <- combn(organs, 2, simplify = FALSE)
  n_pairs <- length(pairs)
  
  results <- data.frame(
    organA = character(n_pairs),
    organB = character(n_pairs),
    p_value = numeric(n_pairs),
    stringsAsFactors = FALSE
  )
  
  for (k in seq_along(pairs)) {
    
    organA <- pairs[[k]][1]
    organB <- pairs[[k]][2]
    
    S_array <- get_pair_covariances(data_list, organA, organB)
    
    if (test_type == "perm") {
      test_res <- Perm_test(S_array, B, seed = 1, regularize = 0)
      pval <- test_res$p_value
    } else {
      test_res <- Asymp_metric_skewness_spd(S_array)
      pval <- test_res$p.value
    }
    
    results$organA[k] <- organA
    results$organB[k] <- organB
    results$p_value[k] <- pval
    
    cat("Done pair:", organA, "-", organB,
        "p =", round(pval, 4), "\n")
  }
  
  results
}






build_graph_from_pvalues <- function(results, alpha = 0.05) {
  
  organs <- unique(c(results$organA, results$organB))
  p <- length(organs)
  
  adj <- matrix(0, p, p)
  rownames(adj) <- colnames(adj) <- organs
  
  for (i in seq_len(nrow(results))) {
    
    if (results$p_value[i] < alpha) {
      
      A <- results$organA[i]
      B <- results$organB[i]
      
      adj[A, B] <- 1
      adj[B, A] <- 1
    }
  }
  
  adj
}


graph_from_adj <- function(adj_matrix) {
  graph_from_adjacency_matrix(adj_matrix,
                              mode = "undirected",
                              diag = FALSE)
}

plot_organ_graph <- function(adj_matrix) {
  
  g <- graph_from_adj(adj_matrix)
  
  plot(g,
       layout = layout_with_fr(g),   # Fruchterman-Reingold layout
       vertex.size = 30,
       vertex.label.cex = 0.8,
       edge.width = 2)
}



### Example


res_all <- pairwise_metric_tests(rest_data, test_type = "asymp")
res_all <- pairwise_metric_tests(rest_data, test_type = "perm", 5)

res_6 <- pairwise_metric_tests(rest_data,
                               n_organs = 6,
                               seed = 123,
                               test_type = "perm", B = 5)

res_subset <- pairwise_metric_tests(
  rest_data,
  organs = c("colon", "RML", "aorta"),
  test_type = "perm", B = 10
)

adj_matrix <- build_graph_from_pvalues(res_subset, alpha = 0.05)
adj_matrix <- build_graph_from_pvalues(res_all, alpha = 0.05)

plot_organ_graph(adj_matrix)




############################################################
pairwise_metric_tests_full <- function(data_list,
                                       organs = NULL,
                                       n_organs = NULL,
                                       B ,
                                       alpha = 0.05,
                                       seed = NULL) {
  
  # ----- Organ selection -----
  all_organs <- data_list[[1]]$organ
  
  if (!is.null(n_organs)) {
    if (!is.null(seed)) set.seed(seed)
    organs <- sample(all_organs, n_organs)
  }
  
  if (is.null(organs)) {
    organs <- all_organs
  }
  
  if (!all(organs %in% all_organs)) {
    stop("Some selected organs not found.")
  }
  
  pairs <- combn(organs, 2, simplify = FALSE)
  n_pairs <- length(pairs)
  
  results <- data.frame(
    organA = character(n_pairs),
    organB = character(n_pairs),
    perm_pvalue = numeric(n_pairs),
    asymp_pvalue = numeric(n_pairs),
    metric_skew_perm = numeric(n_pairs),
    U2 = numeric(n_pairs),
    U1 = numeric(n_pairs),
    effect_size_U2_over_U1 = numeric(n_pairs),
    stringsAsFactors = FALSE
  )
  
  for (k in seq_along(pairs)) {
    
    organA <- pairs[[k]][1]
    organB <- pairs[[k]][2]
    
    # ----- Get covariance matrices -----
    S_array <- get_pair_covariances(data_list, organA, organB)
    
    # ----- Compute distance matrices -----
    D <- distance_logeuclid_cpp(S_array)
    G <- distance_to_inverse_logeuclid_cpp(S_array)
    
    # ----- U-statistics -----
    U1 <- u1_statistic(D)
    U2 <- u2_statistic_rcpp(D, G)
    
    effect_size <- U2 / U1
    
    # ----- Metric skewness used in permutation -----
    metric_skew_value <- metric_skew_fun(S_array)
    
    # ----- Permutation test -----
    perm_res <- Perm_test(S_array, B, seed = 1, regularize = 0)
    perm_p <- perm_res$p_value
    
    # ----- Asymptotic test -----
    asymp_res <- Asymp_metric_skewness_spd(S_array)
    asymp_p <- asymp_res$p.value
    
    # ----- Store results -----
    results[k, ] <- c(
      organA,
      organB,
      perm_p,
      asymp_p,
      metric_skew_value,
      U2,
      U1,
      effect_size
    )
    
    cat("Done:", organA, "-", organB,
        "| perm =", round(perm_p, 4),
        "| asymp =", round(asymp_p, 4),
        "| Metric Skewness", round(metric_skew_value,4),
        "| effect =", round(effect_size, 4), "\n")
  }
  
  results
}







build_graph_from_results <- function(results,
                                     p_col = "perm_pvalue",
                                     alpha = 0.05) {
  
  organs <- unique(c(results$organA, results$organB))
  p <- length(organs)
  
  adj <- matrix(0, p, p)
  rownames(adj) <- colnames(adj) <- organs
  
  for (i in seq_len(nrow(results))) {
    
    if (!is.na(results[[p_col]][i]) &&
        results[[p_col]][i] < alpha) {
      
      A <- results$organA[i]
      B <- results$organB[i]
      
      adj[A, B] <- 1
      adj[B, A] <- 1
    }
  }
  
  graph_from_adjacency_matrix(adj,
                              mode = "undirected",
                              diag = FALSE)
}




build_weighted_graph <- function(results,
                                 weight_col = "effect_size_U2_over_U1",
                                 p_filter = TRUE,
                                 p_col = "perm_pvalue",
                                 alpha = 0.05) {
  
  edges <- results
  
  if (p_filter) {
    edges <- edges[edges[[p_col]] < alpha, ]
  }
  
  edges$weight_numeric <- as.numeric(edges[[weight_col]])
  
  g <- graph_from_data_frame(
    d = data.frame(
      from = edges$organA,
      to   = edges$organB,
      weight = edges$weight_numeric
    ),
    directed = FALSE
  )
  
  g
}


plot_organ_graph <- function(g,
                             vertex_size = 30,
                             show_weights = FALSE) {
  
  if (show_weights) {
    E(g)$label <- round(E(g)$weight, 3)
  }
  
  plot(g,
       layout = layout_with_fr(g),
       vertex.size = vertex_size,
       vertex.label.cex = 0.9,
       edge.width = 2 + 5 * E(g)$weight,
       main = "Organ Connectivity Graph")
}





#### Example

res_subset <- pairwise_metric_tests_full(rest_data, organs = c("heart", "spleen", "LLL"), B =5)


g_perm <- build_graph_from_results(res_subset,
                                   p_col = "perm_pvalue",
                                   alpha = 0.05)

plot_organ_graph(g_perm)


g_weighted <- build_weighted_graph(res_subset)
plot_organ_graph(g_weighted, show_weights = TRUE)




############ Extra ############################
find_maximal_cliques <- function(adj_matrix) {
  
  g <- graph_from_adjacency_matrix(adj_matrix,
                                   mode = "undirected",
                                   diag = FALSE)
  
  cliques(g)
}

find_maximal_cliques(adj_matrix)

largest_clique <- function(adj_matrix) {
  
  g <- graph_from_adjacency_matrix(adj_matrix,
                                   mode = "undirected",
                                   diag = FALSE)
  
  max_cliques(g)
}

largest_clique(adj_matrix)





graph_shape_summary <- function(adj_matrix) {
  
  g <- graph_from_adj(adj_matrix)
  
  list(
    n_nodes = vcount(g),
    n_edges = ecount(g),
    density = edge_density(g),
    n_components = components(g)$no,
    diameter = diameter(g, directed = FALSE),
    avg_clustering = transitivity(g, type = "average"),
    max_clique_size = max(sapply(cliques(g), length))
  )
}

graph_shape_summary(adj_matrix)


centrality_summary <- function(adj_matrix) {
  
  g <- graph_from_adj(adj_matrix)
  
  data.frame(
    organ = V(g)$name,
    degree = degree(g),
    betweenness = betweenness(g),
    closeness = closeness(g)
  )
}


centrality_summary(adj_matrix)

detect_communities <- function(adj_matrix) {
  
  g <- graph_from_adj(adj_matrix)
  cluster_louvain(g)
}


detect_communities(adj_matrix)
