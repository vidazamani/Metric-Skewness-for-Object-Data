## Libraries 
library(dplyr)
Sys.setenv(RGL_USE_NULL=TRUE)
library(CovTools)
library(ICtest)
library(ggplot2)
library(fda)
library(igraph)
library(parallel)
library(shapes)
library(readxl)
library(dplyr)

#####
remotes::install_github('vidazamani/Metric-Skewness-for-Object-Data/LogDis-R-Package@V3')
library(LogDis)




# Invert all matrices function
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

flip_cov <- function(S) {
  S_flipped <- S
  S_flipped[1,2] <- -S_flipped[1,2]
  S_flipped[2,1] <- -S_flipped[2,1]
  S_flipped
}


Perm_test <- function(matrices, iter, seed){
  
  set.seed(seed)
  
  T0 <- metric_skew_fun(matrices)
  
  n <- dim(matrices)[3]
  Tperm <- numeric(iter)
  
  for(b in 1:iter){
    
    flips <- sample(c(TRUE, FALSE), n, replace = TRUE)
    
    perm_sample <- matrices
    
    for(i in 1:n){
      if(flips[i]){
        perm_sample[,,i] <- flip_cov(matrices[,,i])
      }
    }
    
    Tperm[b] <- metric_skew_fun(perm_sample)
  }
  
  pval <- mean(abs(Tperm) >= abs(T0))
  
  list(statistic = T0, p_value = pval, Tperm = Tperm)
}






#################################################################

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
    numeric_mat <- as.matrix(patient_df[,-1])
    
    
    # Set row names to organ names
    rownames(numeric_mat) <- patient_df$organ
    
    # Check organs exist
    if (!all(c(organA, organB) %in% rownames(numeric_mat))) {
      stop("Organ name not found in patient ", i)
    }
    
    Xi <- t(numeric_mat[c(organA, organB), ])
    
    S <- cov(Xi)
    S <- S + 1e-6 * diag(2)
    #S_array[, , i] <- cov(Xi)
    S_array[, , i] <- S
  }
  
  S_array
}


get_pair_covariances(rest_data,'heart' ,'spleen')
get_pair_covariances(stress_data,'heart' ,'spleen')



get_organ_names <- function(data_list) {
  patient_df <- data_list[[1]]
  patient_df$organ
}




  


pairwise_metric_tests_parallel <- function(
    data_list,
    organs = NULL,
    n_organs = NULL,
    B,
    alpha = 0.05,
    seed = 1818
) {
  

  ncores <- detectCores() - 1
  
  all_organs <- get_organ_names(data_list)
  
  # ----- Organ selection logic -----
  
  if (!is.null(organs)) {
    
    # Scenario 1: user provided organ vector
    if (!all(organs %in% all_organs)) {
      stop("Some selected organs not found in dataset.")
    }
    
  } else if (!is.null(n_organs)) {
    
    # Scenario 2: randomly sample organs
    set.seed(seed)
    organs <- sample(all_organs, n_organs)
    
  } else {
    
    # Scenario 3: use all organs
    organs <- all_organs
    
  }
  
  # ----------------------------------
  
  pairs <- combn(organs, 2, simplify = FALSE)
  
  results_list <- mclapply(seq_along(pairs), function(k) {
    
    organA <- pairs[[k]][1]
    organB <- pairs[[k]][2]
    
    S_array <- get_pair_covariances(data_list, organA, organB)
    
    test_res <- Perm_test(S_array, B, seed)
    pval <- test_res$p_value
    
    cat("Done pair:", organA, "-", organB,
        "p =", round(pval, 4), "\n")
    
    data.frame(
      organA = organA,
      organB = organB,
      p_value = pval,
      stringsAsFactors = FALSE
    )
    
  }, mc.cores = ncores)
  
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL
  
  results
}


############# Graph ########


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
       layout = layout_with_fr(g), # Fruchterman-Reingold layout # or layout_in_circle(g)   
       vertex.size = 30,
       vertex.label.cex = 0.8,
       edge.width = 2)
}



######### Examples 

# Scenario 1 — Specific organs
res_subset <- pairwise_metric_tests_parallel(
  rest_data,
  organs = c("heart","colon","LUL"),
  B = 5000
)

# Scenario 2 — Random subset

res_subset <- pairwise_metric_tests_parallel(
  rest_data,
  n_organs = 6,
  seed = 123,
  B = 500
)

# Scenario 3 — All organs

res_all <- pairwise_metric_tests_parallel(
  rest_data,
  B = 1000
)



adj_matrix <- build_graph_from_pvalues(res_subset, alpha = 0.05)
adj_matrix <- build_graph_from_pvalues(res_all, alpha = 0.05)

plot_organ_graph(adj_matrix)

#### Full verjon of permutation test where we have the metric skewness value

pairwise_metric_tests_parallel <- function(
    data_list,
    organs = NULL,
    n_organs = NULL,
    B,
    alpha = 0.05,
    seed = 1818
) {
  

  ncores <- detectCores() - 1
  
  all_organs <- get_organ_names(data_list)
  
  # ----- Organ selection logic -----
  
  if (!is.null(organs)) {
    
    # Scenario 1: user provided organ vector
    if (!all(organs %in% all_organs)) {
      stop("Some selected organs not found in dataset.")
    }
    
  } else if (!is.null(n_organs)) {
    
    # Scenario 2: randomly sample organs
    set.seed(seed)
    organs <- sample(all_organs, n_organs)
    
  } else {
    
    # Scenario 3: use all organs
    organs <- all_organs
    
  }
  
  # ----------------------------------
  
  pairs <- combn(organs, 2, simplify = FALSE)
  
  results_list <- mclapply(seq_along(pairs), function(k) {
    
    organA <- pairs[[k]][1]
    organB <- pairs[[k]][2]
    
    S_array <- get_pair_covariances(data_list, organA, organB)
    
    # weight = metric skewness
    skew_value <- metric_skew_fun(S_array)
    
    # permutation test
    test_res <- Perm_test(S_array, B, seed + k)
    pval <- test_res$p_value
    
    cat("Done pair:", organA, "-", organB,
        "| skew =", round(skew_value,4),
        "| p =", round(pval,4), "\n")
    
    data.frame(
      organA = organA,
      organB = organB,
      metric_skewness = skew_value,
      p_value = pval,
      stringsAsFactors = FALSE
    )
    
  }, mc.cores = ncores)
  
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL
  
  results
}

##### Weighted Graph


build_weighted_graph <- function(results) {
  

  organs <- unique(c(results$organA, results$organB))
  p <- length(organs)
  
  adj <- matrix(0, p, p)
  rownames(adj) <- colnames(adj) <- organs
  
  for (i in seq_len(nrow(results))) {
    
    A <- results$organA[i]
    B <- results$organB[i]
    
    
      w <- as.numeric(results$metric_skewness[i])
    
    
    adj[A, B] <- w
    adj[B, A] <- w
  }
  
  adj
}



#### Example


res_all_rest <- pairwise_metric_tests_parallel(
  rest_data,
  B = 1000
)

res_all_stress <- pairwise_metric_tests_parallel(
  stress_data,
  B = 1000
)




adj_weighted_rest <- build_weighted_graph(res_all_rest)
adj_weighted_stress <- build_weighted_graph(res_all_stress)


thr <- quantile(adj_weighted_rest[adj_weighted_rest > 0], 0.13)

adj_filtered <- adj_weighted_rest
adj_filtered[adj_filtered < thr] <- 0

g_rest <- graph_from_adjacency_matrix(adj_filtered,
                                 mode = "undirected",
                                 weighted = TRUE,
                                 diag = FALSE)

w <- abs(E(g_rest)$weight)

w_scaled <- log1p(w)
w_scaled <- 1 + 5*(w_scaled - min(w_scaled)) /
  (max(w_scaled) - min(w_scaled))

plot(
  g_rest,
  layout = layout_with_fr(g_rest),
  vertex.size = 35,
  edge.width = w_scaled,
  vertex.label.cex = 0.8
)




thr <- quantile(adj_weighted_stress[adj_weighted_stress > 0], 0.13)

adj_filtered <- adj_weighted_stress
adj_filtered[adj_filtered < thr] <- 0

g_stress <- graph_from_adjacency_matrix(adj_filtered,
                                      mode = "undirected",
                                      weighted = TRUE,
                                      diag = FALSE)

w <- abs(E(g_stress)$weight)

w_scaled <- log1p(w)
w_scaled <- 1 + 5*(w_scaled - min(w_scaled)) /
  (max(w_scaled) - min(w_scaled))

plot(
  g_stress,
  layout = layout_with_fr(g_stress),
  vertex.size = 35,
  edge.width = w_scaled,
  vertex.label.cex = 0.8
)




# mod_rest   <- modularity(cluster_louvain(g_rest))
# mod_stress <- modularity(cluster_louvain(g_stress))
# 
# mod_rest > mod_stress




clusters <- cluster_louvain(g_rest)
membership(clusters)
modularity(clusters)


############ Extra ############################
# find_maximal_cliques <- function(adj_matrix) {
#   
#   g <- graph_from_adjacency_matrix(adj_matrix,
#                                    mode = "undirected",
#                                    diag = FALSE)
#   
#   cliques(g)
# }
# 
# find_maximal_cliques(adj_matrix)
# 
# largest_clique <- function(adj_matrix) {
#   
#   g <- graph_from_adjacency_matrix(adj_matrix,
#                                    mode = "undirected",
#                                    diag = FALSE)
#   
#   max_cliques(g)
# }
# 
# largest_clique(adj_matrix)
# 
# 
# 
# 
# 
# graph_shape_summary <- function(adj_matrix) {
#   
#   g <- graph_from_adj(adj_matrix)
#   
#   list(
#     n_nodes = vcount(g),
#     n_edges = ecount(g),
#     density = edge_density(g),
#     n_components = components(g)$no,
#     diameter = diameter(g, directed = FALSE),
#     avg_clustering = transitivity(g, type = "average"),
#     max_clique_size = max(sapply(cliques(g), length))
#   )
# }
# 
# graph_shape_summary(adj_matrix)
# 
# 
# centrality_summary <- function(adj_matrix) {
#   
#   g <- graph_from_adj(adj_matrix)
#   
#   data.frame(
#     organ = V(g)$name,
#     degree = degree(g),
#     betweenness = betweenness(g),
#     closeness = closeness(g)
#   )
# }
# 
# 
# centrality_summary(adj_matrix)
# 
# detect_communities <- function(adj_matrix) {
#   
#   g <- graph_from_adj(adj_matrix)
#   cluster_louvain(g)
# }
# 
# 
# detect_communities(adj_matrix)
# 
# 
# 
# 
# 
# 



###########################


clinical_data <- read_excel("/Users/vizama/Documents/Papers/2nd paper/Dataset/excel_aboutKoveriData.xlsx")


clinical_data <- clinical_data %>%
  mutate(across(everything(),
                ~na_if(., "NA"))) %>%
  mutate(across(everything(),
                ~na_if(., "")))
# Remove rows with ANY NA
clinical_data_clean <- clinical_data %>%
  filter(complete.cases(.))

# Extract patient IDs from filenames
rest_ids <- sub("_rest_15tacs.csv", "", rest_files)
stress_ids <- sub("_stress_15tacs.csv", "", stress_files)

rest_ids

rest_data   <- lapply(rest_files, read_with_organs)
stress_data <- lapply(stress_files, read_with_organs)


names(rest_data)   <- rest_ids
names(stress_data) <- stress_ids


valid_ids <- clinical_data_clean$study_code

rest_data_clean <- rest_data[names(rest_data) %in% valid_ids]
stress_data_clean <- stress_data[names(stress_data) %in% valid_ids]

# clinical_numeric <- clinical_data_clean %>%
#   mutate(
#     weight = as.numeric(gsub(",", ".", weight_in_kg)),
#     height = as.numeric(gsub(",", ".", height_in_cm))
#   )
# 
# clinical_numeric <- clinical_numeric %>%
#   mutate(
#     BMI = weight / (height / 100)^2
#   )
# 
# 
# summary(clinical_numeric$weight)
# summary(clinical_numeric$height)
# summary(clinical_numeric$BMI)
# 
# sd(clinical_numeric$weight, na.rm = TRUE)
# sd(clinical_numeric$BMI, na.rm = TRUE)
# 
# range(clinical_numeric$weight, na.rm = TRUE)
# range(clinical_numeric$BMI, na.rm = TRUE)
# 
# 
# compute_patient_integration <- function(patient_df) {
#   
#   # extract numeric time columns
#   numeric_mat <- as.matrix(
#     patient_df[, sapply(patient_df, is.numeric)]
#   )
#   
#   rownames(numeric_mat) <- patient_df$organ
#   
#   # transpose: time × organs
#   X <- t(numeric_mat)
#   
#   # covariance across organs
#   S <- cov(X)
#   
#   # remove diagonal
#   off_diag <- S[lower.tri(S)]
#   
#   mean(abs(off_diag))
# }
# 
# names(rest_data_clean)
# names(stress_data_clean)
# 
# # rest_data and stress_data are lists of 102 patient data frames
# rest_metric   <- sapply(rest_data_clean, compute_patient_integration)
# stress_metric <- sapply(stress_data_clean, compute_patient_integration)
# 
# delta_metric <- stress_metric - rest_metric
# 
# patient_metrics <- data.frame(
#   study_code = names(rest_metric),
#   rest_integration = rest_metric,
#   stress_integration = stress_metric,
#   delta_integration = delta_metric
# )
# 
# analysis_df <- clinical_numeric %>%
#   inner_join(patient_metrics, by = "study_code")
# 
# model <- lm(delta_integration ~ BMI, data = analysis_df)
# summary(model)
# 
# 
# t.test(analysis_df$stress_integration,
#        analysis_df$rest_integration,
#        paired = TRUE)
# 




######################### Examples ################
# res_rest <- pairwise_metric_tests_full(
#   rest_data_clean,
#   B = 2
# )
# 
# res_stress <- pairwise_metric_tests_full(
#   stress_data_clean,
#   B = 2
# )
# 
# 
# 
# adj_rest   <- build_weighted_graph(res_rest, "effect")
# adj_stress <- build_weighted_graph(res_stress, type="effect")
# 
# g_rest   <- graph_from_adjacency_matrix(adj_rest,
#                                         mode="undirected",
#                                         weighted=TRUE,
#                                         diag=FALSE)
# 
# g_stress <- graph_from_adjacency_matrix(adj_stress,
#                                         mode="undirected",
#                                         weighted=TRUE,
#                                         diag=FALSE)


ischemia_ids <- clinical_data_clean %>%
  filter(`ischemia/YES/NO` == "YES") %>%
  pull(study_code)

# rest_ischemia   <- rest_data_clean[names(rest_data_clean) %in% ischemia_ids]
# rest_noischemia <- rest_data_clean[!names(rest_data_clean) %in% ischemia_ids]


stress_ischemia   <- stress_data_clean[names(stress_data_clean) %in% ischemia_ids]
stress_noischemia <- stress_data_clean[!names(stress_data_clean) %in% ischemia_ids]


# res_rest_isch <- pairwise_metric_tests_full(rest_ischemia, B=2)
# res_rest_no   <- pairwise_metric_tests_full(rest_noischemia, B=2)

st_rest_isch <- pairwise_metric_tests_full(stress_ischemia, B=2)
st_rest_no   <- pairwise_metric_tests_full(stress_noischemia, B=2)

# adj_rest_isch  <- build_weighted_graph(res_rest_isch, "effect")
# adj_rest_no <- build_weighted_graph(res_rest_no, "effect")

adj_st_isch  <- build_weighted_graph(st_rest_isch, "effect")
adj_st_no <- build_weighted_graph(st_rest_no, "effect")

# g_rest_isch  <- graph_from_adjacency_matrix(adj_rest_isch,
#                                         mode="undirected",
#                                         weighted=TRUE,
#                                         diag=FALSE)

# g_rest_no <- graph_from_adjacency_matrix(adj_rest_no,
#                                         mode="undirected",
#                                         weighted=TRUE,
#                                         diag=FALSE)



g_st_isch  <- graph_from_adjacency_matrix(adj_st_isch,
                                            mode="undirected",
                                            weighted=TRUE,
                                            diag=FALSE)

g_st_no <- graph_from_adjacency_matrix(adj_st_no,
                                         mode="undirected",
                                         weighted=TRUE,
                                         diag=FALSE)


w_stisch <- E(g_st_isch)$weight

# Normalize to [1, 10] range
w_scaled_stisch <- 1 + 9 * (w_stisch - min(w_stisch)) / (max(w_stisch) - min(w_stisch))


w_stno <- E(g_st_no)$weight

# Normalize to [1, 10] range
w_scaled_stno <- 1 + 9 * (w_stno - min(w_stno)) / (max(w_stno) - min(w_stno))




# layout_fixed <- layout_with_fr(g_rest_isch)
# plot(g_rest_isch, layout=layout_fixed)
# 
# layout_fixed <- layout_with_fr(g_rest_no)
# plot(g_rest_no, layout=layout_fixed)



layout_fixed <- layout_with_fr(g_st_isch)
plot(g_st_isch, layout=layout_fixed, edge.width =w_scaled_stisch)

layout_fixed <- layout_with_fr(g_st_no)
plot(g_st_no, layout=layout_fixed, edge.width = w_scaled_stno)
