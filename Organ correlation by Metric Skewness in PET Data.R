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
library(patchwork)
library(ggraph)
library(colorRamps)


#####
#remotes::install_github('vidazamani/Metric-Skewness-for-Object-Data/LogDis-R-Package@V3')
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
clinical_data <- read_excel("/Users/vizama/Documents/Papers/2nd paper/Dataset/excel_aboutKoveriData.xlsx")



clinical_data <- clinical_data %>%
  mutate(across(everything(),
                ~na_if(., "NA"))) %>%
  mutate(across(everything(),
                ~na_if(., "")))
# Remove rows with ANY NA
clinical_data_clean <- clinical_data %>%
  filter(complete.cases(.)) %>% filter(study_code != "koveri0081")


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


# Extract patient IDs from filenames
rest_ids <- sub("_rest_15tacs.csv", "", rest_files)
stress_ids <- sub("_stress_15tacs.csv", "", stress_files)




######################################################################
#### In this part, a data frame per each patient is created 

rest_data   <- lapply(rest_files, read_with_organs)
stress_data <- lapply(stress_files, read_with_organs)


names(rest_data)   <- rest_ids
names(stress_data) <- stress_ids


valid_ids <- clinical_data_clean$study_code

rest_data_clean <- rest_data[names(rest_data) %in% valid_ids]
stress_data_clean <- stress_data[names(stress_data) %in% valid_ids]


############## correlation between two different organs ######################

get_pair_corr <- function(data_list, organA, organB) {
  
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
      stop("Organ name not found in patient", i)
    }
    
    Xi <- t(numeric_mat[c(organA, organB), ])
    
    S <- cor(Xi, method = "kendall")
    S <- S + 1e-6 * diag(2)
    #S_array[, , i] <- cov(Xi)
    S_array[, , i] <- S
  }
  
  S_array
}


get_pair_corr(rest_data_clean,'heart' ,'spleen')
get_pair_corr(stress_data_clean,'heart' ,'spleen')



# hist(get_pair_corr(rest_data,'RLL' ,'RML')[1,2,], breaks = 10)

get_organ_names <- function(data_list) {
  patient_df <- data_list[[1]]
  patient_df$organ
}


# ############# Graph ########
# 
# 
# build_graph_from_pvalues <- function(results, alpha = 0.05) {
#   
#   organs <- unique(c(results$organA, results$organB))
#   p <- length(organs)
#   
#   adj <- matrix(0, p, p)
#   rownames(adj) <- colnames(adj) <- organs
#   
#   for (i in seq_len(nrow(results))) {
#     
#     if (results$p_value[i] < alpha) {
#       
#       A <- results$organA[i]
#       B <- results$organB[i]
#       
#       adj[A, B] <- 1
#       adj[B, A] <- 1
#     }
#   }
#   
#   adj
# }
# 
# 
# graph_from_adj <- function(adj_matrix) {
#   graph_from_adjacency_matrix(adj_matrix,
#                               mode = "undirected",
#                               diag = FALSE)
# }
# 
# plot_organ_graph <- function(adj_matrix) {
#   
#   g <- graph_from_adj(adj_matrix)
#   
#   plot(g,
#        layout = layout_with_fr(g), # Fruchterman-Reingold layout # or layout_in_circle(g)   
#        vertex.size = 30,
#        vertex.label.cex = 0.8,
#        edge.width = 2)
# }



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
    
    S_array <- get_pair_corr(data_list, organA, organB)
    
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
    
    
    w <- as.numeric(results$rank_skew[i])
    
    
    adj[A, B] <- w
    adj[B, A] <- w
  }
  
  adj
}


######### Examples 

# Scenario 1 — Specific organs
res_subset <- pairwise_metric_tests_parallel(
  rest_data_clean,
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

res_all_rest <- pairwise_metric_tests_parallel(
  rest_data_clean,
  B = 5000
)


res_all_stress <- pairwise_metric_tests_parallel(
  stress_data_clean,
  B = 5000
)


# adj_matrix <- build_graph_from_pvalues(res_subset, alpha = 0.0001)
# adj_matrix <- build_graph_from_pvalues(res_all, alpha = 0.0001)
#plot_organ_graph(adj_matrix)



res_all_rest <- res_all_rest %>% mutate(rank_skew = rank(metric_skewness))
res_all_stress <- res_all_stress %>% mutate(rank_skew = rank(metric_skewness))






adj_weighted_rest <- build_weighted_graph(res_all_rest)
adj_weighted_stress <- build_weighted_graph(res_all_stress)




# w <- abs(E(g_rest)$weight)
# 
# w_scaled <- log1p(w)
# w_scaled <- 1 + 5*(w_scaled - min(w_scaled)) /
#   (max(w_scaled) - min(w_scaled))
# 
# plot(
#   g_rest,
#   layout = layout_with_fr(g_rest),
#   vertex.size = 35,
#   edge.width = w_scaled,
#   vertex.label.cex = 0.8
# )





par(mfrow = c(1, 2))  


g_rest <- graph_from_adjacency_matrix(adj_weighted_rest,
                                      mode = "undirected",
                                      weighted = TRUE,
                                      diag = FALSE)
set.seed(42)
L <- layout_with_fr(g_rest)          # compute once

plot(
  g_rest,
  layout = L,
  vertex.size = 35,
  edge.width = 2,
  vertex.label.cex = 0.8,
  edge.color = rev(gray.colors(105))[res_all_rest$rank_skew]
)

mtext("Rest", side = 1, line = 0)


g_stress <- graph_from_adjacency_matrix(adj_weighted_stress,
                                        mode = "undirected",
                                        weighted = TRUE,
                                        diag = FALSE)
set.seed(42)
L <- layout_with_fr(g_stress)          # compute once



plot(
  g_stress,
  layout = L,
  vertex.size = 35,
  edge.width = 1.5,
  vertex.label.cex = 0.8,
  edge.color = rev(gray.colors(105))[res_all_stress$rank_skew]
)

mtext("Stress", side = 1, line = 0)




# mod_rest   <- modularity(cluster_louvain(g_rest))
# mod_stress <- modularity(cluster_louvain(g_stress))
# 
# mod_rest > mod_stress
# L <- layout_with_fr(g_stress) 
# 
# p_stress <- ggraph(g_stress, layout = "manual", x = L[,1], y = L[,2]) +
#   geom_edge_link(aes(width = weight, color = res_all_stress$rank_skew),
#                  alpha = 0.8) +
#   geom_node_point(size = 5) +
#   geom_node_text(aes(label = name), repel = TRUE, size = 5) +
#   scale_edge_color_gradientn(colors = rev(gray.colors(105))) +
#   theme_void()
# 
# 
# 
# 
# 
# L2 <- layout_with_fr(g_rest)
# 
# p_rest <- ggraph(g_rest, layout = "manual", x = L2[,1], y = L2[,2]) +
#   geom_edge_link(aes(width = weight, color = res_all_rest$rank_skew),
#                  alpha = 0.8) +
#   geom_node_point(size = 5) +
#   geom_node_text(aes(label = name), repel = TRUE, size = 5) +
#   scale_edge_color_gradientn(colors = rev(gray.colors(105))) +
#   theme_void()
# 
# 
# p_stress + p_rest

# clusters <- cluster_louvain(g_rest)
# membership(clusters)
# modularity(clusters)


ischemia_ids <- clinical_data_clean %>%
  filter(`ischemia/YES/NO` == "YES") %>%
  pull(study_code)

rest_ischemia   <- rest_data_clean[names(rest_data_clean) %in% ischemia_ids]
rest_noischemia <- rest_data_clean[!names(rest_data_clean) %in% ischemia_ids]





res_rest_isch <- pairwise_metric_tests_parallel(rest_ischemia, B=5000)
res_rest_isch <- res_rest_isch %>% mutate(rank_skew = rank(metric_skewness))


res_rest_no   <- pairwise_metric_tests_parallel(rest_noischemia, B=5000)
res_rest_no <- res_rest_no %>% mutate(rank_skew = rank(metric_skewness))





adj_rest_isch  <- build_weighted_graph(res_rest_isch)
adj_rest_no <- build_weighted_graph(res_rest_no)



g_rest_isch  <- graph_from_adjacency_matrix(adj_rest_isch,
                                            mode="undirected",
                                            weighted=TRUE,
                                            diag=FALSE)

g_rest_no <- graph_from_adjacency_matrix(adj_rest_no,
                                         mode="undirected",
                                         weighted=TRUE,
                                         diag=FALSE)






# w_stisch <- E(g_st_isch)$weight
# 
# # Normalize to [1, 10] range
# w_scaled_stisch <- 1 + 9 * (w_stisch - min(w_stisch)) / (max(w_stisch) - min(w_stisch))
# 
# 
# w_stno <- E(g_st_no)$weight
# 
# # Normalize to [1, 10] range
# w_scaled_stno <- 1 + 9 * (w_stno - min(w_stno)) / (max(w_stno) - min(w_stno))
# 



# layout_fixed <- layout_with_fr(g_rest_isch)
# plot(g_rest_isch, layout=layout_fixed)



# layout_fixed <- layout_with_fr(g_rest_no)
# plot(g_rest_no, layout=layout_fixed)


col <- scales::seq_gradient_pal("blue", "white", "Lab")(seq(0,1,length.out=105))
col <- gray.colors(105)

par(mfrow= c(1,2))

set.seed(1212)
layout_fixed <- layout_with_fr(g_rest_no)
# plot(g_st_isch, layout=layout_fixed, edge.width =w_scaled_stisch)
plot(g_rest_isch, layout=layout_fixed,
     vertex.size = 35,
     edge.width = 2,
     vertex.label.cex = 0.8,
     edge.color = col[res_rest_isch$rank_skew])

mtext("Rest Ischemia", side = 1, line = 0)



# layout_fixed <- layout_with_fr(g_rest_isch)
# plot(g_st_no, layout=layout_fixed, edge.width = w_scaled_stno)
plot(g_rest_no, layout=layout_fixed,
     vertex.size = 35,
     edge.width = 2,
     vertex.label.cex = 0.8,
     edge.color = col[res_rest_no$rank_skew])

mtext("Rest No Ischemia", side = 1, line = 0)


#### Stress


stress_ischemia   <- stress_data_clean[names(stress_data_clean) %in% ischemia_ids]
stress_noischemia <- stress_data_clean[!names(stress_data_clean) %in% ischemia_ids]

res_stress_isch <- pairwise_metric_tests_parallel(stress_ischemia, B=5000)
res_stress_isch <- res_stress_isch %>% mutate(rank_skew = rank(metric_skewness))

res_stress_no   <- pairwise_metric_tests_parallel(stress_noischemia, B=5000)
res_stress_no <- res_stress_no %>% mutate(rank_skew = rank(metric_skewness))



adj_st_isch  <- build_weighted_graph(res_stress_isch)

adj_st_no <- build_weighted_graph(res_stress_no)

g_st_isch  <- graph_from_adjacency_matrix(adj_st_isch,
                                          mode="undirected",
                                          weighted=TRUE,
                                          diag=FALSE)

g_st_no <- graph_from_adjacency_matrix(adj_st_no,
                                       mode="undirected",
                                       weighted=TRUE,
                                       diag=FALSE)


par(mfrow= c(2,2), mar = c(2, 2, 2, 2))

set.seed(1111)


###############################################


library(tibble)
library(dplyr)
library(igraph)

organ_coords <- tribble(
  ~name,              ~x,    ~y,
  "brain",             0.0,  10.0,
  "RUL",              -1.2,   7.1,
  "RML",              -1.2,   5.9,
  "RLL",              -1.2,   4.6,
  "LUL",               1.2,   7.0,
  "LLL",               1.2,   5.0,
  "heart",             0.0,   6.0,
  "aorta",             0.0,   4.8,
  "liver",            -1.3,   3.4,
  "spleen",            1.7,   3.3,
  "pancreas",          0.2,   2.8,
  "kidney_right",     -1.4,   2.1,
  "kidney_left",       1.4,   2.1,
  "colon",             0.0,   1.0,
  "urinary_bladder",   0.0,  -0.4
)
# Given an igraph object g whose vertex names are organ names:

make_organ_layout <- function(g, coords = mutate(organ_coords, x = 2 * x, y = 2 * y)) {
  v <- tibble(name = igraph::V(g)$name)
  
  layout_tbl <- v %>%
    left_join(coords, by = "name")
  
  # Optional fallback for unknown nodes
  missing <- which(is.na(layout_tbl$x) | is.na(layout_tbl$y))
  if (length(missing) > 0) {
    fallback <- layout_with_fr(g)
    layout_tbl$x[missing] <- fallback[missing, 1]
    layout_tbl$y[missing] <- fallback[missing, 2]
  }
  
  as.matrix(layout_tbl[, c("x", "y")])
}




layout_fixed <- make_organ_layout(g_st_no)
# plot(g_st_isch, layout=layout_fixed, edge.width =w_scaled_stisch)
plot(g_st_isch, layout=layout_fixed,
     vertex.size = 20,
     edge.width = 2,
     vertex.label.cex = 0.8,
     edge.color = rev(gray.colors(105))[res_stress_isch$rank_skew])

mtext("Stress Ischemia", side = 1, line = 0)

res_stress_isch[order(res_stress_isch$rank_skew, decreasing = TRUE)[1:10],]
xtable::xtable(res_stress_isch[order(res_stress_isch$rank_skew, decreasing = TRUE)[1:10],])




#layout_fixed <- layout_with_fr(g_st_no)
# plot(g_st_no, layout=layout_fixed, edge.width = w_scaled_stno)
plot(g_st_no, layout=layout_fixed,
     vertex.size = 20,
     edge.width = 2,
     vertex.label.cex = 0.8,
     edge.color = rev(gray.colors(105))[res_stress_no$rank_skew])

mtext("Stress No Ischemia", side = 1, line = 0)

res_stress_no[order(res_stress_no$rank_skew, decreasing = TRUE)[1:10],]
xtable::xtable(res_stress_no[order(res_stress_no$rank_skew, decreasing = TRUE)[1:10],])

#layout_fixed <- layout_with_fr(g_rest_isch)
# plot(g_st_isch, layout=layout_fixed, edge.width =w_scaled_stisch)
plot(g_rest_isch, layout=layout_fixed,
     vertex.size = 20,
     edge.width = 2,
     vertex.label.cex = 0.8,
     edge.color = rev(gray.colors(105))[res_rest_isch$rank_skew])

mtext("Rest Ischemia", side = 1, line = 0)

res_rest_isch[order(res_rest_isch$rank_skew, decreasing = TRUE)[1:10],]
xtable::xtable(res_rest_isch[order(res_rest_isch$rank_skew, decreasing = TRUE)[1:10],])



#layout_fixed <- layout_with_fr(g_rest_no)
# plot(g_st_no, layout=layout_fixed, edge.width = w_scaled_stno)
plot(g_rest_no, layout=layout_fixed,
     vertex.size = 20,
     edge.width = 2,
     vertex.label.cex = 0.8,
     edge.color = rev(gray.colors(105))[res_rest_no$rank_skew])

mtext("Rest No Ischemia", side = 1, line = 0)


res_rest_no[order(res_rest_no$rank_skew, decreasing = TRUE)[1:10],]
xtable::xtable(res_rest_no[order(res_rest_no$rank_skew, decreasing = TRUE)[1:10],])


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

