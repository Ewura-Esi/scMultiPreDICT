# ============================================================================
# scMultiPreDICT - Metacell Utilities (ATAC-only Pipeline)
# ============================================================================
# Reusable functions for creating metacells using LSI or PeakVI embeddings
#
# Usage:
#   source("R/utils/metacell_utils.R")
# ============================================================================

#' Build LSI dimensionality reduction (fitted on training data only)
#'
#' @param seurat_obj Seurat object with data_split column
#' @param n_lsi Number of LSI components
#' @return List with lsi_loadings and lsi_sdev
#' @export
build_dimension_reduction <- function(seurat_obj, n_lsi = 50) {
  
  cat("=== Building Dimension Reduction ===\n\n")
  
  # Get training cells
  train_cells <- colnames(seurat_obj)[seurat_obj$data_split == "train"]
  train_obj <- subset(seurat_obj, cells = train_cells)
  
  cat(sprintf("Training on %d cells\n", length(train_cells)))
  
  # ATAC: Run TF-IDF and LSI on training data
  cat("\n--- ATAC LSI ---\n")
  Seurat::DefaultAssay(train_obj) <- "ATAC"
  train_obj <- Signac::RunTFIDF(train_obj, verbose = FALSE)
  train_obj <- Signac::FindTopFeatures(train_obj, min.cutoff = "q0", verbose = FALSE)
  train_obj <- Signac::RunSVD(train_obj, n = n_lsi, verbose = FALSE)
  
  # Extract LSI loadings
  lsi_loadings <- train_obj@reductions$lsi@feature.loadings
  lsi_sdev <- train_obj@reductions$lsi@stdev
  
  cat(sprintf("  LSI components: %d\n", n_lsi))
  
  result <- list(
    lsi_loadings = lsi_loadings,
    lsi_sdev = lsi_sdev,
    n_lsi = n_lsi
  )
  
  cat("\n✓ Dimension reduction complete\n")
  return(result)
}


#' Project cells into LSI space
#'
#' @param seurat_obj Seurat object
#' @param dim_reduction Output from build_dimension_reduction()
#' @return Matrix of LSI embeddings (cells × features)
#' @export
project_to_lsi_space <- function(seurat_obj, dim_reduction) {
  
  cat("Projecting cells to LSI space...\n")
  
  Seurat::DefaultAssay(seurat_obj) <- "ATAC"
  
  # Run TF-IDF on new data
  seurat_obj <- Signac::RunTFIDF(seurat_obj, verbose = FALSE)
  atac_data <- Seurat::GetAssayData(seurat_obj, assay = "ATAC", layer = "data")
  
  # Project using LSI loadings
  common_features <- intersect(rownames(atac_data), rownames(dim_reduction$lsi_loadings))
  lsi_coords <- t(atac_data[common_features, ]) %*% 
    dim_reduction$lsi_loadings[common_features, ]
  
  cat(sprintf("  LSI projection: %d cells × %d dims\n", 
              nrow(lsi_coords), ncol(lsi_coords)))
  
  # Skip first component (captures sequencing depth)
  lsi_coords <- lsi_coords[, 2:ncol(lsi_coords), drop = FALSE]
  
  return(lsi_coords)
}


#' Create metacells using k-nearest neighbors aggregation
#'
#' @param expression_matrix Matrix of expression values (cells × genes)
#' @param embedding_matrix Matrix for neighbor finding (cells × dims)
#' @param k Number of neighbors for aggregation
#' @param aggregation_method "mean" or "sum"
#' @param search_method "kd" or "ann" (approximate)
#' @return Aggregated expression matrix (cells × genes)
#' @export
make_metacells_knn <- function(expression_matrix, embedding_matrix, k = 50,
                                aggregation_method = "mean", search_method = "kd") {
  
  cat("=== Creating Metacells ===\n\n")
  cat(sprintf("k = %d, method = %s, aggregation = %s\n", 
              k, search_method, aggregation_method))
  
  n_cells <- nrow(expression_matrix)
  n_features <- ncol(expression_matrix)
  
  # Ensure matrices are aligned
  if (nrow(embedding_matrix) != n_cells) {
    stop("Expression and embedding matrices have different number of cells")
  }
  
  # Find k nearest neighbors
  cat("Finding nearest neighbors...\n")
  nn_result <- RANN::nn2(
    data = embedding_matrix,
    query = embedding_matrix,
    k = k + 1,  # Include self
    treetype = search_method,
    searchtype = "standard"
  )
  
  # Remove self (first neighbor)
  nn_indices <- nn_result$nn.idx[, 2:(k + 1), drop = FALSE]
  
  # Aggregate expression values
  cat("Aggregating values...\n")
  
  # Initialize output matrix
  metacell_matrix <- matrix(0, nrow = n_cells, ncol = n_features)
  colnames(metacell_matrix) <- colnames(expression_matrix)
  rownames(metacell_matrix) <- rownames(expression_matrix)
  
  pb <- txtProgressBar(min = 0, max = n_cells, style = 3)
  
  for (i in 1:n_cells) {
    neighbor_idx <- nn_indices[i, ]
    neighbor_expr <- expression_matrix[neighbor_idx, , drop = FALSE]
    
    if (aggregation_method == "mean") {
      metacell_matrix[i, ] <- colMeans(neighbor_expr)
    } else if (aggregation_method == "sum") {
      metacell_matrix[i, ] <- colSums(neighbor_expr)
    }
    
    if (i %% 500 == 0) setTxtProgressBar(pb, i)
  }
  setTxtProgressBar(pb, n_cells)
  close(pb)
  
  cat(sprintf("\n✓ Created metacells: %d cells × %d features\n", 
              nrow(metacell_matrix), ncol(metacell_matrix)))
  
  return(metacell_matrix)
}


#' Create metacells using PeakVI latent space
#'
#' @param peakvi_latent Matrix of PeakVI latent embeddings
#' @param expression_matrix Expression matrix to aggregate
#' @param k Number of neighbors
#' @return Aggregated expression matrix
#' @export
make_metacells_peakvi <- function(peakvi_latent, expression_matrix, k = 50) {
  
  cat("=== Creating PeakVI Metacells ===\n\n")
  cat(sprintf("PeakVI latent space: %d dims\n", ncol(peakvi_latent)))
  
  # Normalize latent space
  peakvi_scaled <- scale(peakvi_latent)
  
  # Use the standard KNN function
  metacell_matrix <- make_metacells_knn(
    expression_matrix = expression_matrix,
    embedding_matrix = peakvi_scaled,
    k = k
  )
  
  return(metacell_matrix)
}


#' Save metacell results
#'
#' @param metacell_list List of metacell matrices by split
#' @param output_dir Output directory
#' @param sample_name Sample identifier
#' @param approach Metacell approach name
#' @export
save_metacells <- function(metacell_list, output_dir, sample_name, approach) {
  
  approach_dir <- file.path(output_dir, approach)
  dir.create(approach_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (split_name in names(metacell_list)) {
    mat <- metacell_list[[split_name]]
    file_path <- file.path(approach_dir, 
                           paste0(sample_name, "_", split_name, "_metacells.rds"))
    saveRDS(mat, file_path)
    cat(sprintf("Saved: %s (%d × %d)\n", file_path, nrow(mat), ncol(mat)))
  }
}
