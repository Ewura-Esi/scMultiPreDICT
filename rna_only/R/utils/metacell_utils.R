# ============================================================================
# scMultiPreDICT - Metacell Utilities (RNA-only Pipeline)
# ============================================================================
# Reusable functions for creating metacells using PCA or scVI embeddings
#
# Usage:
#   source("R/utils/metacell_utils.R")
# ============================================================================

#' Build PCA dimensionality reduction (fitted on training data only)
#'
#' @param seurat_obj Seurat object with data_split column
#' @param n_pcs Number of PCs
#' @param n_hvg Number of HVGs to use
#' @return List with pca_rotation and hvg_genes
#' @export
build_dimension_reduction <- function(seurat_obj, n_pcs = 50, n_hvg = 3000) {
  
  cat("=== Building Dimension Reduction ===\n\n")
  
  # Get training cells
  train_cells <- colnames(seurat_obj)[seurat_obj$data_split == "train"]
  train_obj <- subset(seurat_obj, cells = train_cells)
  
  cat(sprintf("Training on %d cells\n", length(train_cells)))
  
  # RNA: Find HVGs and run PCA on training data
  cat("\n--- RNA PCA ---\n")
  Seurat::DefaultAssay(train_obj) <- "RNA"
  train_obj <- Seurat::NormalizeData(train_obj, verbose = FALSE)
  train_obj <- Seurat::FindVariableFeatures(train_obj, nfeatures = n_hvg, verbose = FALSE)
  hvg_genes <- Seurat::VariableFeatures(train_obj)
  
  train_obj <- Seurat::ScaleData(train_obj, features = hvg_genes, verbose = FALSE)
  train_obj <- Seurat::RunPCA(train_obj, features = hvg_genes, npcs = n_pcs, verbose = FALSE)
  
  # Extract PCA rotation matrix
  pca_rotation <- Seurat::Loadings(train_obj, "pca")
  pca_center <- train_obj@reductions$pca@misc$center
  pca_scale <- train_obj@reductions$pca@misc$scale
  
  cat(sprintf("  HVGs: %d\n", length(hvg_genes)))
  cat(sprintf("  PCs: %d\n", n_pcs))
  
  result <- list(
    pca_rotation = pca_rotation,
    pca_center = pca_center,
    pca_scale = pca_scale,
    hvg_genes = hvg_genes,
    n_pcs = n_pcs
  )
  
  cat("\n✓ Dimension reduction complete\n")
  return(result)
}


#' Project cells into PCA space
#'
#' @param seurat_obj Seurat object
#' @param dim_reduction Output from build_dimension_reduction()
#' @return Matrix of PCA embeddings (cells × features)
#' @export
project_to_pca_space <- function(seurat_obj, dim_reduction) {
  
  cat("Projecting cells to PCA space...\n")
  
  # Get RNA data and project to PCA space
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  rna_data <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  
  # Subset to HVGs and scale
  hvgs <- dim_reduction$hvg_genes
  rna_hvg <- as.matrix(t(rna_data[hvgs, , drop = FALSE]))
  
  # Center and scale using training parameters
  if (!is.null(dim_reduction$pca_center)) {
    rna_hvg <- scale(rna_hvg, center = dim_reduction$pca_center, 
                     scale = dim_reduction$pca_scale)
  }
  
  # Project to PCA space
  pca_coords <- rna_hvg %*% dim_reduction$pca_rotation
  
  cat(sprintf("  PCA projection: %d cells × %d dims\n", 
              nrow(pca_coords), ncol(pca_coords)))
  
  return(pca_coords)
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
  cat("Aggregating expression values...\n")
  
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


#' Create metacells using scVI latent space
#'
#' @param scvi_latent Matrix of scVI latent embeddings
#' @param expression_matrix Expression matrix to aggregate
#' @param k Number of neighbors
#' @return Aggregated expression matrix
#' @export
make_metacells_scvi <- function(scvi_latent, expression_matrix, k = 50) {
  
  cat("=== Creating scVI Metacells ===\n\n")
  cat(sprintf("scVI latent space: %d dims\n", ncol(scvi_latent)))
  
  # Normalize latent space
  scvi_scaled <- scale(scvi_latent)
  
  # Use the standard KNN function
  metacell_matrix <- make_metacells_knn(
    expression_matrix = expression_matrix,
    embedding_matrix = scvi_scaled,
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
