# ============================================================================
# scMultiPreDICT - Metacell Utilities
# ============================================================================
# Reusable functions for creating metacells using different approaches
#
# Usage:
#   source("R/utils/metacell_utils.R")
# ============================================================================

#' Build PCA/LSI dimensionality reduction (fitted on training data only)
#'
#' @param seurat_obj Seurat object with data_split column
#' @param n_pcs Number of PCs for RNA
#' @param n_lsi Number of LSI components for ATAC
#' @param n_hvg Number of HVGs to use
#' @return List with pca_rotation, lsi_components, and hvg_genes
#' @export
build_dimension_reduction <- function(seurat_obj, n_pcs = 50, n_lsi = 50,
                                       n_hvg = 3000) {
  
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
  
  # ATAC: Run TF-IDF and LSI on training data
  cat("\n--- ATAC LSI ---\n")
  if ("ATAC" %in% names(seurat_obj@assays)) {
    Seurat::DefaultAssay(train_obj) <- "ATAC"
    train_obj <- Signac::RunTFIDF(train_obj, verbose = FALSE)
    train_obj <- Signac::FindTopFeatures(train_obj, min.cutoff = "q0", verbose = FALSE)
    train_obj <- Signac::RunSVD(train_obj, n = n_lsi, verbose = FALSE)
    
    # Extract LSI loadings
    lsi_loadings <- train_obj@reductions$lsi@feature.loadings
    lsi_sdev <- train_obj@reductions$lsi@stdev
    
    cat(sprintf("  LSI components: %d\n", n_lsi))
  } else {
    lsi_loadings <- NULL
    lsi_sdev <- NULL
    cat("  No ATAC assay found - skipping LSI\n")
  }
  
  result <- list(
    pca_rotation = pca_rotation,
    pca_center = pca_center,
    pca_scale = pca_scale,
    hvg_genes = hvg_genes,
    lsi_loadings = lsi_loadings,
    lsi_sdev = lsi_sdev,
    n_pcs = n_pcs,
    n_lsi = n_lsi
  )
  
  cat("\n✓ Dimension reduction complete\n")
  return(result)
}


#' Project cells into joint PCA/LSI space
#'
#' @param seurat_obj Seurat object
#' @param dim_reduction Output from build_dimension_reduction()
#' @return Matrix of joint embeddings (cells × features)
#' @export
project_to_joint_space <- function(seurat_obj, dim_reduction) {
  
  cat("Projecting cells to joint space...\n")
  
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
  
  # Get ATAC projection if available
  if (!is.null(dim_reduction$lsi_loadings) && "ATAC" %in% names(seurat_obj@assays)) {
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
    
    # Combine (skip LSI component 1 which captures depth)
    joint_coords <- cbind(pca_coords, lsi_coords[, 2:ncol(lsi_coords)])
  } else {
    joint_coords <- pca_coords
  }
  
  cat(sprintf("  Joint space: %d cells × %d dims\n", 
              nrow(joint_coords), ncol(joint_coords)))
  
  return(joint_coords)
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
  
  # Process in chunks for memory efficiency
  chunk_size <- 1000
  n_chunks <- ceiling(n_cells / chunk_size)
  
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


#' Create metacells using WNN (Weighted Nearest Neighbor) approach
#'
#' @param seurat_obj Seurat object with RNA and ATAC
#' @param expression_matrix Expression matrix to aggregate
#' @param k Number of neighbors
#' @param weight_rna Weight for RNA modality (0-1)
#' @param weight_atac Weight for ATAC modality (0-1)
#' @return Aggregated expression matrix
#' @export
make_metacells_wnn <- function(seurat_obj, expression_matrix, k = 50,
                                weight_rna = 0.5, weight_atac = 0.5) {
  
  cat("=== Creating WNN Metacells ===\n\n")
  
  # Run WNN analysis
  cat("Running WNN analysis...\n")
  seurat_obj <- Seurat::FindMultiModalNeighbors(
    seurat_obj,
    reduction.list = list("pca", "lsi"),
    dims.list = list(1:50, 2:50),
    modality.weight.name = "wnn.weight",
    verbose = FALSE
  )
  
  # Get WNN graph
  wnn_graph <- seurat_obj@graphs$wnn
  
  # For each cell, find top k neighbors
  cat("Finding WNN neighbors...\n")
  n_cells <- nrow(expression_matrix)
  metacell_matrix <- matrix(0, nrow = n_cells, ncol = ncol(expression_matrix))
  colnames(metacell_matrix) <- colnames(expression_matrix)
  rownames(metacell_matrix) <- rownames(expression_matrix)
  
  for (i in 1:n_cells) {
    weights <- wnn_graph[i, ]
    top_k <- order(weights, decreasing = TRUE)[1:k]
    neighbor_expr <- expression_matrix[top_k, , drop = FALSE]
    metacell_matrix[i, ] <- colMeans(neighbor_expr)
    
    if (i %% 500 == 0) cat(sprintf("  Processed %d/%d cells\n", i, n_cells))
  }
  
  cat(sprintf("\n✓ Created WNN metacells: %d cells × %d features\n", 
              nrow(metacell_matrix), ncol(metacell_matrix)))
  
  return(metacell_matrix)
}


#' Create metacells using scVI/PeakVI latent space
#'
#' @param rna_latent Matrix of scVI latent embeddings
#' @param atac_latent Matrix of PeakVI latent embeddings
#' @param expression_matrix Expression matrix to aggregate
#' @param k Number of neighbors
#' @param weight_rna Weight for RNA latent space
#' @param weight_atac Weight for ATAC latent space
#' @return Aggregated expression matrix
#' @export
make_metacells_scvi <- function(rna_latent, atac_latent = NULL, 
                                 expression_matrix, k = 50,
                                 weight_rna = 0.5, weight_atac = 0.5) {
  
  cat("=== Creating scVI/PeakVI Metacells ===\n\n")
  
  # Combine latent spaces if ATAC provided
  if (!is.null(atac_latent)) {
    # Normalize each latent space
    rna_scaled <- scale(rna_latent)
    atac_scaled <- scale(atac_latent)
    
    # Weighted combination
    joint_latent <- cbind(weight_rna * rna_scaled, weight_atac * atac_scaled)
    cat(sprintf("Combined latent space: %d dims\n", ncol(joint_latent)))
  } else {
    joint_latent <- scale(rna_latent)
    cat(sprintf("RNA-only latent space: %d dims\n", ncol(joint_latent)))
  }
  
  # Use the standard KNN function
  metacell_matrix <- make_metacells_knn(
    expression_matrix = expression_matrix,
    embedding_matrix = joint_latent,
    k = k
  )
  
  return(metacell_matrix)
}


#' Create metacells using MultiVI joint latent space
#'
#' @param multivi_latent Matrix of MultiVI latent embeddings
#' @param expression_matrix Expression matrix to aggregate
#' @param k Number of neighbors
#' @return Aggregated expression matrix
#' @export
make_metacells_multivi <- function(multivi_latent, expression_matrix, k = 50) {
  
  cat("=== Creating MultiVI Metacells ===\n\n")
  cat(sprintf("MultiVI latent space: %d dims\n", ncol(multivi_latent)))
  
  # Use the standard KNN function
  metacell_matrix <- make_metacells_knn(
    expression_matrix = expression_matrix,
    embedding_matrix = multivi_latent,
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
