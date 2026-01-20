# ============================================================================
# scMultiPreDICT - Feature Extraction Utilities (RNA-only Pipeline)
# ============================================================================
# Reusable functions for extracting RNA features for prediction
#
# Usage:
#   source("R/utils/feature_utils.R")
# ============================================================================

#' Extract RNA features (HVGs) for prediction
#'
#' @param seurat_obj Seurat object
#' @param n_hvg Number of highly variable genes to use
#' @param exclude_genes Genes to exclude (target genes)
#' @param use_metacells Whether to return metacell-aggregated values
#' @param metacell_matrix Metacell expression matrix (if use_metacells = TRUE)
#' @return Matrix of RNA features (cells × features)
#' @export
extract_rna_features <- function(seurat_obj, n_hvg = 1000, exclude_genes = NULL,
                                  use_metacells = FALSE, metacell_matrix = NULL) {
  
  cat("=== Extracting RNA Features ===\n\n")
  
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  
  # Get HVGs (should be pre-computed from training data)
  hvgs <- Seurat::VariableFeatures(seurat_obj)
  
  if (length(hvgs) == 0) {
    cat("Finding HVGs...\n")
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, nfeatures = n_hvg * 2,
                                                verbose = FALSE)
    hvgs <- Seurat::VariableFeatures(seurat_obj)
  }
  
  # Exclude target genes from features
  if (!is.null(exclude_genes)) {
    hvgs <- setdiff(hvgs, exclude_genes)
    cat(sprintf("  Excluded %d target genes from features\n", 
                length(exclude_genes) - length(setdiff(exclude_genes, hvgs))))
  }
  
  # Take top n_hvg
  hvgs <- head(hvgs, n_hvg)
  cat(sprintf("  Using %d HVGs as features\n", length(hvgs)))
  
  if (use_metacells && !is.null(metacell_matrix)) {
    # Use metacell values
    available_hvgs <- intersect(hvgs, colnames(metacell_matrix))
    feature_matrix <- metacell_matrix[, available_hvgs, drop = FALSE]
    cat(sprintf("  Using metacell-aggregated values\n"))
  } else {
    # Use raw/normalized values
    rna_data <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data")
    feature_matrix <- t(as.matrix(rna_data[hvgs, , drop = FALSE]))
  }
  
  cat(sprintf("  RNA feature matrix: %d cells × %d features\n", 
              nrow(feature_matrix), ncol(feature_matrix)))
  
  return(feature_matrix)
}


#' Extract features for all target genes (RNA-only)
#'
#' @param seurat_obj Seurat object
#' @param target_genes Vector of target gene names
#' @param n_hvg Number of HVGs for RNA features
#' @param metacell_rna Metacell RNA matrix (optional)
#' @return List with feature matrices and target values for each gene
#' @export
extract_features <- function(seurat_obj, target_genes, n_hvg = 1000,
                              metacell_rna = NULL) {
  
  cat("=== Extracting RNA-only Features ===\n\n")
  cat(sprintf("Target genes: %d\n", length(target_genes)))
  
  # Get target gene expression
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  rna_data <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  
  # Extract RNA features (HVGs, excluding target genes)
  use_mc_rna <- !is.null(metacell_rna)
  rna_features <- extract_rna_features(
    seurat_obj = seurat_obj,
    n_hvg = n_hvg,
    exclude_genes = target_genes,
    use_metacells = use_mc_rna,
    metacell_matrix = metacell_rna
  )
  
  # Process each target gene
  results <- list()
  n_genes <- length(target_genes)
  
  pb <- txtProgressBar(min = 0, max = n_genes, style = 3)
  
  for (i in seq_along(target_genes)) {
    gene <- target_genes[i]
    
    # Skip if gene not in expression matrix
    if (!gene %in% rownames(rna_data)) {
      setTxtProgressBar(pb, i)
      next
    }
    
    # Get target values (y)
    y <- as.numeric(rna_data[gene, ])
    
    # Store results (RNA features only)
    results[[gene]] <- list(
      X = rna_features,
      y = y,
      n_rna_features = ncol(rna_features),
      n_atac_features = 0
    )
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  cat(sprintf("\n✓ Extracted RNA features for %d genes\n", length(results)))
  
  return(results)
}


#' Save extracted features to disk
#'
#' @param features_list Output from extract_features()
#' @param output_dir Output directory
#' @param filename Output filename
#' @export
save_features <- function(features_list, output_dir, filename = "gene_specific_features.rds") {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(output_dir, filename)
  saveRDS(features_list, output_file)
  
  cat(sprintf("Saved: %s (%d genes)\n", output_file, length(features_list)))
}
