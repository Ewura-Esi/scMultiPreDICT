# ============================================================================
# scMultiPreDICT - Data Utilities
# ============================================================================
# Reusable functions for data loading, splitting, and preprocessing
#
# Usage:
#   source("R/utils/data_utils.R")
# ============================================================================

#' Load 10X Genomics multiome data into a Seurat object
#'
#' @param mtx_path Path to matrix.mtx.gz file
#' @param features_path Path to features.tsv.gz file
#' @param barcodes_path Path to barcodes.tsv.gz file
#' @param fragments_path Path to fragments.tsv.gz file (ATAC)
#' @param genome Genome assembly ("mm10" or "hg38")
#' @param project_name Name for the Seurat project
#' @return Seurat object with RNA and ATAC assays
#' @export
load_multiome_data <- function(mtx_path, features_path, barcodes_path,
                                fragments_path, genome = "mm10",
                                project_name = "multiome") {
  
  cat("=== Loading Multiome Data ===\n\n")
  
  # Validate inputs
  for (f in c(mtx_path, features_path, barcodes_path)) {
    if (!file.exists(f)) stop(paste("File not found:", f))
  }
  
  # Load count matrix
  cat("Loading count matrix...\n")
  counts_all <- Seurat::ReadMtx(
    mtx = mtx_path,
    features = features_path,
    cells = barcodes_path
  )
  cat(sprintf("  Loaded: %d features × %d cells\n", nrow(counts_all), ncol(counts_all)))
  
  # Load features to split RNA and ATAC
  feat <- readr::read_tsv(features_path, col_names = FALSE, show_col_types = FALSE)
  
  # Split by feature type
  rna_rows <- which(feat$X3 %in% c("Gene Expression", "Gene expression", "Gene"))
  atac_rows <- which(feat$X3 %in% c("Peaks", "Peak"))
  
  cat(sprintf("  RNA features: %d\n", length(rna_rows)))
  cat(sprintf("  ATAC features: %d\n", length(atac_rows)))
  
  if (length(rna_rows) == 0) {
    stop("No RNA features found. Check column 3 of features.tsv")
  }
  
  rna_counts <- counts_all[rna_rows, , drop = FALSE]
  atac_counts <- counts_all[atac_rows, , drop = FALSE]
  
  # Create Seurat object
  cat("\nCreating Seurat object...\n")
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA",
    project = project_name
  )
  
  # Add ATAC assay if fragments provided
  if (!is.null(fragments_path) && file.exists(fragments_path)) {
    cat("Adding ATAC assay...\n")
    
    # Create fragment object
    frags <- Signac::CreateFragmentObject(
      path = fragments_path,
      cells = colnames(seurat_obj),
      validate.fragments = TRUE
    )
    
    # Create chromatin assay
    seurat_obj[["ATAC"]] <- Signac::CreateChromatinAssay(
      counts = atac_counts,
      fragments = frags,
      genome = genome,
      sep = c(":", "-")
    )
    
    cat(sprintf("  Added ATAC assay: %d peaks\n", nrow(seurat_obj[["ATAC"]])))
  }
  
  cat(sprintf("\n✓ Created Seurat object: %d cells\n", ncol(seurat_obj)))
  return(seurat_obj)
}


#' Split data into train/validation/test sets
#'
#' @param seurat_obj Seurat object
#' @param train_frac Fraction for training (default: 0.70)
#' @param val_frac Fraction for validation (default: 0.20)
#' @param test_frac Fraction for test (default: 0.10)
#' @param seed Random seed for reproducibility
#' @return Seurat object with 'data_split' column in metadata
#' @export
split_data <- function(seurat_obj, train_frac = 0.70, val_frac = 0.20,
                       test_frac = 0.10, seed = 42) {
  
  cat("=== Splitting Data ===\n\n")
  
  # Validate proportions
  total <- train_frac + val_frac + test_frac
  if (abs(total - 1.0) > 0.001) {
    stop(sprintf("Split fractions must sum to 1.0 (got %.3f)", total))
  }
  
  # Get cell info
  all_cells <- colnames(seurat_obj)
  n_cells <- length(all_cells)
  
  # Set seed and shuffle
  set.seed(seed)
  shuffled_idx <- sample(1:n_cells, size = n_cells, replace = FALSE)
  
  # Calculate split sizes
  n_train <- floor(train_frac * n_cells)
  n_val <- floor(val_frac * n_cells)
  n_test <- n_cells - n_train - n_val
  
  # Assign splits
  train_idx <- shuffled_idx[1:n_train]
  val_idx <- shuffled_idx[(n_train + 1):(n_train + n_val)]
  test_idx <- shuffled_idx[(n_train + n_val + 1):n_cells]
  
  # Add to metadata
  seurat_obj$data_split <- NA
  seurat_obj$data_split[train_idx] <- "train"
  seurat_obj$data_split[val_idx] <- "validation"
  seurat_obj$data_split[test_idx] <- "test"
  
  # Print summary
  cat("Split configuration:\n")
  cat(sprintf("  Random seed: %d\n", seed))
  cat(sprintf("  Train:       %5d cells (%5.1f%%)\n", n_train, 100 * n_train / n_cells))
  cat(sprintf("  Validation:  %5d cells (%5.1f%%)\n", n_val, 100 * n_val / n_cells))
  cat(sprintf("  Test:        %5d cells (%5.1f%%)\n", n_test, 100 * n_test / n_cells))
  cat(sprintf("  Total:       %5d cells\n", n_cells))
  
  # Return split info as well
  attr(seurat_obj, "split_info") <- list(
    train_idx = train_idx,
    val_idx = val_idx,
    test_idx = test_idx,
    seed = seed
  )
  
  cat("\n✓ Data split complete\n")
  return(seurat_obj)
}


#' Apply quality control filters to multiome data
#'
#' @param seurat_obj Seurat object with RNA and ATAC assays
#' @param min_features_rna Minimum genes per cell
#' @param max_features_rna Maximum genes per cell (doublet filter)
#' @param max_percent_mt Maximum mitochondrial percentage
#' @param min_count_atac Minimum ATAC fragments per cell
#' @param max_count_atac Maximum ATAC fragments per cell
#' @param min_tss_enrichment Minimum TSS enrichment score
#' @param max_nucleosome_signal Maximum nucleosome signal
#' @param mt_pattern Pattern for mitochondrial genes ("^mt-" for mouse, "^MT-" for human)
#' @return Filtered Seurat object
#' @export
apply_qc_filters <- function(seurat_obj,
                              min_features_rna = 1000,
                              max_features_rna = 6000,
                              max_percent_mt = 35,
                              min_count_rna = 1000,
                              max_count_rna = 40000,
                              min_count_atac = 1000,
                              max_count_atac = 100000,
                              min_tss_enrichment = 3,
                              max_nucleosome_signal = 1.7,
                              mt_pattern = "^mt-") {
  
  cat("=== Applying QC Filters ===\n\n")
  n_before <- ncol(seurat_obj)
  
  # Compute RNA QC metrics
  cat("Computing RNA metrics...\n")
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  
  # Compute ATAC QC metrics if ATAC assay exists
  if ("ATAC" %in% names(seurat_obj@assays)) {
    cat("Computing ATAC metrics...\n")
    Seurat::DefaultAssay(seurat_obj) <- "ATAC"
    seurat_obj <- Signac::NucleosomeSignal(seurat_obj)
    seurat_obj <- Signac::TSSEnrichment(seurat_obj)
    
    # Apply combined filters (matching Step_031 exactly)
    cat("\nApplying filters...\n")
    seurat_obj <- subset(
      x = seurat_obj,
      subset = nFeature_RNA > min_features_rna &
        nFeature_RNA < max_features_rna &
        percent.mt < max_percent_mt &
        nCount_RNA > min_count_rna &
        nCount_RNA < max_count_rna &
        nCount_ATAC > min_count_atac &
        nCount_ATAC < max_count_atac &
        TSS.enrichment > min_tss_enrichment &
        nucleosome_signal < max_nucleosome_signal
    )
  } else {
    # RNA-only filters
    seurat_obj <- subset(
      x = seurat_obj,
      subset = nFeature_RNA > min_features_rna &
        nFeature_RNA < max_features_rna &
        percent.mt < max_percent_mt &
        nCount_RNA > min_count_rna &
        nCount_RNA < max_count_rna
    )
  }
  
  n_after <- ncol(seurat_obj)
  
  cat(sprintf("\nCells before: %d\n", n_before))
  cat(sprintf("Cells after:  %d\n", n_after))
  cat(sprintf("Removed:      %d (%.1f%%)\n", n_before - n_after, 
              100 * (n_before - n_after) / n_before))
  
  cat("\n✓ QC filtering complete\n")
  return(seurat_obj)
}


#' Select target genes for prediction
#'
#' @param seurat_obj Seurat object (training cells only recommended)
#' @param n_hvg_targets Number of HVG target genes
#' @param n_random_targets Number of random (non-HVG) target genes
#' @param min_detection Minimum detection rate (% of cells)
#' @param max_detection Maximum detection rate (% of cells)
#' @param min_mean_expr Minimum mean expression
#' @param seed Random seed
#' @return List with hvg_genes and random_genes
#' @export
select_target_genes <- function(seurat_obj, n_hvg_targets = 100,
                                 n_random_targets = 100,
                                 min_detection = 5, max_detection = 95,
                                 min_mean_expr = 0.1, seed = 42) {
  
  cat("=== Selecting Target Genes ===\n\n")
  set.seed(seed)
  
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  
  # Find HVGs
  cat("Finding highly variable genes...\n")
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, nfeatures = 3000, verbose = FALSE)
  all_hvgs <- Seurat::VariableFeatures(seurat_obj)
  
  cat(sprintf("  Found %d HVGs\n", length(all_hvgs)))
  
  # Sample HVG targets
  hvg_targets <- sample(all_hvgs, min(n_hvg_targets, length(all_hvgs)))
  cat(sprintf("  Selected %d HVG targets\n", length(hvg_targets)))
  
  # Find random (non-HVG) targets
  cat("\nFinding random target genes...\n")
  all_genes <- rownames(seurat_obj[["RNA"]])
  non_hvgs <- setdiff(all_genes, all_hvgs)
  
  # Compute statistics
  rna_counts <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  rna_counts_non_hvg <- rna_counts[non_hvgs, , drop = FALSE]
  
  detection_rate <- Matrix::rowMeans(rna_counts_non_hvg > 0) * 100
  mean_expression <- Matrix::rowMeans(rna_counts_non_hvg)
  
  # Filter by criteria
  eligible <- non_hvgs[
    detection_rate >= min_detection &
      detection_rate <= max_detection &
      mean_expression >= min_mean_expr
  ]
  
  cat(sprintf("  Eligible non-HVG genes: %d\n", length(eligible)))
  
  # Sample random targets
  random_targets <- sample(eligible, min(n_random_targets, length(eligible)))
  cat(sprintf("  Selected %d random targets\n", length(random_targets)))
  
  result <- list(
    hvg_genes = hvg_targets,
    random_genes = random_targets,
    all_hvgs = all_hvgs
  )
  
  cat("\n✓ Target gene selection complete\n")
  return(result)
}


#' Save split information and cell assignments
#'
#' @param seurat_obj Seurat object with data_split column
#' @param output_dir Output directory
#' @param sample_name Sample identifier
#' @export
save_split_info <- function(seurat_obj, output_dir, sample_name) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save cell-to-split mapping
  split_df <- data.frame(
    cell_barcode = colnames(seurat_obj),
    data_split = seurat_obj$data_split,
    stringsAsFactors = FALSE
  )
  
  csv_file <- file.path(output_dir, paste0(sample_name, "_cell_splits.csv"))
  write.csv(split_df, csv_file, row.names = FALSE)
  cat(sprintf("Saved: %s\n", csv_file))
  
  # Save split indices
  split_info <- attr(seurat_obj, "split_info")
  if (!is.null(split_info)) {
    rds_file <- file.path(output_dir, paste0(sample_name, "_split_indices.rds"))
    saveRDS(split_info, rds_file)
    cat(sprintf("Saved: %s\n", rds_file))
  }
}


#' Filter low-expression genes from Seurat object
#'
#' @param seurat_obj Seurat object with RNA assay
#' @param min_fraction Minimum fraction of cells expressing gene (default: 0.10 = 10%)
#' @return Seurat object with filtered genes
#' @export
filter_low_expression_genes <- function(seurat_obj, min_fraction = 0.10) {
  
  cat("=== Filtering Low-Expression Genes ===\n\n")
  
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  
  n_cells <- ncol(seurat_obj)
  n_genes_before <- nrow(seurat_obj[["RNA"]])
  
  # Calculate fraction of cells each gene is expressed in
  gene_counts <- Seurat::GetAssayData(seurat_obj, layer = "counts")
  gene_frac <- Matrix::rowSums(gene_counts > 0) / n_cells
  
  # Keep genes expressed in at least min_fraction of cells
  genes_keep <- names(gene_frac[gene_frac >= min_fraction])
  
  # Subset only the RNA assay
  seurat_obj[["RNA"]] <- subset(seurat_obj[["RNA"]], features = genes_keep)
  
  n_genes_after <- nrow(seurat_obj[["RNA"]])
  
  cat(sprintf("Filter: Keep genes expressed in >= %.0f%% of cells\n", min_fraction * 100))
  cat(sprintf("Genes before: %d\n", n_genes_before))
  cat(sprintf("Genes after:  %d\n", n_genes_after))
  cat(sprintf("Genes removed: %d (%.1f%%)\n",
              n_genes_before - n_genes_after,
              100 * (n_genes_before - n_genes_after) / n_genes_before))
  
  cat("\n✓ Gene filtering complete\n")
  return(seurat_obj)
}
