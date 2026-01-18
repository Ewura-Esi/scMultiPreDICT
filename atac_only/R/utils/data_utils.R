# ============================================================================
# scMultiPreDICT - Data Utilities (ATAC-only Pipeline)
# ============================================================================
# Reusable functions for ATAC-only data loading, splitting, and preprocessing
#
# Usage:
#   source("R/utils/data_utils.R")
# ============================================================================

#' Load 10X Genomics scATAC-seq data into a Seurat object
#'
#' @param mtx_path Path to matrix.mtx.gz file (peaks × cells)
#' @param features_path Path to features.tsv.gz file (peaks)
#' @param barcodes_path Path to barcodes.tsv.gz file
#' @param fragments_path Path to fragments.tsv.gz file
#' @param genome Genome assembly ("mm10" or "hg38")
#' @param project_name Name for the Seurat project
#' @return Seurat object with ATAC assay
#' @export
load_atac_data <- function(mtx_path, features_path, barcodes_path,
                           fragments_path, genome = "mm10",
                           project_name = "scatac") {
  
  cat("=== Loading ATAC Data ===\n\n")
  
  # Validate inputs
  for (f in c(mtx_path, features_path, barcodes_path)) {
    if (!file.exists(f)) stop(paste("File not found:", f))
  }
  
  # Load count matrix
  cat("Loading peak count matrix...\n")
  counts <- Seurat::ReadMtx(
    mtx = mtx_path,
    features = features_path,
    cells = barcodes_path
  )
  cat(sprintf("  Loaded: %d peaks × %d cells\n", nrow(counts), ncol(counts)))
  
  # Create fragment object
  cat("\nCreating fragment object...\n")
  if (!file.exists(fragments_path)) {
    stop(paste("Fragments file not found:", fragments_path))
  }
  
  frags <- Signac::CreateFragmentObject(
    path = fragments_path,
    cells = colnames(counts),
    validate.fragments = TRUE
  )
  
  # Create chromatin assay
  cat("Creating chromatin assay...\n")
  chrom_assay <- Signac::CreateChromatinAssay(
    counts = counts,
    fragments = frags,
    genome = genome,
    sep = c(":", "-")
  )
  
  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    project = project_name
  )
  
  cat(sprintf("\n✓ Created Seurat object: %d cells, %d peaks\n", 
              ncol(seurat_obj), nrow(seurat_obj)))
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


#' Apply quality control filters to ATAC data
#'
#' @param seurat_obj Seurat object with ATAC assay
#' @param min_count_atac Minimum ATAC fragments per cell
#' @param max_count_atac Maximum ATAC fragments per cell
#' @param min_tss_enrichment Minimum TSS enrichment score
#' @param max_nucleosome_signal Maximum nucleosome signal
#' @return Filtered Seurat object
#' @export
apply_qc_filters <- function(seurat_obj,
                              min_count_atac = 1000,
                              max_count_atac = 100000,
                              min_tss_enrichment = 3,
                              max_nucleosome_signal = 1.7) {
  
  cat("=== Applying QC Filters ===\n\n")
  n_before <- ncol(seurat_obj)
  
  # Compute ATAC QC metrics
  cat("Computing ATAC QC metrics...\n")
  Seurat::DefaultAssay(seurat_obj) <- "ATAC"
  seurat_obj <- Signac::NucleosomeSignal(seurat_obj)
  seurat_obj <- Signac::TSSEnrichment(seurat_obj)
  
  # Apply filters
  cat("\nApplying filters...\n")
  seurat_obj <- subset(
    x = seurat_obj,
    subset = nCount_ATAC > min_count_atac &
      nCount_ATAC < max_count_atac &
      TSS.enrichment > min_tss_enrichment &
      nucleosome_signal < max_nucleosome_signal
  )
  
  n_after <- ncol(seurat_obj)
  
  cat(sprintf("\nCells before: %d\n", n_before))
  cat(sprintf("Cells after:  %d\n", n_after))
  cat(sprintf("Removed:      %d (%.1f%%)\n", n_before - n_after, 
              100 * (n_before - n_after) / n_before))
  
  cat("\n✓ QC filtering complete\n")
  return(seurat_obj)
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


#' Filter low-accessibility peaks from Seurat object
#'
#' @param seurat_obj Seurat object with ATAC assay
#' @param min_fraction Minimum fraction of cells with peak accessible (default: 0.01 = 1%)
#' @return Seurat object with filtered peaks
#' @export
filter_low_accessibility_peaks <- function(seurat_obj, min_fraction = 0.01) {
  
  cat("=== Filtering Low-Accessibility Peaks ===\n\n")
  
  Seurat::DefaultAssay(seurat_obj) <- "ATAC"
  
  n_cells <- ncol(seurat_obj)
  n_peaks_before <- nrow(seurat_obj[["ATAC"]])
  
  # Calculate fraction of cells each peak is accessible in
  peak_counts <- Seurat::GetAssayData(seurat_obj, assay = "ATAC", layer = "counts")
  peak_frac <- Matrix::rowSums(peak_counts > 0) / n_cells
  
  # Keep peaks accessible in at least min_fraction of cells
  peaks_keep <- names(peak_frac[peak_frac >= min_fraction])
  
  # Subset the ATAC assay
  seurat_obj[["ATAC"]] <- subset(seurat_obj[["ATAC"]], features = peaks_keep)
  
  n_peaks_after <- nrow(seurat_obj[["ATAC"]])
  
  cat(sprintf("Filter: Keep peaks accessible in >= %.1f%% of cells\n", min_fraction * 100))
  cat(sprintf("Peaks before: %d\n", n_peaks_before))
  cat(sprintf("Peaks after:  %d\n", n_peaks_after))
  cat(sprintf("Peaks removed: %d (%.1f%%)\n",
              n_peaks_before - n_peaks_after,
              100 * (n_peaks_before - n_peaks_after) / n_peaks_before))
  
  cat("\n✓ Peak filtering complete\n")
  return(seurat_obj)
}
