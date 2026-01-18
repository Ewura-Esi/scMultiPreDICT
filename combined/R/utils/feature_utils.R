# ============================================================================
# scMultiPreDICT - Feature Extraction Utilities
# ============================================================================
# Reusable functions for extracting RNA and ATAC features for prediction
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


#' Get peaks within a window of a gene's TSS
#'
#' @param gene_name Gene name
#' @param peak_names Vector of peak names (chr:start-end format)
#' @param gene_annotations GRanges object with gene annotations
#' @param window_size Distance from TSS in bp (e.g., 250000 for ±250kb)
#' @return Vector of peak names within the window
#' @export
get_peaks_near_gene <- function(gene_name, peak_names, gene_annotations,
                                 window_size = 250000) {
  
  # Find gene in annotations
  gene_idx <- which(gene_annotations$gene_name == gene_name)
  
  if (length(gene_idx) == 0) {
    # Try without version numbers
    gene_annotations$gene_name_clean <- gsub("\\..*", "", gene_annotations$gene_name)
    gene_idx <- which(gene_annotations$gene_name_clean == gene_name)
  }
  
  if (length(gene_idx) == 0) {
    warning(sprintf("Gene '%s' not found in annotations", gene_name))
    return(character(0))
  }
  
  # Get TSS
  gene_gr <- gene_annotations[gene_idx[1]]
  gene_chr <- as.character(GenomicRanges::seqnames(gene_gr))
  gene_strand <- as.character(GenomicRanges::strand(gene_gr))
  
  if (gene_strand == "+") {
    tss <- GenomicRanges::start(gene_gr)
  } else {
    tss <- GenomicRanges::end(gene_gr)
  }
  
  # Define window around TSS
  window_start <- max(0, tss - window_size)
  window_end <- tss + window_size
  
  # Parse peak coordinates
  peak_info <- parse_peak_names(peak_names)
  
  # Find peaks in window
  in_window <- which(
    peak_info$chr == gene_chr &
      peak_info$end >= window_start &
      peak_info$start <= window_end
  )
  
  return(peak_names[in_window])
}


#' Parse peak names into coordinates
#'
#' @param peak_names Vector of peak names (chr:start-end format)
#' @return Data frame with chr, start, end columns
#' @export
parse_peak_names <- function(peak_names) {
  
  # Handle different separators (chr:start-end or chr-start-end)
  parts <- strsplit(peak_names, "[:-]")
  
  df <- data.frame(
    chr = sapply(parts, function(x) x[1]),
    start = as.numeric(sapply(parts, function(x) x[2])),
    end = as.numeric(sapply(parts, function(x) x[3])),
    stringsAsFactors = FALSE
  )
  rownames(df) <- peak_names
  
  return(df)
}


#' Extract ATAC features for a target gene
#'
#' @param seurat_obj Seurat object with ATAC assay
#' @param target_gene Target gene name
#' @param gene_annotations GRanges object with gene annotations
#' @param window_size Distance from TSS in bp
#' @param use_metacells Whether to use metacell-aggregated values
#' @param metacell_matrix Metacell ATAC matrix
#' @return Matrix of ATAC features for this gene (cells × peaks)
#' @export
extract_atac_features_for_gene <- function(seurat_obj, target_gene,
                                            gene_annotations, window_size = 250000,
                                            use_metacells = FALSE,
                                            metacell_matrix = NULL) {
  
  Seurat::DefaultAssay(seurat_obj) <- "ATAC"
  all_peaks <- rownames(seurat_obj[["ATAC"]])
  
  # Get peaks near this gene
  nearby_peaks <- get_peaks_near_gene(
    gene_name = target_gene,
    peak_names = all_peaks,
    gene_annotations = gene_annotations,
    window_size = window_size
  )
  
  if (length(nearby_peaks) == 0) {
    return(NULL)
  }
  
  if (use_metacells && !is.null(metacell_matrix)) {
    # Use metacell values
    available_peaks <- intersect(nearby_peaks, colnames(metacell_matrix))
    if (length(available_peaks) == 0) return(NULL)
    
    feature_matrix <- metacell_matrix[, available_peaks, drop = FALSE]
  } else {
    # Use raw ATAC counts (TF-IDF normalized)
    atac_data <- Seurat::GetAssayData(seurat_obj, assay = "ATAC", layer = "data")
    feature_matrix <- t(as.matrix(atac_data[nearby_peaks, , drop = FALSE]))
  }
  
  return(feature_matrix)
}


#' Extract combined features (RNA HVGs + ATAC peaks) for all target genes
#'
#' @param seurat_obj Seurat object
#' @param target_genes Vector of target gene names
#' @param gene_annotations GRanges object
#' @param n_hvg Number of HVGs for RNA features
#' @param window_size Window size for ATAC peaks (bp from TSS)
#' @param metacell_rna Metacell RNA matrix (optional)
#' @param metacell_atac Metacell ATAC matrix (optional)
#' @return List with feature matrices and target values for each gene
#' @export
extract_features <- function(seurat_obj, target_genes, gene_annotations = NULL,
                              n_hvg = 1000, window_size = 250000,
                              metacell_rna = NULL, metacell_atac = NULL) {
  
  cat("=== Extracting Combined Features (RNA + ATAC) ===\n\n")
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
    
    # Get ATAC features for this gene (peaks within window of TSS)
    use_mc_atac <- !is.null(metacell_atac)
    atac_features <- extract_atac_features_for_gene(
      seurat_obj = seurat_obj,
      target_gene = gene,
      gene_annotations = gene_annotations,
      window_size = window_size,
      use_metacells = use_mc_atac,
      metacell_matrix = metacell_atac
    )
    
    # Combine RNA + ATAC features
    if (!is.null(atac_features) && ncol(atac_features) > 0) {
      X <- cbind(rna_features, atac_features)
    } else {
      X <- rna_features
    }
    
    # Store results
    results[[gene]] <- list(
      X = X,
      y = y,
      n_rna_features = ncol(rna_features),
      n_atac_features = if (!is.null(atac_features)) ncol(atac_features) else 0,
      peak_names = if (!is.null(atac_features)) colnames(atac_features) else NULL
    )
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  cat(sprintf("\n✓ Extracted combined features for %d genes\n", length(results)))
  
  return(results)
}

#' Get gene annotations from Seurat object or external file
#'
#' @param seurat_obj Seurat object (optional)
#' @param annotation_file GTF file path (optional)
#' @param genome Genome assembly ("mm10" or "hg38")
#' @return GRanges object with gene annotations
#' @export
get_gene_annotations <- function(seurat_obj = NULL, annotation_file = NULL,
                                  genome = "mm10") {
  
  cat("=== Loading Gene Annotations ===\n\n")
  
  # Try to get from Seurat ATAC assay
  if (!is.null(seurat_obj) && "ATAC" %in% names(seurat_obj@assays)) {
    annotations <- Signac::Annotation(seurat_obj[["ATAC"]])
    if (!is.null(annotations)) {
      cat(sprintf("Loaded %d gene annotations from Seurat object\n", 
                  length(annotations)))
      return(annotations)
    }
  }
  
  # Load from GTF file
  if (!is.null(annotation_file) && file.exists(annotation_file)) {
    cat(sprintf("Loading annotations from: %s\n", annotation_file))
    annotations <- rtracklayer::import(annotation_file)
    
    # Filter to genes
    annotations <- annotations[annotations$type == "gene"]
    cat(sprintf("Loaded %d gene annotations from file\n", length(annotations)))
    return(annotations)
  }
  
  # Use built-in EnsDb if available
  if (genome == "mm10") {
    if (requireNamespace("EnsDb.Mmusculus.v79", quietly = TRUE)) {
      cat("Using EnsDb.Mmusculus.v79 annotations\n")
      annotations <- GenomicFeatures::genes(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
      return(annotations)
    }
  } else if (genome == "hg38") {
    if (requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
      cat("Using EnsDb.Hsapiens.v86 annotations\n")
      annotations <- GenomicFeatures::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
      return(annotations)
    }
  }
  
  stop("Could not load gene annotations. Please provide annotation_file parameter.")
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
