# ============================================================================
# Step_050: Metacell Creation via k-NN Smoothing using scVI Latent (scRNA-only)
# ============================================================================
# This script creates metacells using k-NN smoothing in scVI latent space.
# Uses pre-computed scVI latent representations from Python autoencoder.
#
# Input: 
#   - Seurat object with data splits (from multiome pipeline)
#   - scVI latent space CSV files (from Step_045)
# Output: Smoothed RNA expression matrices for train/val/test splits
#
# Key Features:
#   - Uses scVI latent space instead of PCA
#   - Z-score normalization fitted on TRAINING data only
#   - Smoothing done WITHIN each split (no data leakage)
#   - All cells preserved (1-to-1 mapping)
#   - Publication-ready diagnostic plots
#
# Usage:
#   1. Edit config.R with your parameters (DIM_REDUCTION_METHOD = "scvi")
#   2. Ensure scVI latent files exist in OUTPUT_LATENT_DIR
#   3. Run: Rscript Step_050.Metacell_Creation_scVI_Automated.R
# ============================================================================

# Get the directory of this script to source config.R
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  return(".")
}
script_dir <- get_script_dir()

# Source configuration file
config_path <- file.path(script_dir, "config.R")
if (!file.exists(config_path)) {
  config_path <- "config.R"
}

if (!file.exists(config_path)) {
  stop("config.R not found! Please ensure config.R is in the same directory as this script.")
}

cat("Loading configuration from:", config_path, "\n")
source(config_path)

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(RANN)         # Fast k-NN search
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(Cairo)
  library(viridis)
})

# ============================================================================
# PUBLICATION-READY PLOTTING SETUP
# ============================================================================

# Colorblind-friendly palette (Wong 2011)
colorblind_palette <- c(
  "#E69F00",  # Orange
  "#56B4E9",  # Sky blue
  "#009E73",  # Bluish green
  "#F0E442",  # Yellow
  "#0072B2",  # Blue
  "#D55E00",  # Vermillion
  "#CC79A7",  # Reddish purple
  "#999999"   # Gray
)

# Split colors for publication
split_colors_pub <- c(
  "train" = "#0072B2",       # Blue
  "validation" = "#E69F00",  # Orange
  "test" = "#009E73"         # Green
)

# Publication theme
theme_publication <- function(base_size = 12, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5, 
                                margin = margin(b = 10)),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "gray40",
                                   margin = margin(b = 10)),
      plot.caption = element_text(size = base_size - 2, hjust = 1, color = "gray50"),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
      strip.text = element_text(size = base_size, face = "bold", color = "black"),
      plot.margin = margin(15, 15, 15, 15)
    )
}

# Save publication plot function
save_publication_plot <- function(plot, filename, width = 10, height = 8, dpi = 600, 
                                   create_pdf = TRUE) {
  ggsave(
    filename = paste0(filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  cat(sprintf("  Saved: %s.png\n", basename(filename)))
  
  if (create_pdf) {
    ggsave(
      filename = paste0(filename, ".pdf"),
      plot = plot,
      width = width,
      height = height,
      device = cairo_pdf,
      bg = "white"
    )
    cat(sprintf("  Saved: %s.pdf\n", basename(filename)))
  }
}

# ============================================================================
# Z-SCORE FUNCTIONS
# ============================================================================

z_fit_cols <- function(M) { 
  mu <- colMeans(M)
  sd <- apply(M, 2, sd)
  sd[!is.finite(sd) | sd < 1e-8] <- 1
  list(mu = mu, sd = sd)
}

z_apply_cols <- function(M, mu, sd) { 
  sweep(sweep(M, 2, mu, "-"), 2, sd, "/")
}

# ============================================================================
# K-NN SMOOTHING FUNCTIONS
# ============================================================================

build_knn_smooth <- function(scores, cell_names, k_smooth = 20, seed = 1) {
  set.seed(seed)
  n <- nrow(scores)
  total_k <- k_smooth
  k_use <- max(1, min(total_k - 1, n - 1))
  
  if (n <= 1) {
    warning(sprintf("Too few cells (%d) for k-NN smoothing. Using identity.", n))
    return(sparseMatrix(
      i = 1, j = 1, x = 1,
      dims = c(n, n),
      dimnames = list(cell_names, cell_names)
    ))
  }
  
  nn <- RANN::nn2(
    data       = scores,
    query      = scores,
    k          = k_use + 1,
    searchtype = "standard"
  )
  
  idx <- nn$nn.idx
  total_k_eff <- k_use + 1
  weight  <- 1 / total_k_eff
  
  i_idx <- as.vector(t(idx))
  j_idx <- rep(seq_len(n), each = total_k)
  x_vals <- rep(weight, length(i_idx))
  
  A <- sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = x_vals,
    dims = c(n, n),
    dimnames = list(cell_names, cell_names)
  )
  
  stopifnot(
    identical(rownames(A), cell_names),
    identical(colnames(A), cell_names)
  )
  A
}

smooth_rna_expression <- function(rna_counts, A_smooth) {
  rna_smoothed <- rna_counts %*% A_smooth
  
  rna_lib <- Matrix::colSums(rna_smoothed)
  rna_cpm <- t(t(rna_smoothed) / (rna_lib + 1e-8)) * 1e6
  rna_log1p <- log1p(rna_cpm)
  
  list(
    rna_log1p = rna_log1p,
    rna_lib = rna_lib,
    rna_counts = rna_smoothed
  )
}

# ============================================================================
# PRINT CONFIGURATION
# ============================================================================
print_config()
print_output_directories()

# Create output directories
create_output_directories()

# Create plots directory
plots_dir <- file.path(OUTPUT_METACELLS_DIR, "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 050: Metacell Creation via k-NN Smoothing (scVI Latent, scRNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat(sprintf("Dimensionality Reduction Method: %s\n", DIM_REDUCTION_METHOD))
cat(sprintf("Output directory: %s\n\n", OUTPUT_METACELLS_DIR))

# ============================================================================
# LOAD INPUT DATA
# ============================================================================
cat("=== Loading Seurat object with data splits ===\n\n")

input_file <- path.expand(INPUT_SEURAT_SPLITS)
cat("Input file:", input_file, "\n")

if (!file.exists(input_file)) {
  stop("ERROR: Input file not found: ", input_file, "\n",
       "Please ensure the multiome preprocessing pipeline was run first.")
}

seurat_obj <- readRDS(input_file)

cat(sprintf("Loaded object: %d cells, %d genes (RNA)\n", 
            ncol(seurat_obj), nrow(seurat_obj[["RNA"]])))

# Check for data splits
if (!"data_split" %in% colnames(seurat_obj@meta.data)) {
  stop("ERROR: 'data_split' column not found in metadata. Run Step_040 first.")
}

split_counts <- table(seurat_obj$data_split)
cat("Data splits:\n")
for (split in names(split_counts)) {
  cat(sprintf("  %s: %d cells\n", split, split_counts[split]))
}

# ============================================================================
# LOAD SCVI LATENT SPACES
# ============================================================================
cat("\n=== Loading scVI latent spaces ===\n\n")

# Define latent file paths (per-split files)
latent_dir <- path.expand(OUTPUT_LATENT_DIR)

# Check for different possible naming conventions
possible_patterns <- list(
  # Pattern 1: latent_scvi_rna_train.csv, latent_scvi_rna_val.csv, latent_scvi_rna_test.csv
  list(train = "latent_scvi_rna_train.csv", 
       validation = "latent_scvi_rna_val.csv", 
       test = "latent_scvi_rna_test.csv"),
  # Pattern 2: scvi_rna_train.csv, etc.
  list(train = "scvi_rna_train.csv", 
       validation = "scvi_rna_val.csv", 
       test = "scvi_rna_test.csv"),
  # Pattern 3: Single file with all cells
  list(all = "latent_scvi_rna_all.csv")
)

latent_files <- list()
use_single_file <- FALSE

# Try to find the files
for (pattern in possible_patterns) {
  if ("all" %in% names(pattern)) {
    all_file <- file.path(latent_dir, pattern$all)
    if (file.exists(all_file)) {
      cat("Found single latent file:", all_file, "\n")
      use_single_file <- TRUE
      latent_files$all <- all_file
      break
    }
  } else {
    train_file <- file.path(latent_dir, pattern$train)
    if (file.exists(train_file)) {
      latent_files$train <- train_file
      latent_files$validation <- file.path(latent_dir, pattern$validation)
      latent_files$test <- file.path(latent_dir, pattern$test)
      cat("Found per-split latent files:\n")
      cat("  Train:", latent_files$train, "\n")
      cat("  Val:", latent_files$validation, "\n")
      cat("  Test:", latent_files$test, "\n")
      break
    }
  }
}

if (length(latent_files) == 0) {
  stop("ERROR: Could not find scVI latent files in: ", latent_dir, "\n",
       "Expected one of:\n",
       "  - latent_scvi_rna_train.csv, latent_scvi_rna_val.csv, latent_scvi_rna_test.csv\n",
       "  - latent_scvi_rna_all.csv\n",
       "Please run Step_045 to generate scVI embeddings first.")
}

# Load latent spaces
if (use_single_file) {
  cat("\n[INFO] Loading single latent file...\n")
  latent_all <- read.csv(latent_files$all, row.names = 1, check.names = FALSE)
  
  # Identify latent columns
  latent_cols <- grep("^scVI|^latent|^z_", colnames(latent_all), value = TRUE)
  if (length(latent_cols) == 0) {
    # Use all numeric columns except data_split
    latent_cols <- setdiff(colnames(latent_all)[sapply(latent_all, is.numeric)], "data_split")
  }
  
  cat(sprintf("  Found %d latent dimensions: %s\n", 
              length(latent_cols), 
              paste(head(latent_cols, 3), collapse = ", ")))
}

# Set random seed
set.seed(SEED_METACELL)

# ============================================================================
# PROCESS TRAINING DATA
# ============================================================================
cat("\n=== Processing TRAINING data ===\n\n")

train_cells <- colnames(seurat_obj)[seurat_obj$data_split == "train"]
seurat_train <- seurat_obj[, train_cells]

cat(sprintf("Training cells: %d\n", length(train_cells)))

# Ensure RNA counts are sparse
DefaultAssay(seurat_train) <- "RNA"
rc <- GetAssayData(seurat_train, assay = "RNA", layer = "counts")
if (!inherits(rc, "dgCMatrix")) {
  seurat_train <- SetAssayData(seurat_train, assay = "RNA", layer = "counts", 
                               new.data = as(rc, "dgCMatrix"))
}

cat(sprintf("Genes: %d\n", nrow(seurat_train[["RNA"]])))

# Load training scVI latent
cat("\n[INFO] Loading scVI latent for training...\n")
if (use_single_file) {
  latent_train <- latent_all[train_cells, latent_cols, drop = FALSE]
} else {
  latent_train_df <- read.csv(latent_files$train, row.names = 1, check.names = FALSE)
  latent_cols <- grep("^scVI|^latent|^z_", colnames(latent_train_df), value = TRUE)
  if (length(latent_cols) == 0) {
    latent_cols <- setdiff(colnames(latent_train_df)[sapply(latent_train_df, is.numeric)], "data_split")
  }
  latent_train <- latent_train_df[, latent_cols, drop = FALSE]
}

rna_scores_train <- as.matrix(latent_train)
rna_scores_train <- rna_scores_train[train_cells, , drop = FALSE]

cat(sprintf("  scVI latent: %d cells × %d dims\n",
            nrow(rna_scores_train), ncol(rna_scores_train)))

# Z-score standardization (fit on training data)
cat("[INFO] Computing z-score parameters from training...\n")
rna_zs <- z_fit_cols(rna_scores_train)
rna_train_z <- z_apply_cols(rna_scores_train, rna_zs$mu, rna_zs$sd)

# Build k-NN smoothing matrix for training
cat(sprintf("\n[INFO] Building k-NN smoothing matrix (k=%d)...\n", K_NEIGHBORS))
A_smooth_train <- build_knn_smooth(rna_train_z, train_cells, 
                                    k_smooth = K_NEIGHBORS, seed = SEED_METACELL)

neighbors_per_cell <- Matrix::rowSums(A_smooth_train > 0)
cat(sprintf("  k-NN stats: mean=%.1f, min=%d, max=%d neighbors\n",
            mean(neighbors_per_cell), min(neighbors_per_cell), max(neighbors_per_cell)))

# Smooth training expression
cat("[INFO] Smoothing training expression...\n")
rna_counts_train <- GetAssayData(seurat_train, assay = "RNA", layer = "counts")
smoothed_train <- smooth_rna_expression(rna_counts_train, A_smooth_train)
cat(sprintf("[OK] Training: %d cells smoothed\n", ncol(smoothed_train$rna_log1p)))

# Save training transformation parameters
train_transforms <- list(
  rna_zs = rna_zs,
  k_neighbors = K_NEIGHBORS,
  latent_dims = ncol(rna_scores_train),
  latent_cols = latent_cols
)
saveRDS(train_transforms, file.path(OUTPUT_METACELLS_DIR, "train_transforms.rds"))
cat("[OK] Saved train_transforms.rds\n")

# Save smoothed training data
train_output <- list(
  rna_log1p = smoothed_train$rna_log1p,
  rna_lib = smoothed_train$rna_lib,
  rna_counts = smoothed_train$rna_counts,
  scvi_latent = rna_scores_train,
  scvi_latent_z = rna_train_z,
  cells = train_cells,
  split_name = "train",
  processing_seed = SEED_METACELL,
  k_smooth = K_NEIGHBORS
)
saveRDS(train_output, file.path(OUTPUT_METACELLS_DIR, "smoothed_train.rds"))
cat("[OK] Saved smoothed_train.rds\n")

# Store for combined plots
all_latent_scores <- list(train = rna_scores_train)
all_latent_z <- list(train = rna_train_z)
all_splits_info <- list(train = list(n_cells = length(train_cells)))

# ============================================================================
# PROCESS VALIDATION AND TEST DATA
# ============================================================================
cat("\n=== Processing VALIDATION and TEST data ===\n\n")

for (split_name in c("validation", "test")) {
  cat(sprintf("\n--- Processing %s ---\n", toupper(split_name)))
  
  split_cells <- colnames(seurat_obj)[seurat_obj$data_split == split_name]
  seurat_split <- seurat_obj[, split_cells]
  
  cat(sprintf("%s cells: %d\n", split_name, length(split_cells)))
  
  # Ensure sparse
  DefaultAssay(seurat_split) <- "RNA"
  rc <- GetAssayData(seurat_split, assay = "RNA", layer = "counts")
  if (!inherits(rc, "dgCMatrix")) {
    seurat_split <- SetAssayData(seurat_split, assay = "RNA", layer = "counts",
                                 new.data = as(rc, "dgCMatrix"))
  }
  
  # Load scVI latent for this split
  cat(sprintf("[INFO] Loading scVI latent for %s...\n", split_name))
  if (use_single_file) {
    latent_split <- latent_all[split_cells, latent_cols, drop = FALSE]
  } else {
    latent_split_df <- read.csv(latent_files[[split_name]], row.names = 1, check.names = FALSE)
    latent_split <- latent_split_df[, latent_cols, drop = FALSE]
  }
  
  rna_scores_split <- as.matrix(latent_split)
  rna_scores_split <- rna_scores_split[split_cells, , drop = FALSE]
  
  cat(sprintf("  scVI latent: %d cells × %d dims\n",
              nrow(rna_scores_split), ncol(rna_scores_split)))
  
  # Apply TRAINING z-score parameters
  cat("[INFO] Applying training z-score normalization...\n")
  rna_split_z <- z_apply_cols(rna_scores_split, train_transforms$rna_zs$mu, train_transforms$rna_zs$sd)
  
  # Build k-NN smoothing matrix WITHIN this split only
  cat(sprintf("[INFO] Building k-NN within %s cells (k=%d)...\n", split_name, K_NEIGHBORS))
  A_smooth_split <- build_knn_smooth(rna_split_z, split_cells, 
                                      k_smooth = K_NEIGHBORS, seed = SEED_METACELL)
  
  neighbors_per_cell <- Matrix::rowSums(A_smooth_split > 0)
  cat(sprintf("  k-NN stats: mean=%.1f, min=%d, max=%d neighbors\n",
              mean(neighbors_per_cell), min(neighbors_per_cell), max(neighbors_per_cell)))
  
  # Smooth expression
  cat(sprintf("[INFO] Smoothing %s expression...\n", split_name))
  rna_counts_split <- GetAssayData(seurat_split, assay = "RNA", layer = "counts")
  smoothed_split <- smooth_rna_expression(rna_counts_split, A_smooth_split)
  cat(sprintf("[OK] %s: %d cells smoothed\n", split_name, ncol(smoothed_split$rna_log1p)))
  
  # Save smoothed data
  split_output <- list(
    rna_log1p = smoothed_split$rna_log1p,
    rna_lib = smoothed_split$rna_lib,
    rna_counts = smoothed_split$rna_counts,
    scvi_latent = rna_scores_split,
    scvi_latent_z = rna_split_z,
    cells = split_cells,
    split_name = split_name,
    processing_seed = SEED_METACELL,
    k_smooth = K_NEIGHBORS
  )
  saveRDS(split_output, file.path(OUTPUT_METACELLS_DIR, paste0("smoothed_", split_name, ".rds")))
  cat(sprintf("[OK] Saved smoothed_%s.rds\n", split_name))
  
  # Store for plots
  all_latent_scores[[split_name]] <- rna_scores_split
  all_latent_z[[split_name]] <- rna_split_z
  all_splits_info[[split_name]] <- list(n_cells = length(split_cells))
}

# ============================================================================
# GENERATE PUBLICATION-READY PLOTS
# ============================================================================
cat("\n=== Generating publication-ready plots ===\n\n")

# Determine number of latent dimensions
n_latent_dims <- ncol(rna_scores_train)

# ----------------------------------------------------------------------------
# 1. LATENT SPACE EMBEDDING BY SPLIT (Dim 1 vs Dim 2)
# ----------------------------------------------------------------------------
latent_plot_df <- do.call(rbind, lapply(names(all_latent_scores), function(split) {
  scores <- all_latent_scores[[split]]
  data.frame(
    Dim1 = scores[, 1],
    Dim2 = scores[, 2],
    Dim3 = if(ncol(scores) >= 3) scores[, 3] else NA,
    split = split,
    cell = rownames(scores)
  )
}))
latent_plot_df$split <- factor(latent_plot_df$split, levels = c("train", "validation", "test"))

# Dim1 vs Dim2
p_latent12 <- ggplot(latent_plot_df, aes(x = Dim1, y = Dim2, color = split)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = split_colors_pub, name = "Data Split",
                     labels = c("Train", "Validation", "Test")) +
  labs(
    x = "scVI Dimension 1",
    y = "scVI Dimension 2",
    title = "scVI Latent Space Embedding by Data Split",
    subtitle = sprintf("%s | %d dimensions | k=%d neighbors", SAMPLE_NAME, n_latent_dims, K_NEIGHBORS)
  ) +
  theme_publication() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

save_publication_plot(p_latent12, 
                       file.path(plots_dir, paste0(SAMPLE_NAME, "_scvi_dim1_dim2_by_split")),
                       width = 10, height = 8)

# Dim1 vs Dim3
if (!all(is.na(latent_plot_df$Dim3))) {
  p_latent13 <- ggplot(latent_plot_df, aes(x = Dim1, y = Dim3, color = split)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_color_manual(values = split_colors_pub, name = "Data Split",
                       labels = c("Train", "Validation", "Test")) +
    labs(
      x = "scVI Dimension 1",
      y = "scVI Dimension 3",
      title = "scVI Latent Space (Dim1 vs Dim3)",
      subtitle = sprintf("%s", SAMPLE_NAME)
    ) +
    theme_publication() +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  save_publication_plot(p_latent13, 
                         file.path(plots_dir, paste0(SAMPLE_NAME, "_scvi_dim1_dim3_by_split")),
                         width = 10, height = 8)
}

# ----------------------------------------------------------------------------
# 2. SPLIT-FACETED LATENT PLOT
# ----------------------------------------------------------------------------
p_latent_facet <- ggplot(latent_plot_df, aes(x = Dim1, y = Dim2, color = split)) +
  geom_point(alpha = 0.5, size = 0.6) +
  scale_color_manual(values = split_colors_pub, name = "Data Split",
                     labels = c("Train", "Validation", "Test")) +
  facet_wrap(~ split, labeller = labeller(split = c(
    "train" = sprintf("Train (n=%d)", all_splits_info$train$n_cells),
    "validation" = sprintf("Validation (n=%d)", all_splits_info$validation$n_cells),
    "test" = sprintf("Test (n=%d)", all_splits_info$test$n_cells)
  ))) +
  labs(
    x = "scVI Dimension 1",
    y = "scVI Dimension 2",
    title = "scVI Latent Space by Data Split (Faceted)",
    subtitle = sprintf("%s | k=%d neighbors", SAMPLE_NAME, K_NEIGHBORS)
  ) +
  theme_publication() +
  theme(legend.position = "none")

save_publication_plot(p_latent_facet, 
                       file.path(plots_dir, paste0(SAMPLE_NAME, "_scvi_faceted_by_split")),
                       width = 14, height = 5)

# ----------------------------------------------------------------------------
# 3. DENSITY PLOTS FOR LATENT DIMENSIONS BY SPLIT
# ----------------------------------------------------------------------------
p_density_d1 <- ggplot(latent_plot_df, aes(x = Dim1, fill = split, color = split)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  scale_fill_manual(values = split_colors_pub, name = "Split") +
  scale_color_manual(values = split_colors_pub, name = "Split") +
  labs(
    x = "scVI Dimension 1",
    y = "Density",
    title = "Dim 1 Distribution by Split"
  ) +
  theme_publication()

p_density_d2 <- ggplot(latent_plot_df, aes(x = Dim2, fill = split, color = split)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  scale_fill_manual(values = split_colors_pub, name = "Split") +
  scale_color_manual(values = split_colors_pub, name = "Split") +
  labs(
    x = "scVI Dimension 2",
    y = "Density",
    title = "Dim 2 Distribution by Split"
  ) +
  theme_publication()

p_density_combined <- p_density_d1 + p_density_d2 +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "scVI Latent Dimension Distributions by Data Split",
    subtitle = sprintf("%s | Validation of split randomization", SAMPLE_NAME),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40")
    )
  )

save_publication_plot(p_density_combined, 
                       file.path(plots_dir, paste0(SAMPLE_NAME, "_scvi_density_by_split")),
                       width = 14, height = 6)

# ----------------------------------------------------------------------------
# 4. LATENT DIMENSION VARIANCE DISTRIBUTION
# ----------------------------------------------------------------------------
latent_var <- apply(rna_scores_train, 2, var)
var_df <- data.frame(
  Dimension = 1:length(latent_var),
  Variance = latent_var
)

p_var <- ggplot(var_df, aes(x = Dimension, y = Variance)) +
  geom_bar(stat = "identity", fill = colorblind_palette[3], color = "black", alpha = 0.7) +
  labs(
    x = "scVI Latent Dimension",
    y = "Variance",
    title = "scVI Latent Dimension Variance",
    subtitle = sprintf("%s | %d dimensions (training data)", SAMPLE_NAME, n_latent_dims)
  ) +
  theme_publication()

save_publication_plot(p_var, 
                       file.path(plots_dir, paste0(SAMPLE_NAME, "_scvi_variance_per_dim")),
                       width = 10, height = 6)

# ----------------------------------------------------------------------------
# 5. SMOOTHING EFFECT VISUALIZATION
# ----------------------------------------------------------------------------
raw_train <- as.matrix(GetAssayData(seurat_train, assay = "RNA", layer = "counts")[1:min(100, nrow(seurat_train[["RNA"]])), ])
smoothed_train_mat <- as.matrix(smoothed_train$rna_counts[1:min(100, nrow(smoothed_train$rna_counts)), ])

var_raw <- apply(raw_train, 1, var)
var_smooth <- apply(smoothed_train_mat, 1, var)

var_df <- data.frame(
  gene = names(var_raw),
  raw = var_raw,
  smoothed = var_smooth
) %>%
  mutate(var_ratio = smoothed / (raw + 1e-10))

p_var_reduction <- ggplot(var_df, aes(x = log10(raw + 1), y = log10(smoothed + 1))) +
  geom_point(alpha = 0.5, color = colorblind_palette[3], size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    x = "Log10 Variance (Raw Counts + 1)",
    y = "Log10 Variance (Smoothed Counts + 1)",
    title = "k-NN Smoothing Effect on Gene Variance (scVI)",
    subtitle = sprintf("%s | k=%d neighbors | First 100 genes", SAMPLE_NAME, K_NEIGHBORS),
    caption = sprintf("Points below diagonal = variance reduction\nMedian reduction: %.1f%%",
                      (1 - median(var_df$var_ratio, na.rm = TRUE)) * 100)
  ) +
  theme_publication() +
  coord_fixed()

save_publication_plot(p_var_reduction, 
                       file.path(plots_dir, paste0(SAMPLE_NAME, "_scvi_smoothing_variance_reduction")),
                       width = 8, height = 8)

# ----------------------------------------------------------------------------
# 6. SUMMARY STATISTICS BAR PLOT
# ----------------------------------------------------------------------------
summary_df <- data.frame(
  Split = c("Train", "Validation", "Test"),
  Cells = c(all_splits_info$train$n_cells, 
            all_splits_info$validation$n_cells, 
            all_splits_info$test$n_cells),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Percentage = sprintf("%.1f%%", Cells / sum(Cells) * 100),
    k_Neighbors = K_NEIGHBORS
  )

p_summary <- ggplot(summary_df, aes(x = Split, y = Cells, fill = Split)) +
  geom_col(color = "black", alpha = 0.8, width = 0.6) +
  geom_text(aes(label = format(Cells, big.mark = ",")), 
            vjust = -0.5, size = 4, fontface = "bold") +
  geom_text(aes(label = Percentage), vjust = 1.5, size = 3.5, color = "white") +
  scale_fill_manual(values = split_colors_pub) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "",
    y = "Number of Cells",
    title = "Data Split Distribution",
    subtitle = sprintf("%s | scVI-based Metacell Creation", SAMPLE_NAME),
    caption = sprintf("Method: scVI latent + k-NN smoothing (k=%d, %d dims)", K_NEIGHBORS, n_latent_dims)
  ) +
  theme_publication() +
  theme(legend.position = "none")

save_publication_plot(p_summary, 
                       file.path(plots_dir, paste0(SAMPLE_NAME, "_scvi_split_summary")),
                       width = 8, height = 6)

# ============================================================================
# COMPLETION SUMMARY
# ============================================================================
cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 050 COMPLETE (scVI Latent, scRNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("k-NN smoothing completed for all data splits (scVI Latent, RNA-ONLY)\n\n")
cat("Parameters:\n")
cat(sprintf("  - k = %d neighbors\n", K_NEIGHBORS))
cat(sprintf("  - scVI latent dims = %d\n", n_latent_dims))
cat(sprintf("  - Random seed = %d\n", SEED_METACELL))
cat(sprintf("  - Method: %s\n", DIMRED_METHOD_SUFFIX))

cat("\nKey features:\n")
cat("  ✓ Used scVI latent space for k-NN geometry\n")
cat("  ✓ Z-score fitted on training, applied to val/test\n")
cat("  ✓ Val/Test: Smoothed only within each split (no leakage)\n")
cat("  ✓ All cells preserved (1-to-1 mapping)\n")
cat("  ✓ RNA-ONLY: No ATAC data used\n")

cat("\nFiles saved:\n")
cat(sprintf("  - %s/train_transforms.rds\n", basename(OUTPUT_METACELLS_DIR)))
cat(sprintf("  - %s/smoothed_train.rds\n", basename(OUTPUT_METACELLS_DIR)))
cat(sprintf("  - %s/smoothed_validation.rds\n", basename(OUTPUT_METACELLS_DIR)))
cat(sprintf("  - %s/smoothed_test.rds\n", basename(OUTPUT_METACELLS_DIR)))

cat("\nPublication-ready plots saved:\n")
cat(sprintf("  - %s_scvi_dim1_dim2_by_split.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_scvi_dim1_dim3_by_split.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_scvi_faceted_by_split.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_scvi_density_by_split.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_scvi_variance_per_dim.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_scvi_smoothing_variance_reduction.pdf/png\n", SAMPLE_NAME))
cat(sprintf("  - %s_scvi_split_summary.pdf/png\n", SAMPLE_NAME))

cat("\nNext: Run Step_060.Feature_Extraction_Automated.R\n")
