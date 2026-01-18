#!/usr/bin/env Rscript
# ============================================================================
# Step_081: SHAP Feature Importance for Deep Neural Networks (scRNA-only)
# ============================================================================
# This script computes SHAP values for a trained DeepNN model to explain
# feature importance in gene expression prediction.
#
# For each gene:
#   - Loads gene-specific features (train/val/test)
#   - Re-applies min-max scaling (same as training)
#   - Loads saved Keras model
#   - Computes SHAP values on TEST set using KernelExplainer
#   - Saves feature importance CSV and summary plot
#
# Note: This is the RNA-only version (HVG expression features only)
#
# Usage:
#   Rscript Step_081.SHAP_Feature_Importance_Automated.R <GENE_NAME> [GENE_SET]
#   
# Examples:
#   Rscript Step_081.SHAP_Feature_Importance_Automated.R Larp1b HVG
#   Rscript Step_081.SHAP_Feature_Importance_Automated.R Sox2
#
# For batch processing, use with SLURM array jobs
# ============================================================================

# Record start time
start_time <- Sys.time()

# ============================================================================
# PARSE COMMAND LINE ARGUMENTS
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript Step_081.SHAP_Feature_Importance_Automated.R <GENE_NAME> [GENE_SET]\n",
       "Example: Rscript Step_081.SHAP_Feature_Importance_Automated.R Larp1b HVG\n")
}

gene_name <- trimws(args[1])
gene_set_arg <- if (length(args) >= 2) trimws(args[2]) else NULL

# ============================================================================
# SOURCE CONFIGURATION
# ============================================================================
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  return(".")
}
script_dir <- get_script_dir()

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
# LOAD LIBRARIES
# ============================================================================
library(reticulate)
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(Cairo)
})

# ============================================================================
# PUBLICATION-READY PLOTTING SETUP
# ============================================================================

# Colorblind-friendly palette (Wong 2011)
colorblind_palette <- c(

"#0072B2",  # Blue
  "#56B4E9",  # Sky blue
  "#009E73",  # Bluish green
  "#E69F00",  # Orange
  "#D55E00",  # Vermillion
  "#CC79A7",  # Reddish purple
  "#F0E442",  # Yellow
  "#999999"   # Gray
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

# Save publication-quality plots (PNG + PDF)
save_publication_plot <- function(plot, filename, width = 10, height = 8, dpi = 600, 
                                   create_pdf = TRUE) {
  # Save PNG
  ggsave(
    filename = paste0(filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  cat(sprintf("  Saved: %s.png\n", basename(filename)))
  
  # Save PDF (vector format)
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
# PYTHON ENVIRONMENT SETUP
# ============================================================================
use_condaenv(CONDA_ENV_NAME, required = TRUE)
cat("Python configuration:\n")
print(py_config())

# CPU-friendly thread settings for cluster computing
Sys.setenv(TF_NUM_INTRAOP_THREADS = "1", TF_NUM_INTEROP_THREADS = "1")

# Import Python modules
tf    <- import("tensorflow")
keras <- tf$keras
np    <- import("numpy", convert = FALSE)
shap  <- import("shap", convert = FALSE)
plt   <- import("matplotlib.pyplot", convert = FALSE)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Determine gene set from argument or config
if (!is.null(gene_set_arg)) {
  gene_set_name <- gene_set_arg
} else if (MODEL_GENE_SET == "HVG") {
  gene_set_name <- "HVG"
} else if (MODEL_GENE_SET %in% c("Random", "Custom")) {
  gene_set_name <- "Custom"
} else {
  gene_set_name <- "HVG"  # Default
}

# SHAP-specific parameters
SHAP_N_BACKGROUND <- 20          # Number of kmeans clusters for background
SHAP_TOP_N_FEATURES <- 30        # Number of features to show in bar plot

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 081: SHAP Feature Importance for DeepNN (scRNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat("Configuration:\n")
cat(sprintf("  Gene: %s\n", gene_name))
cat(sprintf("  Gene set: %s\n", gene_set_name))
cat(sprintf("  Sample: %s\n", SAMPLE_NAME))
cat(sprintf("  Feature type: HVG expression only (RNA-only)\n"))

# ============================================================================
# SET PATHS
# ============================================================================

# Features file
features_rds <- file.path(
  path.expand(OUTPUT_FEATURES_DIR),
  gene_set_name,
  "gene_specific_features.rds"
)

# Model directory and file
model_dir <- file.path(
  path.expand(OUTPUT_MODELS_NN_DIR),
  gene_set_name,
  gene_name
)
model_path <- file.path(model_dir, "best_model.keras")

# Output directories
shap_out_dir <- file.path(model_dir, "shap_feature_importance")
dir.create(shap_out_dir, showWarnings = FALSE, recursive = TRUE)

# Figure output directory
fig_dir <- file.path(
  path.expand(OUTPUT_FIGURES_DIR),
  "NEURAL_NETWORKS",
  gene_set_name,
  "SHAP",
  gene_name
)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Output file paths
imp_csv <- file.path(shap_out_dir, "shap_feature_importance.csv")
shap_plot_path <- file.path(fig_dir, paste0("shap_summary_", gene_name, ".png"))

cat(sprintf("  Features file: %s\n", features_rds))
cat(sprintf("  Model path: %s\n", model_path))
cat(sprintf("  SHAP output dir: %s\n", shap_out_dir))
cat(sprintf("  Figure output dir: %s\n", fig_dir))
cat("\n")

# ============================================================================
# HELPER FUNCTION: MIN-MAX SCALING
# ============================================================================

scale_mm <- function(M, mn, rg, zv = integer(0)) {
  M <- as.matrix(M)
  M <- sweep(M, 2, mn, "-")
  M <- sweep(M, 2, rg, "/")
  if (length(zv)) M[, zv] <- 0
  M
}

# ============================================================================
# VALIDATE INPUT FILES
# ============================================================================

if (!file.exists(features_rds)) {
  stop("Feature file not found: ", features_rds)
}

if (!file.exists(model_path)) {
  stop("Model file not found: ", model_path, 
       "\nMake sure Step_071 has been run for this gene.")
}

# ============================================================================
# LOAD FEATURES
# ============================================================================

cat("Loading features...\n")
gene_features <- readRDS(features_rds)

if (!gene_name %in% names(gene_features)) {
  stop("Gene '", gene_name, "' not found in gene_specific_features.rds\n",
       "Available genes: ", paste(head(names(gene_features), 10), collapse = ", "), "...")
}

gene_data <- gene_features[[gene_name]]

X_tr <- gene_data$train$X
y_tr <- gene_data$train$y
cells_tr <- gene_data$train$cells

X_va <- gene_data$validation$X
y_va <- gene_data$validation$y
cells_va <- gene_data$validation$cells

X_te <- gene_data$test$X
y_te <- gene_data$test$y
cells_te <- gene_data$test$cells

cat(sprintf("Features: %d | Train: %d | Val: %d | Test: %d\n",
            ncol(X_tr), nrow(X_tr), nrow(X_va), nrow(X_te)))

# ============================================================================
# REMOVE TARGET GENE FROM PREDICTORS (if present)
# ============================================================================

hvg_feature_names <- colnames(X_tr)
target_in_features <- grepl(paste0("^", gene_name, "$"), hvg_feature_names)

if (any(target_in_features)) {
  cat("Target gene found in predictors – removing it\n")
  X_tr <- X_tr[, !target_in_features, drop = FALSE]
  X_va <- X_va[, !target_in_features, drop = FALSE]
  X_te <- X_te[, !target_in_features, drop = FALSE]
}

# ============================================================================
# APPLY MIN-MAX SCALING (same as training)
# ============================================================================

cat("Re-scaling features with train min/max...\n")
mins   <- apply(as.matrix(X_tr), 2, min, na.rm = TRUE)
maxs   <- apply(as.matrix(X_tr), 2, max, na.rm = TRUE)
ranges <- maxs - mins
zerov  <- which(!is.finite(ranges) | ranges == 0)

if (length(zerov)) {
  ranges[zerov] <- 1
  cat("Zero-variance features:", length(zerov), "\n")
}

X_tr_s <- scale_mm(X_tr, mins, ranges, zerov)
X_va_s <- scale_mm(X_va, mins, ranges, zerov)
X_te_s <- scale_mm(X_te, mins, ranges, zerov)

# Convert to plain matrices
X_tr_s <- as.matrix(X_tr_s)
X_va_s <- as.matrix(X_va_s)
X_te_s <- as.matrix(X_te_s)

# ============================================================================
# LOAD KERAS MODEL
# ============================================================================

model_path_abs <- normalizePath(model_path, mustWork = TRUE)
cat("Loading model from:", model_path_abs, "\n")

model <- keras$models$load_model(model_path_abs)

# ============================================================================
# BUILD NUMPY ARRAYS FOR SHAP
# ============================================================================

cat("Preparing data for SHAP...\n")

# Test set for explanation
X_test <- X_te_s
feature_names <- colnames(X_te)

# Convert to NumPy arrays
X_test_np <- np$array(X_test)

# Background sample for KernelExplainer – small subset of training data using kmeans
bg_kmeans <- shap$kmeans(np$array(X_tr_s), as.integer(SHAP_N_BACKGROUND))
X_bg_np <- bg_kmeans$data

cat(sprintf("Test set shape: %s\n", paste(py_to_r(X_test_np$shape), collapse = " x ")))
cat(sprintf("Background samples: %d\n", SHAP_N_BACKGROUND))

# ============================================================================
# DEFINE PREDICTION FUNCTION FOR SHAP
# ============================================================================

# Function that takes a NumPy array and returns predictions
f <- function(x) {
  preds <- model$predict(x)
  preds
}

# ============================================================================
# COMPUTE SHAP VALUES
# ============================================================================

cat("Fitting SHAP KernelExplainer...\n")
explainer <- shap$KernelExplainer(f, X_bg_np)

cat("Computing SHAP values on test set (this can take a while)...\n")
shap_values <- explainer$shap_values(X_test_np)

# For single-output model, shap_values may be a list of length 1
if (inherits(shap_values, "python.builtin.list")) {
  shap_values <- shap_values[[1]]
}

# Handle 3D output (samples x features x 1)
shape <- py_to_r(shap_values$shape)
if (length(shape) == 3L && shape[[3]] == 1L) {
  shap_values_2d <- shap_values$reshape(as.integer(c(shape[[1]], shape[[2]])))
} else {
  shap_values_2d <- shap_values
}

cat("SHAP values shape:", paste(py_to_r(shap_values_2d$shape), collapse = " x "), "\n")
cat("X_test_np shape:", paste(py_to_r(X_test_np$shape), collapse = " x "), "\n")

# ============================================================================
# COMPUTE FEATURE IMPORTANCE (Mean |SHAP|)
# ============================================================================

cat("Computing mean(|SHAP|) per feature...\n")

# Convert shap_values into R matrix
shap_vals_r <- py_to_r(shap_values_2d)

# Calculate mean absolute SHAP value per feature
mean_abs_shap <- colMeans(abs(shap_vals_r), na.rm = TRUE)

# Create feature importance data frame
imp_df <- data.frame(
  Gene = gene_name,
  Feature = feature_names,
  MeanAbsSHAP = mean_abs_shap,
  stringsAsFactors = FALSE
) %>%
  arrange(desc(MeanAbsSHAP))

# Save CSV
write_csv(imp_df, imp_csv)

cat("\nTop 10 features by mean |SHAP|:\n")
print(head(imp_df %>% select(Feature, MeanAbsSHAP), 10))

cat("\nFeature importance CSV saved to:", imp_csv, "\n")

# ============================================================================
# PLOT SHAP FEATURE IMPORTANCE (PUBLICATION-READY)
# ============================================================================

cat("\nPlotting SHAP feature importance (publication-ready)...\n")

# Show top N features
top_n <- min(SHAP_TOP_N_FEATURES, nrow(imp_df))
imp_df_plot <- imp_df[1:top_n, ]

# Truncate long feature names for display
imp_df_plot$Feature_Short <- sapply(imp_df_plot$Feature, function(x) {
  if (nchar(x) > 35) paste0(substr(x, 1, 32), "...") else x
})

# Create publication-ready ggplot
p_shap <- ggplot(imp_df_plot, aes(x = MeanAbsSHAP, y = reorder(Feature_Short, MeanAbsSHAP))) +
  geom_col(fill = "#0072B2", alpha = 0.85, color = "black", linewidth = 0.3) +
  labs(
    title = sprintf("SHAP Feature Importance: %s", gene_name),
    subtitle = sprintf("%s | Top %d features | RNA-only (HVG)", SAMPLE_NAME, top_n),
    x = "Mean |SHAP Value|",
    y = NULL,
    caption = "Higher values indicate greater contribution to prediction"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_publication(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 10))
  )

# Save using publication function (PNG + PDF)
shap_plot_base <- file.path(fig_dir, paste0("shap_summary_", gene_name))
save_publication_plot(p_shap, shap_plot_base, width = 12, height = 10)

# Also save in the model directory
shap_plot_model_base <- file.path(shap_out_dir, paste0("shap_feature_importance_", gene_name))
save_publication_plot(p_shap, shap_plot_model_base, width = 12, height = 10)

cat("SHAP summary plots saved\n")

# ============================================================================
# SAVE RAW SHAP VALUES (for further analysis)
# ============================================================================

shap_raw_file <- file.path(shap_out_dir, "shap_values_raw.rds")
saveRDS(list(
  gene = gene_name,
  shap_values = shap_vals_r,
  feature_names = feature_names,
  cell_ids = cells_te,
  mean_abs_shap = mean_abs_shap
), shap_raw_file)

cat("Raw SHAP values saved to:", shap_raw_file, "\n")

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 081 COMPLETE for gene:", gene_name, "(scRNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat(sprintf("Runtime: %.2f minutes\n", as.numeric(runtime)))
cat("\nOutputs (publication-ready):\n")
cat(sprintf("  - Feature importance CSV: %s\n", imp_csv))
cat(sprintf("  - SHAP summary plot (PNG 600dpi): %s.png\n", shap_plot_base))
cat(sprintf("  - SHAP summary plot (PDF vector): %s.pdf\n", shap_plot_base))
cat(sprintf("  - Raw SHAP values: %s\n", shap_raw_file))

cat("\nNext step: Run Step_082 to aggregate SHAP results across all genes\n")

cat("\n", paste(rep("=", 70), collapse = ""), "\n\n")
