#!/usr/bin/env Rscript
# ============================================================================
# Step_080: Plotting Model Results (scRNA-only, Automated)
# ============================================================================
# This script creates comprehensive visualizations comparing:
#   - Linear models (OLS, Ridge, Lasso, ElasticNet)
#   - Random Forest
#   - Deep Neural Network (DeepNN)
#
# All results are combined for easy comparison across model types.
# This is the RNA-only version (no ATAC/peak features).
#
# Input: Model outputs from Step_070 and Step_071
# Output: Publication-ready PNG/PDF figures and summary CSVs
#
# Usage:
#   1. Ensure config.R is properly configured
#   2. Run: Rscript Step_080.Plotting_Results_Automated.R
# ============================================================================

# Record start time
start_time <- Sys.time()

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
# LOAD LIBRARIES
# ============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(viridis)
  library(ggrepel)
  library(scales)
  library(forcats)
  library(tidyr)
  library(dplyr)
  library(Cairo)
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
# PLOTTING SETTINGS (Colorblind-Friendly)
# ============================================================================

# Model color palette (Wong 2011 inspired)
model_palette <- c(
  "OLS" = "#0072B2",           # Blue
  "Ridge" = "#56B4E9",         # Sky blue
  "Lasso" = "#009E73",         # Bluish green
  "ElasticNet" = "#E69F00",    # Orange
  "RandomForest" = "#D55E00",  # Vermillion
  "DeepNN" = "#CC79A7"         # Reddish purple
)

# Model type grouping
model_type_palette <- c(
  "Linear" = "#0072B2",        # Blue
  "Tree-Based" = "#D55E00",    # Vermillion
  "Deep Learning" = "#CC79A7"  # Reddish purple
)

# Split colors
split_palette <- c(
  "Train" = "#0072B2",         # Blue
  "Val" = "#E69F00",           # Orange
  "Test" = "#D55E00"           # Vermillion
)

# Set default theme to publication
theme_set(theme_publication(base_size = 12))

# Model order for consistent plotting
model_order <- c("OLS", "Ridge", "Lasso", "ElasticNet", "RandomForest", "DeepNN")

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 080: Plotting Model Results (scRNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

# Create output directories
create_output_directories()

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input directories from config
linear_results_dir <- path.expand(OUTPUT_MODELS_LINEAR_DIR)
nn_results_dir <- path.expand(OUTPUT_MODELS_NN_DIR)
out_dir <- path.expand(OUTPUT_FIGURES_DIR)

# Gene set to process
GENE_SET <- MODEL_GENE_SET

cat("Configuration:\n")
cat(sprintf("  Sample: %s\n", SAMPLE_NAME))
cat(sprintf("  Dimensionality Reduction: %s\n", DIM_REDUCTION_METHOD))
cat(sprintf("  Gene set: %s\n", GENE_SET))
cat(sprintf("  Linear models dir: %s\n", linear_results_dir))
cat(sprintf("  Neural network dir: %s\n", nn_results_dir))
cat(sprintf("  Output figures dir: %s\n", out_dir))
cat(sprintf("  Feature type: HVG expression only (RNA-only)\n"))

# Create output directory
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Safely read CSV files
safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
}

# List all gene directories
list_gene_dirs <- function(base) {
  dirs <- list.dirs(base, recursive = FALSE, full.names = TRUE)
  dirs[file.exists(dirs)]
}

# Outlier removal helper (removes values beyond 1.5*IQR)
remove_outliers <- function(df, value_col) {
  col_sym <- rlang::sym(value_col)
  df %>%
    group_by(Model) %>%
    mutate(
      Q1 = quantile(!!col_sym, 0.25, na.rm = TRUE),
      Q3 = quantile(!!col_sym, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      lower_bound = Q1 - 1.5 * IQR,
      upper_bound = Q3 + 1.5 * IQR
    ) %>%
    filter(!!col_sym >= lower_bound & !!col_sym <= upper_bound) %>%
    select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound) %>%
    ungroup()
}

# Quantile normalization helper
quantile_normalize <- function(x) {
  n <- length(x)
  if (n == 0 || all(is.na(x))) return(x)
  ranks <- rank(x, ties.method = "average", na.last = "keep")
  qnorm((ranks - 0.5) / n)
}

# ============================================================================
# LOAD DATA - Process each gene set
# ============================================================================

# Auto-detect available gene sets
available_gene_sets <- c()
for (potential_set in c("HVG", "Random_genes")) {
  linear_path <- file.path(linear_results_dir, potential_set)
  nn_path <- file.path(nn_results_dir, potential_set)
  if (dir.exists(linear_path) || dir.exists(nn_path)) {
    available_gene_sets <- c(available_gene_sets, potential_set)
  }
}

# Also check for nn subdirectories (One_hidden_layer, Two_hidden_layer, etc.)
nn_arch_dirs <- list.dirs(nn_results_dir, recursive = FALSE, full.names = FALSE)
nn_arch_dirs <- nn_arch_dirs[grepl("hidden_layer", nn_arch_dirs, ignore.case = TRUE)]
if (length(nn_arch_dirs) > 0) {
  for (arch_dir in nn_arch_dirs) {
    for (potential_set in c("HVG", "Random_genes")) {
      nn_path <- file.path(nn_results_dir, arch_dir, potential_set)
      if (dir.exists(nn_path) && !(potential_set %in% available_gene_sets)) {
        available_gene_sets <- c(available_gene_sets, potential_set)
      }
    }
  }
}

# Filter based on config or use all available
if (GENE_SET == "HVG") {
  gene_sets_to_process <- intersect(available_gene_sets, "HVG")
} else if (GENE_SET == "Random_genes") {
  gene_sets_to_process <- intersect(available_gene_sets, "Random_genes")
} else {
  gene_sets_to_process <- available_gene_sets
}

if (length(gene_sets_to_process) == 0) {
  # Try without gene set subdirectory
  gene_sets_to_process <- c("default")
  cat("No gene set subdirectories found, using default paths\n")
}

cat(sprintf("Gene sets to process: %s\n", paste(gene_sets_to_process, collapse = ", ")))

for (gene_set_name in gene_sets_to_process) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf(">>> Processing Gene Set: %s <<<\n", gene_set_name))
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  # Set paths for this gene set
  if (gene_set_name == "default") {
    linear_dir <- linear_results_dir
    # Find nn directory (may have architecture subdirectory)
    if (length(nn_arch_dirs) > 0) {
      nn_dir <- file.path(nn_results_dir, nn_arch_dirs[1])
    } else {
      nn_dir <- nn_results_dir
    }
    set_out_dir <- out_dir
  } else {
    linear_dir <- file.path(linear_results_dir, gene_set_name)
    # Find nn directory with architecture
    if (length(nn_arch_dirs) > 0) {
      nn_dir <- file.path(nn_results_dir, nn_arch_dirs[1], gene_set_name)
    } else {
      nn_dir <- file.path(nn_results_dir, gene_set_name)
    }
    set_out_dir <- file.path(out_dir, gene_set_name)
  }
  
  dir.create(set_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ============================================================================
  # LOAD METRICS FROM BOTH MODEL TYPES
  # ============================================================================
  
  cat("\n=== Loading Metrics ===\n")
  
  # Load linear/RF metrics
  linear_metrics_file <- file.path(linear_dir, "aggregated_metrics.csv")
  if (!file.exists(linear_metrics_file)) {
    linear_metrics_file <- file.path(linear_dir, "all_genes_combined_metrics.csv")
  }
  
  linear_metrics <- NULL
  if (file.exists(linear_metrics_file)) {
    linear_metrics <- safe_read_csv(linear_metrics_file)
    cat(sprintf("Loaded %d entries from linear/RF models\n", nrow(linear_metrics)))
  } else {
    cat("Note: Linear model metrics not found at:", linear_metrics_file, "\n")
  }
  
  # Load DeepNN metrics
  nn_metrics_file <- file.path(nn_dir, "all_genes_combined_metrics.csv")
  nn_metrics <- NULL
  if (file.exists(nn_metrics_file)) {
    nn_metrics <- safe_read_csv(nn_metrics_file)
    cat(sprintf("Loaded %d entries from DeepNN models\n", nrow(nn_metrics)))
  } else {
    cat("Note: Neural network metrics not found at:", nn_metrics_file, "\n")
  }
  
  # Combine metrics
  if (!is.null(linear_metrics) && !is.null(nn_metrics)) {
    # Ensure consistent column names
    common_cols <- intersect(names(linear_metrics), names(nn_metrics))
    metrics <- bind_rows(
      linear_metrics %>% select(all_of(common_cols)),
      nn_metrics %>% select(all_of(common_cols))
    )
    cat(sprintf("Combined metrics: %d total entries\n", nrow(metrics)))
  } else if (!is.null(linear_metrics)) {
    metrics <- linear_metrics
    cat("Using linear/RF metrics only\n")
  } else if (!is.null(nn_metrics)) {
    metrics <- nn_metrics
    cat("Using DeepNN metrics only\n")
  } else {
    cat("ERROR: No metrics found! Skipping this gene set.\n")
    next
  }
  
  # ============================================================================
  # LOAD COEFFICIENTS (LINEAR MODELS ONLY)
  # ============================================================================
  
  cat("\n=== Loading Coefficients ===\n")
  
  gene_dirs_linear <- list_gene_dirs(linear_dir)
  cat(sprintf("Found %d gene directories for linear models\n", length(gene_dirs_linear)))
  
  coeffs_list <- list()
  
  for (gd in gene_dirs_linear) {
    gene_name <- basename(gd)
    
    # Load OLS coefficients
    ols_file <- file.path(gd, "coefficients_ols.csv")
    if (file.exists(ols_file)) {
      ols <- safe_read_csv(ols_file)
      if (!is.null(ols)) {
        ols$Gene <- gene_name
        ols$Model <- "OLS"
        coeffs_list[[length(coeffs_list) + 1]] <- ols
      }
    }
    
    # Load Ridge coefficients
    ridge_file <- file.path(gd, "coefficients_ridge.csv")
    if (file.exists(ridge_file)) {
      ridge <- safe_read_csv(ridge_file)
      if (!is.null(ridge)) {
        ridge$Gene <- gene_name
        ridge$Model <- "Ridge"
        coeffs_list[[length(coeffs_list) + 1]] <- ridge
      }
    }
    
    # Load Lasso coefficients
    lasso_file <- file.path(gd, "coefficients_lasso.csv")
    if (file.exists(lasso_file)) {
      lasso <- safe_read_csv(lasso_file)
      if (!is.null(lasso)) {
        lasso$Gene <- gene_name
        lasso$Model <- "Lasso"
        coeffs_list[[length(coeffs_list) + 1]] <- lasso
      }
    }
    
    # Load ElasticNet coefficients
    enet_file <- file.path(gd, "coefficients_enet.csv")
    if (file.exists(enet_file)) {
      enet <- safe_read_csv(enet_file)
      if (!is.null(enet)) {
        enet$Gene <- gene_name
        enet$Model <- "ElasticNet"
        coeffs_list[[length(coeffs_list) + 1]] <- enet
      }
    }
    
    # Load RF importance
    rf_file <- file.path(gd, "rf_importance.csv")
    if (file.exists(rf_file)) {
      rf <- safe_read_csv(rf_file)
      if (!is.null(rf)) {
        rf$Gene <- gene_name
        rf$Model <- "RandomForest"
        if ("Importance" %in% names(rf)) {
          rf$Coefficient <- rf$Importance
        }
        if (!"Original_Feature" %in% names(rf) && "Feature" %in% names(rf)) {
          rf$Original_Feature <- rf$Feature
        }
        coeffs_list[[length(coeffs_list) + 1]] <- rf
      }
    }
  }
  
  # Combine coefficients
  coeffs_df <- bind_rows(coeffs_list)
  if (nrow(coeffs_df) > 0) {
    cat(sprintf("Loaded %d coefficient entries from %d genes\n", 
                nrow(coeffs_df), n_distinct(coeffs_df$Gene)))
    
    coeffs_df <- coeffs_df %>%
      mutate(
        Coefficient = as.numeric(Coefficient),
        abs_coeff = abs(Coefficient)
      ) %>%
      filter(!is.na(Coefficient))
  }
  
  # ============================================================================
  # PLOT 1: METRICS DISTRIBUTIONS (ALL MODELS COMBINED)
  # ============================================================================
  
  cat("\n=== Creating metrics distribution plots (All Models) ===\n")
  
  if (nrow(metrics) > 0) {
    
    # Add model type classification
    metrics <- metrics %>%
      mutate(ModelType = case_when(
        Model %in% c("OLS", "Ridge", "Lasso", "ElasticNet") ~ "Linear",
        Model == "RandomForest" ~ "Tree-Based",
        Model == "DeepNN" ~ "Deep Learning",
        TRUE ~ "Other"
      ))
    
    # Filter for test set only
    test_metrics <- metrics %>% filter(Split == "Test")
    
    # Ensure model order for plotting
    test_metrics <- test_metrics %>%
      mutate(Model = factor(Model, levels = model_order))
    
    # --- R² Distribution ---
    r2_clean <- test_metrics %>%
      group_by(Model) %>%
      filter(R2 > quantile(R2, 0.05, na.rm = TRUE) & 
               R2 < quantile(R2, 0.95, na.rm = TRUE)) %>%
      ungroup()
    
    r2_medians <- r2_clean %>%
      group_by(Model) %>%
      summarise(median_val = median(R2, na.rm = TRUE), .groups = "drop")
    
    p_r2 <- r2_clean %>%
      left_join(r2_medians, by = "Model") %>%
      ggplot(aes(x = Model, y = R2, fill = Model)) +
      geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.3) +
      geom_boxplot(width = 0.25, alpha = 0.9, outlier.shape = NA, 
                   color = "black", linewidth = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      geom_text(aes(y = median_val, label = sprintf("%.2f", median_val)),
                vjust = -0.8, size = 3.5, fontface = "bold") +
      scale_fill_manual(values = model_palette) +
      coord_cartesian(ylim = c(-1, 1)) +
      labs(title = expression(bold(R^2 ~ "Distribution by Model (Test Set)")),
           subtitle = sprintf("%s | %s | RNA-only (HVG features)", SAMPLE_NAME, DIMRED_METHOD_SUFFIX),
           x = NULL, y = expression(R^2)) +
      theme_publication() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
    
    save_publication_plot(p_r2, file.path(set_out_dir, "01a_metrics_R2_all_models"),
                          width = 10, height = 7)
    
    # --- RMSE Distribution ---
    rmse_clean <- test_metrics %>%
      group_by(Model) %>%
      filter(RMSE < quantile(RMSE, 0.95, na.rm = TRUE)) %>%
      ungroup()
    
    rmse_medians <- rmse_clean %>%
      group_by(Model) %>%
      summarise(median_val = median(RMSE, na.rm = TRUE), .groups = "drop")
    
    p_rmse <- rmse_clean %>%
      left_join(rmse_medians, by = "Model") %>%
      ggplot(aes(x = Model, y = RMSE, fill = Model)) +
      geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.3) +
      geom_boxplot(width = 0.25, alpha = 0.9, outlier.shape = NA, 
                   color = "black", linewidth = 0.5) +
      geom_text(aes(y = median_val, label = sprintf("%.2f", median_val)),
                vjust = -0.8, size = 3.5, fontface = "bold") +
      scale_fill_manual(values = model_palette) +
      labs(title = "RMSE Distribution by Model (Test Set)",
           subtitle = sprintf("%s | %s | Lower is better", SAMPLE_NAME, DIMRED_METHOD_SUFFIX),
           x = NULL, y = "Root Mean Squared Error") +
      theme_publication() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
    
    save_publication_plot(p_rmse, file.path(set_out_dir, "01b_metrics_RMSE_all_models"),
                          width = 10, height = 7)
    
    # --- Spearman Correlation ---
    spear_clean <- test_metrics %>%
      filter(!is.na(Spearman)) %>%
      group_by(Model) %>%
      ungroup()
    
    spear_medians <- spear_clean %>%
      group_by(Model) %>%
      summarise(median_val = median(Spearman, na.rm = TRUE), .groups = "drop")
    
    p_spearman <- spear_clean %>%
      left_join(spear_medians, by = "Model") %>%
      ggplot(aes(x = Model, y = Spearman, fill = Model)) +
      geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.3) +
      geom_boxplot(width = 0.25, alpha = 0.9, outlier.shape = NA, 
                   color = "black", linewidth = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      geom_text(aes(y = median_val, label = sprintf("%.2f", median_val)),
                vjust = -0.8, size = 3.5, fontface = "bold") +
      scale_fill_manual(values = model_palette) +
      labs(title = "Spearman Correlation by Model (Test Set)",
           subtitle = sprintf("%s | %s | Rank-based correlation", SAMPLE_NAME, DIMRED_METHOD_SUFFIX),
           x = NULL, y = expression("Spearman " * rho)) +
      theme_publication() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
    
    save_publication_plot(p_spearman, file.path(set_out_dir, "01c_metrics_Spearman_all_models"),
                          width = 10, height = 7)
    
    # --- Combined metrics plot ---
    combined_metrics <- (p_r2 | p_rmse) / p_spearman +
      plot_annotation(
        title = "Model Performance Comparison (RNA-only)",
        subtitle = sprintf("%s | %s | All Models on Test Set", SAMPLE_NAME, DIMRED_METHOD_SUFFIX),
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"))
      )
    save_publication_plot(combined_metrics, file.path(set_out_dir, "01_metrics_all_models_combined"),
                          width = 16, height = 14)
    
    cat("  ✓ Saved metrics distribution plots\n")
    
    # ============================================================================
    # PLOT 2: MODEL TYPE COMPARISON (Grouped)
    # ============================================================================
    
    cat("\n=== Creating model type comparison ===\n")
    
    # Summarize by model type
    type_summary <- test_metrics %>%
      group_by(ModelType) %>%
      summarise(
        N_Predictions = n(),
        Mean_R2 = mean(R2, na.rm = TRUE),
        Median_R2 = median(R2, na.rm = TRUE),
        SD_R2 = sd(R2, na.rm = TRUE),
        Mean_RMSE = mean(RMSE, na.rm = TRUE),
        Mean_Spearman = mean(Spearman, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(Mean_R2))
    
    write_csv(type_summary, file.path(set_out_dir, "model_type_summary.csv"))
    
    # Plot by model type
    p_type_r2 <- test_metrics %>%
      ggplot(aes(x = ModelType, y = R2, fill = ModelType)) +
      geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.3) +
      geom_boxplot(width = 0.25, alpha = 0.9, outlier.shape = NA, color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      scale_fill_manual(values = model_type_palette) +
      labs(title = expression(bold(R^2 ~ "by Model Type")),
           subtitle = sprintf("%s | Linear vs Tree-Based vs Deep Learning", SAMPLE_NAME),
           x = NULL, y = expression(R^2)) +
      theme_publication() +
      theme(legend.position = "none")
    
    save_publication_plot(p_type_r2, file.path(set_out_dir, "02_model_type_comparison"),
                          width = 9, height = 7)
    
    cat("  ✓ Saved model type comparison\n")
    
    # ============================================================================
    # PLOT 3: MODEL WINS (WHICH MODEL IS BEST PER GENE?)
    # ============================================================================
    
    cat("\n=== Creating model wins analysis ===\n")
    
    # Find best model per gene based on R² on test set
    best_models <- test_metrics %>%
      group_by(Gene) %>%
      arrange(desc(R2)) %>%
      slice_head(n = 1) %>%
      ungroup()
    
    # Count wins per model
    model_wins <- best_models %>%
      count(Model, name = "N_Genes_Won") %>%
      arrange(desc(N_Genes_Won)) %>%
      mutate(Model = factor(Model, levels = model_order))
    
    p_model_wins <- model_wins %>%
      ggplot(aes(x = reorder(Model, N_Genes_Won), y = N_Genes_Won, fill = Model)) +
      geom_col(alpha = 0.8, color = "black", linewidth = 0.5) +
      geom_text(aes(label = N_Genes_Won), hjust = -0.3, size = 4, fontface = "bold") +
      scale_fill_manual(values = model_palette) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      coord_flip() +
      labs(title = "Best Model per Gene",
           subtitle = sprintf("%s | Winner = highest R² on test set", SAMPLE_NAME),
           x = NULL, y = "Number of Genes Won") +
      theme_publication() +
      theme(legend.position = "none")
    
    save_publication_plot(p_model_wins, file.path(set_out_dir, "03a_model_wins_per_gene"),
                          width = 10, height = 6)
    
    # Distribution of winning R² values
    p_winning_r2 <- best_models %>%
      mutate(Model = factor(Model, levels = model_order)) %>%
      ggplot(aes(x = Model, y = R2, fill = Model)) +
      geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.3) +
      geom_boxplot(width = 0.25, alpha = 0.9, outlier.shape = NA, color = "black") +
      geom_jitter(width = 0.12, alpha = 0.3, size = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      scale_fill_manual(values = model_palette) +
      labs(title = expression(bold(R^2 ~ "When Each Model Wins")),
           subtitle = sprintf("%s | Distribution of best R² values", SAMPLE_NAME),
           x = NULL, y = expression(R^2 ~ "(Test Set)")) +
      theme_publication() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    save_publication_plot(p_winning_r2, file.path(set_out_dir, "03b_winning_r2_distribution"),
                          width = 10, height = 7)
    
    # Combined wins plot
    combined_wins <- p_model_wins | p_winning_r2
    save_publication_plot(combined_wins, file.path(set_out_dir, "03_model_wins_combined"),
                          width = 16, height = 7)
    
    write_csv(model_wins, file.path(set_out_dir, "model_wins_summary.csv"))
    write_csv(best_models, file.path(set_out_dir, "best_model_per_gene.csv"))
    
    cat("  ✓ Saved model wins analysis\n")
    
    # ============================================================================
    # PLOT 4: DEEP LEARNING vs BEST LINEAR/RF MODEL
    # ============================================================================
    
    if ("DeepNN" %in% unique(test_metrics$Model)) {
      
      cat("\n=== Creating DeepNN vs Traditional Models comparison ===\n")
      
      # Get DeepNN performance and best non-DL model per gene
      dl_comparison <- test_metrics %>%
        mutate(is_DL = Model == "DeepNN") %>%
        group_by(Gene, is_DL) %>%
        summarise(
          Best_R2 = max(R2, na.rm = TRUE),
          Best_Model = Model[which.max(R2)],
          .groups = "drop"
        ) %>%
        pivot_wider(
          id_cols = Gene,
          names_from = is_DL,
          values_from = c(Best_R2, Best_Model),
          names_glue = "{.value}_{ifelse(is_DL, 'DL', 'Traditional')}"
        )

      # Normalize pivoted column names to stable identifiers (handles True/False or DL/Traditional variants)
      if ("Best_R2_DL" %in% names(dl_comparison)) {
        dl_comparison <- dl_comparison %>% rename(DeepNN_R2 = Best_R2_DL)
      } else if ("Best_R2_TRUE" %in% names(dl_comparison)) {
        dl_comparison <- dl_comparison %>% rename(DeepNN_R2 = Best_R2_TRUE)
      } else {
        dl_comparison$DeepNN_R2 <- NA_real_
      }

      if ("Best_R2_Traditional" %in% names(dl_comparison)) {
        dl_comparison <- dl_comparison %>% rename(Traditional_R2 = Best_R2_Traditional)
      } else if ("Best_R2_FALSE" %in% names(dl_comparison)) {
        dl_comparison <- dl_comparison %>% rename(Traditional_R2 = Best_R2_FALSE)
      } else {
        dl_comparison$Traditional_R2 <- NA_real_
      }

      if ("Best_Model_Traditional" %in% names(dl_comparison)) {
        dl_comparison <- dl_comparison %>% rename(Best_Traditional = Best_Model_Traditional)
      } else if ("Best_Model_FALSE" %in% names(dl_comparison)) {
        dl_comparison <- dl_comparison %>% rename(Best_Traditional = Best_Model_FALSE)
      } else if ("Best_Model_Traditional" %in% names(dl_comparison) == FALSE && !("Best_Traditional" %in% names(dl_comparison))) {
        dl_comparison$Best_Traditional <- NA_character_
      }

      dl_comparison <- dl_comparison %>%
        mutate(
          R2_Diff = DeepNN_R2 - Traditional_R2,
          DL_Wins = ifelse(is.na(R2_Diff), FALSE, R2_Diff > 0)
        )
      
      # Scatter plot: DeepNN R² vs Best Traditional R²
      p_dl_vs_trad <- dl_comparison %>%
        ggplot(aes(x = Traditional_R2, y = DeepNN_R2, color = DL_Wins)) +
        geom_point(alpha = 0.6, size = 2) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
        scale_color_manual(values = c("FALSE" = "#0072B2", "TRUE" = "#CC79A7"),
                           labels = c("Traditional Better", "DeepNN Better"),
                           name = "Winner") +
        labs(title = "DeepNN vs Best Traditional Model",
             subtitle = sprintf("%s | DeepNN wins for %.1f%% of genes", SAMPLE_NAME,
                               100 * mean(dl_comparison$DL_Wins, na.rm = TRUE)),
             x = expression("Best Traditional Model " * R^2),
             y = expression("DeepNN " * R^2)) +
        theme_publication() +
        theme(legend.position = "bottom") +
        coord_fixed()
      
      save_publication_plot(p_dl_vs_trad, file.path(set_out_dir, "04a_deepnn_vs_traditional"),
                            width = 9, height = 9)
      
      # R² improvement distribution
      p_r2_improvement <- dl_comparison %>%
        ggplot(aes(x = R2_Diff)) +
        geom_histogram(bins = 30, fill = "#CC79A7", color = "black", alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "#D55E00", linewidth = 1) +
        geom_vline(xintercept = median(dl_comparison$R2_Diff, na.rm = TRUE),
                   linetype = "solid", color = "#0072B2", linewidth = 1) +
        annotate("text", x = median(dl_comparison$R2_Diff, na.rm = TRUE), 
                 y = Inf, label = sprintf("Median: %.3f", median(dl_comparison$R2_Diff, na.rm = TRUE)),
                 vjust = 2, hjust = -0.1, color = "#0072B2", fontface = "bold", size = 4) +
        labs(title = expression(bold(R^2 ~ "Improvement: DeepNN - Best Traditional")),
             subtitle = sprintf("%s | Positive = DeepNN better", SAMPLE_NAME),
             x = expression(R^2 ~ "Difference (DeepNN - Traditional)"),
             y = "Number of Genes",
             caption = "Red dashed = zero, Blue = median") +
        theme_publication()
      
      save_publication_plot(p_r2_improvement, file.path(set_out_dir, "04b_deepnn_improvement_distribution"),
                            width = 10, height = 7)
      
      # Combined DL comparison
      combined_dl <- p_dl_vs_trad | p_r2_improvement
      save_publication_plot(combined_dl, file.path(set_out_dir, "04_deepnn_comparison_combined"),
                            width = 16, height = 8)
      
      write_csv(dl_comparison, file.path(set_out_dir, "deepnn_vs_traditional_comparison.csv"))
      
      cat("  ✓ Saved DeepNN comparison plots\n")
    }
    
    # ============================================================================
    # PLOT 5: OVERFITTING ANALYSIS (ALL MODELS)
    # ============================================================================
    
    cat("\n=== Creating overfitting analysis ===\n")
    
    train_test_gap <- metrics %>%
      filter(Split %in% c("Train", "Test")) %>%
      select(Gene, Model, Split, R2, RMSE) %>%
      pivot_wider(names_from = Split, values_from = c(R2, RMSE)) %>%
      mutate(
        R2_Gap = R2_Train - R2_Test,
        RMSE_Gap = RMSE_Test - RMSE_Train,
        Model = factor(Model, levels = model_order)
      )
    
    # R² train vs test scatter
    p_r2_traintest <- train_test_gap %>%
      ggplot(aes(x = R2_Train, y = R2_Test, color = Model)) +
      geom_point(alpha = 0.5, size = 1.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
      scale_color_manual(values = model_palette) +
      facet_wrap(~ Model, ncol = 3) +
      labs(title = expression(bold("Train vs Test " * R^2 * ": Overfitting Check")),
           subtitle = sprintf("%s | Points below diagonal = overfitting", SAMPLE_NAME),
           x = expression(R^2 ~ "(Train)"), y = expression(R^2 ~ "(Test)")) +
      theme_publication(base_size = 11) +
      theme(legend.position = "none")
    
    save_publication_plot(p_r2_traintest, file.path(set_out_dir, "05a_train_vs_test_r2"),
                          width = 14, height = 10)
    
    # Overfitting gap by model
    gap_medians <- train_test_gap %>%
      group_by(Model) %>%
      summarise(median_gap = median(R2_Gap, na.rm = TRUE), .groups = "drop")
    
    p_overfit_gap <- train_test_gap %>%
      left_join(gap_medians, by = "Model") %>%
      ggplot(aes(x = Model, y = R2_Gap, fill = Model)) +
      geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.3) +
      geom_boxplot(width = 0.25, alpha = 0.9, outlier.shape = NA, color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
      geom_text(aes(y = median_gap, label = sprintf("%.2f", median_gap)),
                vjust = -0.8, size = 3, fontface = "bold") +
      scale_fill_manual(values = model_palette) +
      labs(title = expression(bold("Overfitting Gap: " * R^2 * " Train - " * R^2 * " Test")),
           subtitle = sprintf("%s | Higher values = more overfitting (0 = perfect generalization)", SAMPLE_NAME),
           x = NULL, y = expression(R^2 ~ "Gap (Train - Test)")) +
      theme_publication() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
    
    save_publication_plot(p_overfit_gap, file.path(set_out_dir, "05b_overfitting_gap"),
                          width = 10, height = 7)
    
    # Combined
    combined_overfit <- p_r2_traintest / p_overfit_gap + plot_layout(heights = c(2, 1))
    save_publication_plot(combined_overfit, file.path(set_out_dir, "05_overfitting_analysis_combined"),
                          width = 14, height = 14)
    
    write_csv(train_test_gap, file.path(set_out_dir, "train_test_gap_analysis.csv"))
    
    cat("  ✓ Saved overfitting analysis\n")
    
    # ============================================================================
    # PLOT 6: GENE-LEVEL PERFORMANCE HEATMAP
    # ============================================================================
    
    cat("\n=== Creating gene-level heatmap ===\n")
    
    heatmap_data <- test_metrics %>%
      select(Gene, Model, R2) %>%
      mutate(Model = factor(Model, levels = model_order)) %>%
      pivot_wider(names_from = Model, values_from = R2)
    
    if (nrow(heatmap_data) <= 100) {
      heatmap_long <- heatmap_data %>%
        pivot_longer(-Gene, names_to = "Model", values_to = "R2") %>%
        mutate(Model = factor(Model, levels = model_order))
      
      p_heatmap <- heatmap_long %>%
        ggplot(aes(x = Model, y = Gene, fill = R2)) +
        geom_tile(color = "white", linewidth = 0.3) +
        scale_fill_viridis_c(option = "plasma", name = expression(R^2), limits = c(-0.5, 1)) +
        labs(title = "Gene-Level Model Performance (Test Set)",
             subtitle = sprintf("%s | All models compared side-by-side", SAMPLE_NAME),
             x = NULL, y = NULL) +
        theme_publication(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
              axis.text.y = element_text(size = 6),
              panel.grid = element_blank())
      
      save_publication_plot(p_heatmap, file.path(set_out_dir, "06_gene_level_heatmap"),
                            width = 9, height = max(7, nrow(heatmap_data) * 0.15))
      
      cat("  ✓ Saved gene-level heatmap\n")
    } else {
      cat("  ⊘ Too many genes for heatmap (", nrow(heatmap_data), "), skipping\n")
    }
    
    # ============================================================================
    # PLOT 7: COEFFICIENT/IMPORTANCE ANALYSIS (LINEAR MODELS + RF)
    # ============================================================================
    
    if (nrow(coeffs_df) > 0) {
      
      cat("\n=== Creating coefficient analysis ===\n")
      
      # Distribution by model
      coeffs_df_hist <- coeffs_df %>%
        group_by(Model) %>%
        filter(abs_coeff <= quantile(abs_coeff, 0.99, na.rm = TRUE)) %>%
        ungroup()
      
      p_coef_dist <- coeffs_df_hist %>%
        ggplot(aes(x = abs_coeff, fill = Model)) +
        geom_histogram(bins = 50, alpha = 0.7, color = "white", linewidth = 0.2) +
        facet_wrap(~Model, scales = "free") +
        scale_x_log10(labels = label_scientific()) +
        scale_fill_manual(values = model_palette) +
        labs(title = "Coefficient/Importance Distribution (Linear + RF)",
             subtitle = sprintf("%s | Log-scaled x-axis | RNA-only features", SAMPLE_NAME),
             x = "Absolute Coefficient/Importance (log scale)", y = "Count") +
        theme_publication(base_size = 11) +
        theme(legend.position = "none")
      
      save_publication_plot(p_coef_dist, file.path(set_out_dir, "07a_coefficient_distribution"),
                            width = 14, height = 9)
      
      # Top features overall
      top_features_overall <- coeffs_df %>%
        group_by(Model) %>%
        slice_max(order_by = abs_coeff, n = 20) %>%
        ungroup()
      
      write_csv(top_features_overall, file.path(set_out_dir, "top20_features_by_model.csv"))
      
      cat("  ✓ Saved coefficient analysis plots\n")
    }
    
    # ============================================================================
    # PLOT 8: SCATTER PLOTS (ALL MODELS)
    # ============================================================================
    
    cat("\n=== Creating scatter plots ===\n")
    
    scatter_dir <- file.path(set_out_dir, "scatter_plots")
    dir.create(scatter_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Select random genes
    set.seed(42)
    all_genes <- unique(metrics$Gene)
    selected_genes <- sample(all_genes, min(10, length(all_genes)))
    
    for (gene in selected_genes) {
      # Load predictions from both sources
      all_preds <- list()
      
      # Linear/RF predictions
      linear_pred_file <- file.path(linear_dir, gene, "predictions_test.csv")
      if (file.exists(linear_pred_file)) {
        linear_preds <- safe_read_csv(linear_pred_file)
        if (!is.null(linear_preds)) all_preds[["linear"]] <- linear_preds
      }
      
      # DeepNN predictions
      nn_pred_file <- file.path(nn_dir, gene, "predictions_test.csv")
      if (file.exists(nn_pred_file)) {
        nn_preds <- safe_read_csv(nn_pred_file)
        if (!is.null(nn_preds)) all_preds[["nn"]] <- nn_preds
      }
      
      if (length(all_preds) == 0) next
      
      # Combine predictions
      if (length(all_preds) == 2) {
        # Robustly determine join column between linear and nn prediction tables
        lin_names <- names(all_preds$linear)
        nn_names <- names(all_preds$nn)

        # Prefer exact common names
        common <- intersect(lin_names, nn_names)
        if (length(common) > 0) {
          id_col <- common[1]
        } else {
          # Try case-insensitive matching
          lin_low <- tolower(lin_names)
          nn_low <- tolower(nn_names)
          shared_low <- intersect(lin_low, nn_low)
          if (length(shared_low) > 0) {
            # Use the linear column name as canonical id_col and rename nn column to match
            matched_low <- shared_low[1]
            id_col <- lin_names[which(lin_low == matched_low)[1]]
            nn_idx <- which(nn_low == matched_low)[1]
            names(all_preds$nn)[nn_idx] <- id_col
          } else {
            # Try common candidate names
            candidates <- c("Cell_ID", "cell", "cell_id", "Cell", "cellID", "cellId")
            cand_lin <- intersect(candidates, lin_names)
            cand_nn <- intersect(candidates, nn_names)
            if (length(intersect(cand_lin, cand_nn)) > 0) {
              id_col <- intersect(cand_lin, cand_nn)[1]
            } else if (length(cand_lin) > 0 && length(cand_nn) > 0) {
              # Rename nn candidate to linear candidate
              id_col <- cand_lin[1]
              names(all_preds$nn)[which(nn_names == cand_nn[1])[1]] <- id_col
            } else {
              warning(sprintf("Cannot determine join column between linear and nn predictions for gene %s — skipping scatter plot.", gene))
              next
            }
          }
        }

        pred_data <- all_preds$linear %>%
          left_join(all_preds$nn %>% select(any_of(c(id_col, "DeepNN"))), by = id_col)
      } else if ("linear" %in% names(all_preds)) {
        pred_data <- all_preds$linear
      } else {
        pred_data <- all_preds$nn
      }
      
      # Find actual column
      actual_col <- intersect(names(pred_data), c("Actual", "actual", "y"))
      if (length(actual_col) == 0) next
      
      pred_data <- pred_data %>% rename(Actual = !!actual_col[1])
      
      # Get model columns
      model_cols <- intersect(names(pred_data), c("OLS", "Ridge", "Lasso", "ElasticNet", "RandomForest", "DeepNN"))
      if (length(model_cols) == 0) next
      
      # Convert to long format
      pred_long <- pred_data %>%
        select(Actual, any_of(model_cols)) %>%
        pivot_longer(cols = -Actual, names_to = "Model", values_to = "Predicted") %>%
        filter(!is.na(Actual) & !is.na(Predicted)) %>%
        mutate(Model = factor(Model, levels = model_order))
      
      # Apply quantile normalization
      preds_norm <- pred_long %>%
        group_by(Model) %>%
        mutate(
          Actual_QN = quantile_normalize(Actual),
          Predicted_QN = quantile_normalize(Predicted)
        ) %>%
        ungroup()
      
      # Get metrics for this gene
      gene_metrics <- metrics %>%
        filter(Gene == gene, Split == "Test") %>%
        select(Model, R2, RMSE, Spearman) %>%
        mutate(Model = factor(Model, levels = model_order))
      
      preds_with_metrics <- preds_norm %>%
        left_join(gene_metrics, by = "Model")
      
      # Create scatter plot
      p_scatter <- preds_with_metrics %>%
        ggplot(aes(x = Actual_QN, y = Predicted_QN)) +
        geom_point(alpha = 0.2, size = 0.5, color = "#0072B2") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_smooth(method = "lm", se = TRUE, linewidth = 0.6, 
                    alpha = 0.15, color = "#D55E00", fill = "#D55E00") +
        geom_label(aes(label = sprintf("R² = %.3f\nρ = %.3f", R2, Spearman)),
                   x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05,
                   size = 2.5, fontface = "plain", label.size = 0.1,
                   fill = "white", alpha = 0.9,
                   data = . %>% distinct(Model, R2, Spearman)) +
        facet_wrap(~ Model, scales = "free", ncol = 3) +
        labs(title = sprintf("Gene: %s (All Models)", gene),
             subtitle = sprintf("%s | Predicted vs Actual Expression (Quantile Normalized)", SAMPLE_NAME),
             x = "Actual Expression (QN)", y = "Predicted Expression (QN)") +
        theme_publication(base_size = 10) +
        theme(panel.grid = element_line(color = "gray95"))
      
      safe_gene_name <- str_replace_all(gene, "[^A-Za-z0-9_-]", "_")
      save_publication_plot(p_scatter, file.path(scatter_dir, sprintf("scatter_%s", safe_gene_name)),
                            width = 11, height = 9)
      
      cat(sprintf("  ✓ Scatter plot: %s\n", gene))
    }
    
    # ============================================================================
    # PLOT 9: SUMMARY STATISTICS TABLE
    # ============================================================================
    
    cat("\n=== Creating summary statistics ===\n")
    
    summary_stats <- test_metrics %>%
      group_by(Model) %>%
      summarise(
        N_Genes = n_distinct(Gene),
        Mean_R2 = mean(R2, na.rm = TRUE),
        Median_R2 = median(R2, na.rm = TRUE),
        SD_R2 = sd(R2, na.rm = TRUE),
        Mean_RMSE = mean(RMSE, na.rm = TRUE),
        Median_RMSE = median(RMSE, na.rm = TRUE),
        Mean_Spearman = mean(Spearman, na.rm = TRUE),
        Median_Spearman = median(Spearman, na.rm = TRUE),
        Pct_Positive_R2 = 100 * mean(R2 > 0, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(Median_R2))
    
    write_csv(summary_stats, file.path(set_out_dir, "summary_statistics.csv"))
    
    # Print summary
    cat("\n=== Performance Summary (Test Set) ===\n")
    print(summary_stats, n = Inf)
    
    cat("\n  ✓ Saved summary statistics\n")
  }
  
  cat(sprintf("\n>>> Completed %s gene set <<<\n", gene_set_name))
}

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 080 COMPLETE - PUBLICATION-READY PLOTS (scRNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat(sprintf("Total runtime: %.2f minutes\n\n", as.numeric(runtime)))

cat("Generated publication-ready plots (PNG 600dpi + PDF vector):\n")
cat("  01 - Metrics distribution (all models)\n")
cat("  02 - Model type comparison\n")
cat("  03 - Model wins per gene\n")
cat("  04 - DeepNN vs Traditional comparison\n")
cat("  05 - Overfitting analysis\n")
cat("  06 - Gene-level heatmap\n")
cat("  07 - Coefficient analysis (Linear + RF only)\n")
cat("  08 - Scatter plots (sample genes)\n")

cat("\nPlot features:\n")
cat("  - Colorblind-friendly palette (Wong 2011)\n")
cat("  - High resolution PNG (600 dpi)\n")
cat("  - Vector PDF for publication\n")
cat("  - Consistent publication theme\n")
cat("  - RNA-only (HVG features, no ATAC)\n")

cat("\nOutput directory:", path.expand(OUTPUT_FIGURES_DIR), "\n\n")

cat("=", rep("=", 70), "\n\n", sep = "")
