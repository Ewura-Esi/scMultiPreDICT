#!/usr/bin/env Rscript
# ============================================================================
# Step_082: Aggregate SHAP Feature Importance Across All Genes (scRNA-only)
# ============================================================================
# This script aggregates SHAP feature importance results from Step_081
# across all genes and creates publication-ready summary visualizations.
#
# Features:
#   - Combines per-gene SHAP importance into a single dataset
#   - Identifies globally important features (appear in many genes' top features)
#   - Creates publication-ready summary plots and heatmaps
#   - Outputs both PNG (600 dpi) and PDF (vector) formats
#
# Note: This is the RNA-only version (no ATAC vs RNA comparison)
#
# Input: Per-gene SHAP results from Step_081
# Output: Aggregated CSV, summary statistics, and publication-ready plots
#
# Usage:
#   Rscript Step_082.Aggregate_SHAP_Results_Automated.R [GENE_SET]
#
# Examples:
#   Rscript Step_082.Aggregate_SHAP_Results_Automated.R HVG
#   Rscript Step_082.Aggregate_SHAP_Results_Automated.R
# ============================================================================

# Record start time
start_time <- Sys.time()

# ============================================================================
# PARSE COMMAND LINE ARGUMENTS
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)
gene_set_arg <- if (length(args) >= 1) trimws(args[1]) else NULL

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
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  library(stringr)
  library(Cairo)
  library(scales)
})

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 082: Aggregate SHAP Feature Importance (scRNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

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

# Primary color for RNA-only plots
PRIMARY_COLOR <- "#0072B2"
SECONDARY_COLOR <- "#56B4E9"
ACCENT_COLOR <- "#D55E00"

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

# Set default theme
theme_set(theme_publication(base_size = 12))

# ============================================================================
# CONFIGURATION
# ============================================================================

# Determine gene sets to process
gene_sets_to_process <- c()

if (!is.null(gene_set_arg)) {
  gene_sets_to_process <- gene_set_arg
} else if (MODEL_GENE_SET == "both") {
  gene_sets_to_process <- c("HVG", "Custom")
} else if (MODEL_GENE_SET == "HVG") {
  gene_sets_to_process <- "HVG"
} else if (MODEL_GENE_SET %in% c("Random", "Custom")) {
  gene_sets_to_process <- "Custom"
}

cat("Configuration:\n")
cat(sprintf("  Sample: %s\n", SAMPLE_NAME))
cat(sprintf("  Gene sets to process: %s\n", paste(gene_sets_to_process, collapse = ", ")))
cat(sprintf("  Feature type: HVG expression only (RNA-only)\n"))
cat("\n")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Safely read CSV
safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
}

# Truncate feature names for display
truncate_name <- function(x, max_len = 35) {
  ifelse(nchar(x) > max_len, paste0(substr(x, 1, max_len - 3), "..."), x)
}

# ============================================================================
# PROCESS EACH GENE SET
# ============================================================================

for (gene_set_name in gene_sets_to_process) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf(">>> Processing Gene Set: %s <<<\n", gene_set_name))
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  # Paths - check for architecture subdirectory
  nn_base_dir <- path.expand(OUTPUT_MODELS_NN_DIR)
  
  # Find architecture directory
  arch_dirs <- list.dirs(nn_base_dir, recursive = FALSE, full.names = FALSE)
  arch_dirs <- arch_dirs[grepl("hidden_layer", arch_dirs, ignore.case = TRUE)]
  
  if (length(arch_dirs) > 0) {
    nn_results_dir <- file.path(nn_base_dir, arch_dirs[1], gene_set_name)
  } else {
    nn_results_dir <- file.path(nn_base_dir, gene_set_name)
  }
  
  output_dir <- file.path(path.expand(OUTPUT_FIGURES_DIR), "NEURAL_NETWORKS", gene_set_name, "SHAP_Aggregated")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat(sprintf("  NN results dir: %s\n", nn_results_dir))
  cat(sprintf("  Output dir: %s\n", output_dir))
  
  # ============================================================================
  # COLLECT SHAP RESULTS FROM ALL GENES
  # ============================================================================
  
  cat("\n=== Collecting SHAP results from all genes ===\n")
  
  # Find all gene directories with SHAP results
  gene_dirs <- list.dirs(nn_results_dir, recursive = FALSE, full.names = TRUE)
  
  shap_results <- list()
  genes_processed <- c()
  genes_missing <- c()
  
  for (gd in gene_dirs) {
    gene_name <- basename(gd)
    shap_file <- file.path(gd, "shap_feature_importance", "shap_feature_importance.csv")
    
    if (file.exists(shap_file)) {
      shap_df <- safe_read_csv(shap_file)
      
      if (!is.null(shap_df) && nrow(shap_df) > 0) {
        if (!"Gene" %in% names(shap_df)) {
          shap_df$Gene <- gene_name
        }
        shap_results[[gene_name]] <- shap_df
        genes_processed <- c(genes_processed, gene_name)
      }
    } else {
      genes_missing <- c(genes_missing, gene_name)
    }
  }
  
  cat(sprintf("  Genes with SHAP results: %d\n", length(genes_processed)))
  cat(sprintf("  Genes missing SHAP results: %d\n", length(genes_missing)))
  
  if (length(shap_results) == 0) {
    cat("  WARNING: No SHAP results found. Skipping this gene set.\n")
    cat("  Run Step_081 first to generate SHAP values.\n")
    next
  }
  
  # ============================================================================
  # COMBINE ALL SHAP RESULTS
  # ============================================================================
  
  cat("\n=== Combining SHAP results ===\n")
  
  all_shap <- bind_rows(shap_results)
  
  cat(sprintf("  Total entries: %d\n", nrow(all_shap)))
  cat(sprintf("  Unique features: %d\n", n_distinct(all_shap$Feature)))
  cat(sprintf("  Unique genes: %d\n", n_distinct(all_shap$Gene)))
  
  # Save combined results
  combined_file <- file.path(output_dir, "all_genes_shap_combined.csv")
  write_csv(all_shap, combined_file)
  cat(sprintf("  Saved combined results to: %s\n", combined_file))
  
  # ============================================================================
  # COMPUTE GLOBAL FEATURE IMPORTANCE
  # ============================================================================
  
  cat("\n=== Computing global feature importance ===\n")
  
  global_importance <- all_shap %>%
    group_by(Feature) %>%
    summarise(
      Mean_SHAP = mean(MeanAbsSHAP, na.rm = TRUE),
      Median_SHAP = median(MeanAbsSHAP, na.rm = TRUE),
      SD_SHAP = sd(MeanAbsSHAP, na.rm = TRUE),
      Max_SHAP = max(MeanAbsSHAP, na.rm = TRUE),
      N_Genes = n_distinct(Gene),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_SHAP)) %>%
    mutate(Rank = row_number())
  
  global_file <- file.path(output_dir, "global_feature_importance.csv")
  write_csv(global_importance, global_file)
  
  cat("  Top 10 globally important features:\n")
  print(head(global_importance %>% select(Rank, Feature, Mean_SHAP, N_Genes), 10))
  
  # ============================================================================
  # COMPUTE FEATURE FREQUENCY IN TOP-20
  # ============================================================================
  
  cat("\n=== Computing feature frequency in top-20 ===\n")
  
  top20_per_gene <- all_shap %>%
    group_by(Gene) %>%
    slice_max(order_by = MeanAbsSHAP, n = 20) %>%
    ungroup()
  
  feature_frequency <- top20_per_gene %>%
    count(Feature, name = "Frequency") %>%
    arrange(desc(Frequency)) %>%
    mutate(
      Pct_Genes = 100 * Frequency / n_distinct(all_shap$Gene),
      Rank = row_number()
    )
  
  freq_file <- file.path(output_dir, "feature_frequency_in_top20.csv")
  write_csv(feature_frequency, freq_file)
  
  cat("  Top 10 most frequent features in top-20:\n")
  print(head(feature_frequency, 10))
  
  # ============================================================================
  # PUBLICATION-READY PLOTS
  # ============================================================================
  
  cat("\n=== Creating publication-ready plots ===\n")
  
  n_genes <- n_distinct(all_shap$Gene)
  
  # --- Plot 1: Global Feature Importance (Top 30) ---
  p1_data <- global_importance %>%
    slice_head(n = 30) %>%
    mutate(Feature_Short = truncate_name(Feature, 35))
  
  p1 <- ggplot(p1_data, aes(x = Mean_SHAP, y = reorder(Feature_Short, Mean_SHAP))) +
    geom_col(fill = PRIMARY_COLOR, alpha = 0.85, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.4f", Mean_SHAP)), hjust = -0.1, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Global Feature Importance (Top 30)",
      subtitle = sprintf("%s | Mean |SHAP| across %d genes | RNA-only (HVG)", SAMPLE_NAME, n_genes),
      x = "Mean |SHAP Value|",
      y = NULL,
      caption = "Features ranked by average importance across all target genes"
    ) +
    theme_publication(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  
  save_publication_plot(p1, file.path(output_dir, "01_global_feature_importance"),
                        width = 12, height = 10)
  
  # --- Plot 2: Feature Frequency in Top-20 ---
  p2_data <- feature_frequency %>%
    slice_head(n = 30) %>%
    mutate(Feature_Short = truncate_name(Feature, 35))
  
  p2 <- ggplot(p2_data, aes(x = Frequency, y = reorder(Feature_Short, Frequency))) +
    geom_col(fill = ACCENT_COLOR, alpha = 0.85, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.0f%%", Pct_Genes)), hjust = -0.1, size = 3, fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Most Frequent Features in Top-20 (Across Genes)",
      subtitle = sprintf("%s | How often each feature appears in top-20 SHAP | %d genes", SAMPLE_NAME, n_genes),
      x = "Number of Genes",
      y = NULL,
      caption = "Percentage indicates the proportion of genes where this feature is in top-20"
    ) +
    theme_publication(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  
  save_publication_plot(p2, file.path(output_dir, "02_feature_frequency_top20"),
                        width = 12, height = 10)
  
  # --- Plot 3: SHAP Distribution Histogram ---
  p3 <- ggplot(all_shap, aes(x = MeanAbsSHAP)) +
    geom_histogram(bins = 50, fill = PRIMARY_COLOR, color = "white", alpha = 0.85) +
    geom_vline(xintercept = median(all_shap$MeanAbsSHAP), 
               linetype = "dashed", color = ACCENT_COLOR, linewidth = 1) +
    annotate("text", x = median(all_shap$MeanAbsSHAP), y = Inf,
             label = sprintf("Median: %.4f", median(all_shap$MeanAbsSHAP)),
             vjust = 2, hjust = -0.1, color = ACCENT_COLOR, fontface = "bold") +
    scale_x_log10(labels = label_scientific()) +
    labs(
      title = "Distribution of Feature Importance (Mean |SHAP|)",
      subtitle = sprintf("%s | Log-scaled x-axis | %d features x %d genes", SAMPLE_NAME, 
                         n_distinct(all_shap$Feature), n_genes),
      x = "Mean |SHAP Value| (log scale)",
      y = "Count"
    ) +
    theme_publication(base_size = 12)
  
  save_publication_plot(p3, file.path(output_dir, "03_shap_distribution"),
                        width = 10, height = 7)
  
  # --- Plot 4: Boxplot of SHAP by Feature Rank ---
  top10_features <- global_importance$Feature[1:min(10, nrow(global_importance))]
  
  p4_data <- all_shap %>%
    filter(Feature %in% top10_features) %>%
    mutate(
      Feature = factor(Feature, levels = rev(top10_features)),
      Feature_Short = truncate_name(as.character(Feature), 30)
    )
  
  # Get median values for labels
  p4_medians <- p4_data %>%
    group_by(Feature) %>%
    summarise(median_val = median(MeanAbsSHAP, na.rm = TRUE), .groups = "drop")
  
  p4 <- ggplot(p4_data, aes(x = Feature, y = MeanAbsSHAP)) +
    geom_boxplot(fill = SECONDARY_COLOR, alpha = 0.7, color = "black", 
                 outlier.size = 0.8, outlier.alpha = 0.5) +
    geom_point(data = p4_medians, aes(x = Feature, y = median_val),
               color = ACCENT_COLOR, size = 2) +
    coord_flip() +
    labs(
      title = "SHAP Distribution Across Genes (Top 10 Features)",
      subtitle = sprintf("%s | Variability of feature importance | RNA-only", SAMPLE_NAME),
      x = NULL,
      y = "Mean |SHAP Value|",
      caption = "Each box shows the distribution of SHAP values across all target genes"
    ) +
    theme_publication(base_size = 12) +
    theme(axis.text.y = element_text(size = 10))
  
  save_publication_plot(p4, file.path(output_dir, "04_shap_boxplot_top_features"),
                        width = 12, height = 8)
  
  # --- Plot 5: Importance vs Consistency Scatter ---
  p5_data <- global_importance %>%
    slice_head(n = 50) %>%
    mutate(Feature_Short = truncate_name(Feature, 20))
  
  p5 <- ggplot(p5_data, aes(x = N_Genes, y = Mean_SHAP)) +
    geom_point(size = 3, alpha = 0.7, color = PRIMARY_COLOR) +
    geom_text(aes(label = Feature_Short), hjust = -0.1, vjust = 0.5, size = 2.5, 
              check_overlap = TRUE, color = "gray30") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.4))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    labs(
      title = "Feature Importance vs Consistency",
      subtitle = sprintf("%s | Top 50 features | More genes = more consistent predictor", SAMPLE_NAME),
      x = "Number of Genes Where Feature Has SHAP Value",
      y = "Mean |SHAP Value|",
      caption = "Features in upper-right are both important and consistently predictive"
    ) +
    theme_publication(base_size = 12)
  
  save_publication_plot(p5, file.path(output_dir, "05_importance_vs_consistency"),
                        width = 12, height = 9)
  
  cat("  ✓ Saved basic plots (1-5)\n")
  
  # ============================================================================
  # GENE-FEATURE HEATMAP
  # ============================================================================
  
  cat("\n=== Creating gene-feature heatmap ===\n")
  
  top_features <- global_importance %>% slice_head(n = 30) %>% pull(Feature)
  
  gene_importance <- all_shap %>%
    group_by(Gene) %>%
    summarise(MeanImportance = mean(MeanAbsSHAP, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(MeanImportance))
  
  n_genes_for_heatmap <- min(50, n_genes)
  top_genes <- gene_importance %>% slice_head(n = n_genes_for_heatmap) %>% pull(Gene)
  
  heatmap_data <- all_shap %>%
    filter(Feature %in% top_features, Gene %in% top_genes) %>%
    select(Gene, Feature, MeanAbsSHAP) %>%
    pivot_wider(names_from = Feature, values_from = MeanAbsSHAP, values_fill = 0)
  
  heatmap_long <- heatmap_data %>%
    pivot_longer(-Gene, names_to = "Feature", values_to = "SHAP") %>%
    mutate(
      Feature = factor(Feature, levels = top_features),
      Gene = factor(Gene, levels = rev(top_genes)),
      Feature_Short = truncate_name(as.character(Feature), 25)
    )
  
  heatmap_title <- sprintf("SHAP Values Heatmap (Top %d Genes x Top 30 Features)", n_genes_for_heatmap)
  
  p_heat <- ggplot(heatmap_long, aes(x = Feature, y = Gene, fill = SHAP)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_viridis_c(option = "plasma", name = "|SHAP|",
                         labels = label_scientific()) +
    labs(
      title = heatmap_title,
      subtitle = sprintf("%s | RNA-only (HVG features)", SAMPLE_NAME),
      x = NULL, 
      y = NULL,
      caption = "Brighter colors indicate higher feature importance for that gene"
    ) +
    theme_publication(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
  save_publication_plot(p_heat, file.path(output_dir, "06_shap_heatmap"),
                        width = 16, height = max(8, n_genes_for_heatmap * 0.2))
  
  cat(sprintf("  ✓ Saved heatmap (%d genes x 30 features)\n", n_genes_for_heatmap))
  
  # ============================================================================
  # ADDITIONAL ANALYSIS PLOTS
  # ============================================================================
  
  cat("\n=== Creating additional analysis plots ===\n")
  
  # --- Plot 7: Cumulative importance ---
  cumulative_data <- global_importance %>%
    mutate(
      Cumulative_SHAP = cumsum(Mean_SHAP),
      Cumulative_Pct = 100 * Cumulative_SHAP / sum(Mean_SHAP)
    )
  
  # Find key thresholds
  n_50 <- sum(cumulative_data$Cumulative_Pct <= 50) + 1
  n_80 <- sum(cumulative_data$Cumulative_Pct <= 80) + 1
  n_95 <- sum(cumulative_data$Cumulative_Pct <= 95) + 1
  
  p7 <- ggplot(cumulative_data, aes(x = Rank, y = Cumulative_Pct)) +
    geom_line(color = PRIMARY_COLOR, linewidth = 1.2) +
    geom_hline(yintercept = c(50, 80, 95), linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_point(data = cumulative_data %>% filter(Rank %in% c(n_50, n_80, n_95)),
               color = ACCENT_COLOR, size = 3) +
    annotate("text", x = c(n_50, n_80, n_95), y = c(50, 80, 95) + 3,
             label = c(sprintf("%d features", n_50), sprintf("%d features", n_80), sprintf("%d features", n_95)),
             hjust = 0, size = 3.5, color = ACCENT_COLOR, fontface = "bold") +
    scale_x_log10() +
    scale_y_continuous(breaks = seq(0, 100, 20)) +
    labs(
      title = "Cumulative Feature Importance",
      subtitle = sprintf("%s | How many features explain X%% of total importance", SAMPLE_NAME),
      x = "Number of Features (log scale)",
      y = "Cumulative % of Total SHAP",
      caption = sprintf("50%%: %d features | 80%%: %d features | 95%%: %d features", n_50, n_80, n_95)
    ) +
    theme_publication(base_size = 12)
  
  save_publication_plot(p7, file.path(output_dir, "07_cumulative_importance"),
                        width = 10, height = 7)
  
  # --- Plot 8: Per-gene variability ---
  gene_variability <- all_shap %>%
    group_by(Gene) %>%
    summarise(
      Mean_SHAP = mean(MeanAbsSHAP, na.rm = TRUE),
      SD_SHAP = sd(MeanAbsSHAP, na.rm = TRUE),
      Max_SHAP = max(MeanAbsSHAP, na.rm = TRUE),
      Top_Feature = Feature[which.max(MeanAbsSHAP)],
      .groups = "drop"
    ) %>%
    mutate(CV = SD_SHAP / Mean_SHAP)
  
  p8 <- ggplot(gene_variability, aes(x = Mean_SHAP, y = Max_SHAP)) +
    geom_point(alpha = 0.6, color = PRIMARY_COLOR, size = 2.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_smooth(method = "lm", se = TRUE, color = ACCENT_COLOR, alpha = 0.2, linewidth = 0.8) +
    labs(
      title = "Gene-Level SHAP Summary",
      subtitle = sprintf("%s | Mean vs Max feature importance per gene", SAMPLE_NAME),
      x = "Mean |SHAP| (across all features)",
      y = "Max |SHAP| (top feature)",
      caption = "Points above diagonal: genes with one dominant predictive feature"
    ) +
    theme_publication(base_size = 12)
  
  save_publication_plot(p8, file.path(output_dir, "08_gene_level_summary"),
                        width = 9, height = 8)
  
  write_csv(gene_variability, file.path(output_dir, "gene_level_shap_summary.csv"))
  
  cat("  ✓ Saved additional plots (7-8)\n")
  
  # ============================================================================
  # COMBINED PUBLICATION PLOT
  # ============================================================================
  
  cat("\n=== Creating combined publication figure ===\n")
  
  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "SHAP Feature Importance Analysis (scRNA-only)",
      subtitle = sprintf("%s | %d genes | %d HVG features", SAMPLE_NAME, n_genes, n_distinct(all_shap$Feature)),
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40")
      )
    )
  
  save_publication_plot(combined, file.path(output_dir, "combined_shap_analysis"),
                        width = 20, height = 16)
  
  cat("  ✓ Saved combined publication figure\n")
  
  # ============================================================================
  # SUMMARY STATISTICS
  # ============================================================================
  
  cat("\n=== Summary Statistics ===\n")
  
  summary_stats <- list(
    sample_name = SAMPLE_NAME,
    gene_set = gene_set_name,
    n_genes = n_genes,
    n_features = n_distinct(all_shap$Feature),
    mean_shap_overall = mean(all_shap$MeanAbsSHAP, na.rm = TRUE),
    median_shap_overall = median(all_shap$MeanAbsSHAP, na.rm = TRUE),
    top_feature = global_importance$Feature[1],
    top_feature_shap = global_importance$Mean_SHAP[1],
    top_feature_n_genes = global_importance$N_Genes[1],
    features_for_50pct = n_50,
    features_for_80pct = n_80,
    features_for_95pct = n_95
  )
  
  summary_df <- as.data.frame(summary_stats)
  write_csv(summary_df, file.path(output_dir, "summary_statistics.csv"))
  
  cat(sprintf("  Total genes: %d\n", summary_stats$n_genes))
  cat(sprintf("  Total features: %d\n", summary_stats$n_features))
  cat(sprintf("  Overall mean SHAP: %.4f\n", summary_stats$mean_shap_overall))
  cat(sprintf("  Top feature: %s (SHAP: %.4f, in %d genes)\n", 
              summary_stats$top_feature, summary_stats$top_feature_shap, 
              summary_stats$top_feature_n_genes))
  cat(sprintf("  Features for 50%% importance: %d\n", summary_stats$features_for_50pct))
  cat(sprintf("  Features for 80%% importance: %d\n", summary_stats$features_for_80pct))
  cat(sprintf("  Features for 95%% importance: %d\n", summary_stats$features_for_95pct))
  
  cat(sprintf("\n>>> Completed %s gene set <<<\n", gene_set_name))
}

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("STEP 082 COMPLETE - AGGREGATED SHAP ANALYSIS (scRNA-only)\n")
cat("=", rep("=", 70), "\n\n", sep = "")

cat(sprintf("Total runtime: %.2f minutes\n\n", as.numeric(runtime)))

cat("Outputs per gene set:\n")
cat("  - all_genes_shap_combined.csv (all SHAP values)\n")
cat("  - global_feature_importance.csv\n")
cat("  - feature_frequency_in_top20.csv\n")
cat("  - gene_level_shap_summary.csv\n")
cat("  - summary_statistics.csv\n")

cat("\nPublication-ready plots (PNG 600dpi + PDF vector):\n")
cat("  01 - Global feature importance (top 30)\n")
cat("  02 - Feature frequency in top-20\n")
cat("  03 - SHAP distribution histogram\n")
cat("  04 - Boxplot of top features across genes\n")
cat("  05 - Importance vs consistency scatter\n")
cat("  06 - Gene-feature SHAP heatmap\n")
cat("  07 - Cumulative importance curve\n")
cat("  08 - Gene-level summary scatter\n")
cat("  combined_shap_analysis (multi-panel figure)\n")

cat("\nPlot features:\n")
cat("  - Colorblind-friendly palette (Wong 2011)\n")
cat("  - High resolution PNG (600 dpi)\n")
cat("  - Vector PDF for publication\n")
cat("  - Consistent publication theme\n")

cat("\n", paste(rep("=", 70), collapse = ""), "\n\n")
