#!/usr/bin/env Rscript
# ============================================================================
# Multi-Dataset Combined Performance Plots (Publication Ready)
# ============================================================================

start_time <- Sys.time()

# ============================================================================
# 1. LOAD LIBRARIES
# ============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(viridis)
  library(ggrepel)
  library(scales)
  library(ggpubr)    # For statistical significance
  library(gghalves)  # For Raincloud plots (Half-violins)
  library(cowplot)   # For arranging complex grids
  library(Cairo)
  library(gridExtra)
})

# ============================================================================
# 2. CONFIGURATION - EDIT THIS SECTION
# ============================================================================
# Update BASE_RESULTS_DIR to point to your scMultiPreDICT output directory
# The expected folder structure is:
#   BASE_RESULTS_DIR/
#     LINEAR_AND_TREE_BASED/{sample}/{method}/{gene_set}/
#     NEURAL_NETWORKS/{sample}/{method}/Three_hidden_layer/{gene_set}/

BASE_RESULTS_DIR <- "~/scMultiPreDICT_output/results/models"

# Helper function to construct paths
make_paths <- function(base, sample, method, gene_set = "HVG") {
  list(
    linear_path = file.path(base, "LINEAR_AND_TREE_BASED", sample, method, gene_set),
    nn_path = file.path(base, "NEURAL_NETWORKS", sample, method, "Three_hidden_layer", "Three_hidden_layer", gene_set)
  )
}

# Define datasets and methods
SAMPLES <- c("E7.5_rep1", "E7.5_rep2", "T_Cells")
SAMPLE_LABELS <- c("E7.5_REP1", "E7.5_REP2", "T_Cells")
GENE_SET <- "HVG"  # or "Random_genes" for non-HVG analysis

# Build path configurations for each method
pca_lsi <- setNames(lapply(SAMPLES, function(s) make_paths(BASE_RESULTS_DIR, s, "pca_lsi", GENE_SET)), SAMPLE_LABELS)
wnn <- setNames(lapply(SAMPLES, function(s) make_paths(BASE_RESULTS_DIR, s, "wnn", GENE_SET)), SAMPLE_LABELS)
scvi_peakvi <- setNames(lapply(SAMPLES, function(s) make_paths(BASE_RESULTS_DIR, s, "scvi_peakvi", GENE_SET)), SAMPLE_LABELS)
multivi <- setNames(lapply(SAMPLES, function(s) make_paths(BASE_RESULTS_DIR, s, "multivi", GENE_SET)), SAMPLE_LABELS)

# Output directory for figures
output_dir <- "~/scMultiPreDICT_output/paper_figures/integrated_analysis/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================================
# 3. FIGURES THEME & PALETTES
# ============================================================================

# Dataset palette 
dataset_cols <- c("pca_lsi" = "#E69F00", "wnn" = "#56B4E9", "scvi_peakvi" = "#009E73", "multivi" = "#CC79A7")
# Ash color used for jitter points ('glittering')
ash_col <- "#3F3837"
model_cols   <- c("OLS" = "#999999", "Ridge" = "#0072B2", "Lasso" = "#009E73", 
                  "ElasticNet" = "#E69F00", "RandomForest" = "#D55E00", "DeepNN" = "#CC79A7")

theme_nature <- function(base_size = 12, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle  = element_text(size = base_size, hjust = 0.5, color = "gray35"),
      plot.caption = element_text(size = base_size - 2, color = "gray40", hjust = 1),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text  = element_text(size = base_size - 1, color = "black"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text  = element_text(size = base_size - 1),
      panel.grid.major= element_line(color = "gray92", linewidth = 0.3),
      panel.grid.minor= element_blank(),
      strip.background= element_rect(fill = "gray95", color = "black", linewidth = 0.4),
      strip.text = element_text(face = "bold"),
      plot.margin  = margin(12, 12, 12, 12)
    )
}
theme_set(theme_nature())

# ============================================================================
# 4. DATA LOADING FUNCTIONS
# ============================================================================

find_file <- function(base_path) {
  # Priorities: specific combined metrics -> generic summary -> recursive search
  patterns <- c("all_genes_combined_metrics.csv", "metrics_summary.csv", "combined_metrics.csv", "aggregated_metrics.csv")
  for (p in patterns) {
    f <- file.path(base_path, p)
    if (file.exists(f)) return(f)
  }
  # Recursive search if standard path fails
  fs <- list.files(base_path, pattern = "combined_metrics.csv", recursive = TRUE, full.names = TRUE)
  if (length(fs) > 0) return(fs[1])
  return(NULL)
}

load_data_safe <- function(dataset_name, config) {
  # Load Linear
  f_lin <- find_file(config$linear_path)
  d_lin <- if (!is.null(f_lin)) read_csv(f_lin, show_col_types = FALSE) else NULL
  
  # Load NN
  f_nn <- find_file(config$nn_path)
  d_nn <- if (!is.null(f_nn)) read_csv(f_nn, show_col_types = FALSE) else NULL
  
  # Combine
  if (is.null(d_lin) & is.null(d_nn)) return(NULL)
  combined <- bind_rows(d_lin, d_nn) %>%
    mutate(Dataset = dataset_name) %>%
    filter(Split == "Test") # Ensure we only use Test data
  
  # Clean Model Names
  combined <- combined %>%
    rename_with(str_to_title) %>% # Fix column case issues (model -> Model)
    # Ensure consistent metric column names (e.g. RMSE -> RMSE not Rmse)
    rename_with(~ ifelse(. == "Rmse", "RMSE", .)) %>%
    mutate(Model = case_when(
      str_detect(Model, "(?i)ols") ~ "OLS",
      str_detect(Model, "(?i)ridge") ~ "Ridge",
      str_detect(Model, "(?i)lasso") ~ "Lasso",
      str_detect(Model, "(?i)elastic") ~ "ElasticNet",
      str_detect(Model, "(?i)random|rf") ~ "RandomForest",
      str_detect(Model, "(?i)deep|nn") ~ "DeepNN",
      TRUE ~ Model
    )) %>%
    filter(Model %in% names(model_cols)) # Filter only known models
  
  return(combined)
}

# ============================================================================
# 5. MAIN EXECUTION
# ============================================================================

# PCA + LSI
cat("Loading Datasets...\n")
pca_lsi_data <- map_dfr(names(pca_lsi), ~load_data_safe(., pca_lsi[[.]]))

# Factor Ordering (Crucial for plotting order)
pca_lsi_data$Model <- factor(pca_lsi_data$Model, levels = c("OLS", "Lasso", "Ridge", "ElasticNet", "RandomForest", "DeepNN"))
pca_lsi_data$Dataset <- factor(pca_lsi_data$Dataset, levels = names(pca_lsi))

cat(sprintf("Loaded %d rows across %d datasets.\n", nrow(pca_lsi_data), n_distinct(pca_lsi_data$Dataset)))

# WNN
cat("Loading Datasets...\n")
wnn_data <- map_dfr(names(wnn), ~load_data_safe(., wnn[[.]]))

# Factor Ordering (Crucial for plotting order)
wnn_data$Model <- factor(wnn_data$Model, levels = c("OLS", "Lasso", "Ridge", "ElasticNet", "RandomForest", "DeepNN"))
wnn_data$Dataset <- factor(wnn_data$Dataset, levels = names(wnn))

cat(sprintf("Loaded %d rows across %d datasets.\n", nrow(wnn_data), n_distinct(wnn_data$Dataset)))

# SCVI + PEAKVI
cat("Loading Datasets...\n")
scvi_peakvi_data <- map_dfr(names(scvi_peakvi), ~load_data_safe(., scvi_peakvi[[.]]))

# Factor Ordering (Crucial for plotting order)
scvi_peakvi_data$Model <- factor(scvi_peakvi_data$Model, levels = c("OLS", "Lasso", "Ridge", "ElasticNet", "RandomForest", "DeepNN"))
scvi_peakvi_data$Dataset <- factor(scvi_peakvi_data$Dataset, levels = names(scvi_peakvi))

cat(sprintf("Loaded %d rows across %d datasets.\n", nrow(scvi_peakvi_data), n_distinct(scvi_peakvi_data$Dataset)))

# MULTIVI
cat("Loading Datasets...\n")
multivi_data <- map_dfr(names(multivi), ~load_data_safe(., multivi[[.]]))

# Factor Ordering (Crucial for plotting order)
multivi_data$Model <- factor(multivi_data$Model, levels = c("OLS", "Lasso", "Ridge", "ElasticNet", "RandomForest", "DeepNN"))
multivi_data$Dataset <- factor(multivi_data$Dataset, levels = names(multivi))

cat(sprintf("Loaded %d rows across %d datasets.\n", nrow(multivi_data), n_distinct(multivi_data$Dataset)))

# ============================================================================
# 6. Integrate all the dataset into a big dataframe
# ============================================================================
 # Add an Integration column to each dimensionality strategy

pca_lsi_data$Integrated_Analysis <- "PCA_LSI"
scvi_peakvi_data$Integrated_Analysis <- "SCVI_PEAKVI"
multivi_data$Integrated_Analysis <- "MULTIVI"
wnn_data$Integrated_Analysis <- "WNN"


# Combine all integration tables into one master table
df_all <- bind_rows(
  pca_lsi_data,
  wnn_data,
  scvi_peakvi_data,
  multivi_data
)

# Select columns needed for plotting
df_main <- df_all %>% 
  filter(
    Split == "Test",
    Geneset == "HVG"
  )

# Collapse across models
df_gene_level <- df_main %>% 
  group_by(Dataset, Integrated_Analysis, Gene) %>% 
  summarise(
    spearman = median(Spearman, na.rm = TRUE),
    r2 = median(R2, na.rm = TRUE),
    rmse = median(RMSE, na.rm = TRUE),
    .groups = "drop"
  )

# Keep Integration method ordered
df_gene_level$Integrated_Analysis <- factor(df_gene_level$Integrated_Analysis,
                                    levels = c("PCA_LSI", "WNN", "SCVI_PEAKVI","MULTIVI"))

# Pool across datasets
df_pooled <- df_gene_level %>% 
  select(Integrated_Analysis, spearman)

#df_heatmap <- df_gene_level %>% 
 # group_by(Dataset, Integrated_Analysis) %>% 
  # summarise(
  #  median_spearman = median(spearman, na.rm = TRUE),
  #  median_r2 = median(r2, na.rm = TRUE),
   # median_rmse = median(rmse, na.rm = TRUE),
  #  .groups = "drop"
  # )


# ============================================================================
# 6. GENERATE FIGURES
# ============================================================================

pA <- ggplot(df_gene_level, aes(x = Integrated_Analysis, y = spearman, fill = Integrated_Analysis)) + 
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.6) + 
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(shape = 21, fill = ash_col, color = "black", stroke = 0.08, width = 0.12, size = 0.6, alpha = 0.25) +
  scale_y_continuous(limits = c(-0.1, 1.0), breaks = seq(0, 1, 0.2)) +
  labs(title = "Pooled performance across integrated analysis strategies",x = NULL, y = "Spearman (gene-level median across models)") + 
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

pA

ggsave(file.path(output_dir, "Fig4_Integrated_Analysis_Comparison.pdf"), pA, device = cairo_pdf,
       width = 10, height = 5 + length(metrics_present) * 0.3, units = "in")
ggsave(file.path(output_dir, "Fig4_Integrated_Analysis_Comparison.png"), pA, width = 10, height = 5 + length(metrics_present) * 0.3, units = "in", dpi = 600)


pB <- ggplot(df_heatmap , aes(x = Integrated_Analysis, y = Dataset, fill = median_spearman)) + 
  geom_tile(color = "white", linewidth = 0.4) + 
  labs(x = NULL, y = NULL, fill = "Median\nSpearman") + 
  theme_minimal(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  geom_text(aes(label = round(median_spearman, 2)), size = 3)

pB


# Ensure metrics exist and compute simple summary table for heatmaps
metrics_present <- c("R2", "Spearman", "RMSE")


# Prepare long-form mean values per Metric/Dataset/Model
heatmap_dat <- df_main %>%
  pivot_longer(cols = all_of(metrics_present), names_to = "Metric", values_to = "Value") %>%
  group_by(Metric, Dataset, Integrated_Analysis) %>%
  summarise(Median = median(Value, na.rm = TRUE), .groups = "drop") %>%
  group_by(Metric) %>%
  mutate(Thresh = (max(Median, na.rm = TRUE) + min(Median, na.rm = TRUE)) / 2,
         TextColor = ifelse(is.na(Median), "black", ifelse(Median >= Thresh, "white", "black"))) %>%
  ungroup()

# Keep Integration method ordered
heatmap_dat$Integrated_Analysis <- factor(heatmap_dat$Integrated_Analysis,
                                            levels = c("PCA_LSI", "WNN", "SCVI_PEAKVI","MULTIVI"))

heatmap_dat$Dataset <- factor(heatmap_dat$Dataset)
heatmap_dat$Metric <- factor(heatmap_dat$Metric, 
                             levels = c("Spearman","R2","RMSE"))

p_heat_all <- ggplot(heatmap_dat, aes(x = Integrated_Analysis, y = Dataset, fill = Median)) +
  geom_tile(color = "white", height = 0.92) +
  geom_text(aes(label = sprintf("%.2f", Median), color = TextColor), fontface = "bold", size = 4.0, show.legend = FALSE) +
  facet_wrap(~Metric, nrow = 1, scales = "free_x") +
  scale_fill_viridis_c(option = "maco", direction = 1, na.value = "grey80", name = "Median") +
  scale_color_identity() +
  labs(title = "Heatmap summarizing pooled multimodal performance across metrics", x = NULL, y = NULL) +
  theme_nature(base_size = 12) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

p_heat_all

ggsave(file.path(output_dir, "Fig4_Heatmap_AllMetrics_Integrated_Analysis_Comparison.pdf"), p_heat_all, device = cairo_pdf,
       width = 8, height = 3 + length(metrics_present) * 0.2, units = "in")
ggsave(file.path(output_dir, "Fig4_Heatmap_AllMetrics_Integrated_Analysis_Comparison.png"), p_heat_all, width = 10, height = 5 + length(metrics_present) * 0.3, units = "in", dpi = 600)
