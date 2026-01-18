# ============================================================================
# scMultiPreDICT - Plotting Utilities
# ============================================================================
# Reusable functions for visualization and publication-ready figures
#
# Usage:
#   source("R/utils/plot_utils.R")
# ============================================================================

library(ggplot2)

#' Publication-ready ggplot theme
#'
#' @param base_size Base font size
#' @param base_family Base font family
#' @return ggplot2 theme
#' @export
theme_publication <- function(base_size = 12, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      # Panel
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      
      # Axes
      axis.text = element_text(color = "black", size = base_size * 0.9),
      axis.title = element_text(color = "black", size = base_size, face = "bold"),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      
      # Legend
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = base_size * 0.85),
      legend.title = element_text(size = base_size * 0.9, face = "bold"),
      
      # Title
      plot.title = element_text(size = base_size * 1.2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      
      # Strip (facets)
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(size = base_size * 0.9, face = "bold")
    )
}


#' Color palettes for models
#'
#' @return Named vector of colors
#' @export
get_model_colors <- function() {
  c(
    "OLS" = "#1f77b4",
    "Ridge" = "#ff7f0e",
    "Lasso" = "#2ca02c",
    "Elastic Net" = "#d62728",
    "Random Forest" = "#9467bd",
    "Neural Network" = "#8c564b",
    "DeepNN" = "#8c564b"
  )
}


#' Color palettes for modalities
#'
#' @return Named vector of colors
#' @export
get_modality_colors <- function() {
  c(
    "RNA-only" = "#1f77b4",
    "ATAC-only" = "#ff7f0e",
    "Combined" = "#2ca02c",
    "rna_only" = "#1f77b4",
    "atac_only" = "#ff7f0e",
    "combined" = "#2ca02c"
  )
}


#' Color palettes for metacell approaches
#'
#' @return Named vector of colors
#' @export
get_approach_colors <- function() {
  c(
    "PCA+LSI" = "#1f77b4",
    "WNN" = "#ff7f0e",
    "scVI+PeakVI" = "#2ca02c",
    "MultiVI" = "#d62728"
  )
}


#' Plot model performance comparison
#'
#' @param results_df Data frame with model performance results
#' @param metric Metric to plot ("r2", "pearson_r", "rmse", etc.)
#' @param group_by Grouping variable (default: "model")
#' @param title Plot title
#' @return ggplot object
#' @export
plot_model_comparison <- function(results_df, metric = "pearson_r",
                                   group_by = "model", title = NULL) {
  
  # Clean model names
  model_order <- c("OLS", "Ridge", "Lasso", "Elastic Net", 
                   "Random Forest", "Neural Network")
  
  results_df$model_clean <- dplyr::case_when(
    grepl("ols", tolower(results_df[[group_by]])) ~ "OLS",
    grepl("ridge", tolower(results_df[[group_by]])) ~ "Ridge",
    grepl("lasso", tolower(results_df[[group_by]])) ~ "Lasso",
    grepl("elastic", tolower(results_df[[group_by]])) ~ "Elastic Net",
    grepl("random|rf", tolower(results_df[[group_by]])) ~ "Random Forest",
    grepl("neural|nn|deep", tolower(results_df[[group_by]])) ~ "Neural Network",
    TRUE ~ results_df[[group_by]]
  )
  
  results_df$model_clean <- factor(results_df$model_clean, levels = model_order)
  
  # Get metric label
  metric_labels <- c(
    "r2" = expression(R^2),
    "pearson_r" = "Pearson r",
    "spearman_rho" = "Spearman ρ",
    "rmse" = "RMSE",
    "mae" = "MAE"
  )
  
  y_label <- if (metric %in% names(metric_labels)) metric_labels[metric] else metric
  
  p <- ggplot(results_df, aes(x = model_clean, y = .data[[metric]], 
                               fill = model_clean)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    scale_fill_manual(values = get_model_colors()) +
    labs(
      x = NULL,
      y = y_label,
      title = title
    ) +
    theme_publication() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}


#' Plot modality comparison
#'
#' @param results_df Data frame with results
#' @param metric Metric to plot
#' @param model_name Model to show (or NULL for all)
#' @param title Plot title
#' @return ggplot object
#' @export
plot_modality_comparison <- function(results_df, metric = "pearson_r",
                                      model_name = NULL, title = NULL) {
  
  if (!is.null(model_name)) {
    results_df <- results_df[grepl(model_name, results_df$model, ignore.case = TRUE), ]
  }
  
  # Clean modality names
  results_df$modality_clean <- dplyr::case_when(
    grepl("rna.*only|rna_only", tolower(results_df$modality)) ~ "RNA-only",
    grepl("atac.*only|atac_only", tolower(results_df$modality)) ~ "ATAC-only",
    grepl("combined|both", tolower(results_df$modality)) ~ "Combined",
    TRUE ~ results_df$modality
  )
  
  modality_order <- c("RNA-only", "ATAC-only", "Combined")
  results_df$modality_clean <- factor(results_df$modality_clean, levels = modality_order)
  
  p <- ggplot(results_df, aes(x = modality_clean, y = .data[[metric]],
                               fill = modality_clean)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    scale_fill_manual(values = get_modality_colors()) +
    labs(
      x = NULL,
      y = metric,
      title = title
    ) +
    theme_publication() +
    theme(legend.position = "none")
  
  return(p)
}


#' Plot prediction scatter with correlation
#'
#' @param y_true True values
#' @param y_pred Predicted values
#' @param gene_name Gene name for title
#' @param model_name Model name for title
#' @return ggplot object
#' @export
plot_prediction_scatter <- function(y_true, y_pred, gene_name = NULL,
                                     model_name = NULL) {
  
  df <- data.frame(
    true = y_true,
    pred = y_pred
  )
  
  # Calculate correlation
  r <- cor(y_true, y_pred, use = "complete.obs")
  r2 <- r^2
  
  # Build title
  title_parts <- c()
  if (!is.null(gene_name)) title_parts <- c(title_parts, gene_name)
  if (!is.null(model_name)) title_parts <- c(title_parts, model_name)
  title <- paste(title_parts, collapse = " - ")
  
  p <- ggplot(df, aes(x = true, y = pred)) +
    geom_point(alpha = 0.3, size = 0.5, color = "#1f77b4") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm", se = FALSE, color = "darkblue", linewidth = 0.8) +
    annotate("text", x = Inf, y = -Inf, 
             label = sprintf("r = %.3f\nR² = %.3f", r, r2),
             hjust = 1.1, vjust = -0.5, size = 4) +
    labs(
      x = "Observed Expression",
      y = "Predicted Expression",
      title = title
    ) +
    theme_publication() +
    coord_equal()
  
  return(p)
}


#' Plot SHAP feature importance
#'
#' @param shap_df Data frame with feature, importance columns
#' @param top_n Number of top features to show
#' @param title Plot title
#' @return ggplot object
#' @export
plot_shap_importance <- function(shap_df, top_n = 20, title = "Feature Importance") {
  
  # Get top features
  shap_df <- shap_df[order(-abs(shap_df$importance)), ]
  shap_df <- head(shap_df, top_n)
  
  # Order by importance
  shap_df$feature <- factor(shap_df$feature, 
                             levels = shap_df$feature[order(abs(shap_df$importance))])
  
  # Color by positive/negative
  shap_df$direction <- ifelse(shap_df$importance >= 0, "Positive", "Negative")
  
  p <- ggplot(shap_df, aes(x = importance, y = feature, fill = direction)) +
    geom_col(alpha = 0.8) +
    scale_fill_manual(values = c("Positive" = "#d62728", "Negative" = "#1f77b4")) +
    labs(
      x = "Mean |SHAP value|",
      y = NULL,
      title = title
    ) +
    theme_publication() +
    theme(legend.position = "bottom")
  
  return(p)
}


#' Plot heatmap of model x modality performance
#'
#' @param results_df Data frame with results
#' @param metric Metric for fill
#' @param title Plot title
#' @return ggplot object
#' @export
plot_performance_heatmap <- function(results_df, metric = "pearson_r",
                                      title = "Model Performance") {
  
  # Summarize by model and modality
  summary_df <- results_df %>%
    dplyr::group_by(model, modality) %>%
    dplyr::summarize(
      mean_metric = mean(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    )
  
  p <- ggplot(summary_df, aes(x = modality, y = model, fill = mean_metric)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", mean_metric)), 
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#1f77b4", mid = "white", high = "#d62728",
                         midpoint = mean(summary_df$mean_metric)) +
    labs(
      x = NULL,
      y = NULL,
      fill = metric,
      title = title
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  return(p)
}


#' Save publication-quality figure
#'
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi Resolution
#' @export
save_publication_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  
  # Determine format from extension
  ext <- tolower(tools::file_ext(filename))
  
  # Create directory if needed
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  
  if (ext == "pdf") {
    ggsave(filename, plot, width = width, height = height, device = "pdf")
  } else if (ext == "png") {
    ggsave(filename, plot, width = width, height = height, dpi = dpi, device = "png")
  } else if (ext == "svg") {
    ggsave(filename, plot, width = width, height = height, device = "svg")
  } else if (ext %in% c("tiff", "tif")) {
    ggsave(filename, plot, width = width, height = height, dpi = dpi, 
           device = "tiff", compression = "lzw")
  } else {
    ggsave(filename, plot, width = width, height = height, dpi = dpi)
  }
  
  cat(sprintf("Saved: %s\n", filename))
}


#' Create summary results table
#'
#' @param results_df Data frame with results
#' @param group_vars Variables to group by
#' @param metrics Metrics to summarize
#' @return Summary data frame
#' @export
create_summary_table <- function(results_df, 
                                  group_vars = c("model", "modality"),
                                  metrics = c("pearson_r", "r2", "rmse")) {
  
  summary_df <- results_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarize(
      n_genes = dplyr::n(),
      dplyr::across(
        dplyr::all_of(metrics),
        list(
          mean = ~mean(.x, na.rm = TRUE),
          sd = ~sd(.x, na.rm = TRUE),
          median = ~median(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
  
  return(summary_df)
}
