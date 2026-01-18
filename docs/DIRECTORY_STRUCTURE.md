# Output Directory Structure

## Overview

scMultiPreDICT organizes outputs into two separate base directories to maintain clear separation between processed data and analysis results:

1. **`BASE_OUTPUT_DIR`**: Intermediate preprocessing outputs (Seurat objects, metacells, features)
2. **`BASE_RESULTS_DIR`**: Final analysis results (model outputs, figures)

---

## Complete Directory Hierarchy

```
${BASE_OUTPUT_DIR}/                    # Preprocessing data outputs
│
├── seurat_obj/${SAMPLE}/
│   ├── ${SAMPLE}_seurat_multiome_processed.rds    # QC-filtered Seurat object
│   └── ${SAMPLE}_QC_plots.pdf                     # Quality control visualizations
│
├── splits/${SAMPLE}/
│   ├── ${SAMPLE}_seurat_obj_with_splits.rds       # Seurat object with split assignments
│   ├── ${SAMPLE}_target_genes_hvg.txt             # HVG target gene list
│   └── ${SAMPLE}_target_genes_random.txt          # Random target gene list
│
├── target_genes/${SAMPLE}/
│   └── target_genes_random_100.txt                # Selected random genes
│
├── latent_spaces/${SAMPLE}/                       # Autoencoder latent representations
│   ├── latent_multivi_all.csv                     # (if using MultiVI)
│   ├── latent_scvi_rna_all.csv                    # (if using scVI)
│   └── latent_peakvi_atac_all.csv                 # (if using PeakVI)
│
├── metacells/${SAMPLE}/${METHOD}/                 # e.g., pca_lsi, wnn, multivi
│   ├── smoothed_train.rds                         # Metacell-aggregated training data
│   ├── smoothed_val.rds                           # Metacell-aggregated validation data
│   └── smoothed_test.rds                          # Metacell-aggregated test data
│
└── features_extracted/${SAMPLE}/${METHOD}/
    ├── HVG/
    │   ├── gene_specific_features.rds             # All HVG gene features (single file)
    │   ├── gene_features_metadata.csv             # Feature extraction metadata
    │   └── feature_extraction_params.rds          # Extraction parameters
    └── Random_genes/
        ├── gene_specific_features.rds             # All random gene features (single file)
        ├── gene_features_metadata.csv
        └── feature_extraction_params.rds


${BASE_RESULTS_DIR}/                   # Analysis results
│
├── models/
│   │
│   ├── LINEAR_AND_TREE_BASED/${SAMPLE}/${METHOD}/
│   │   ├── HVG/
│   │   │   ├── ${GENE_NAME}/                      # Per-gene results folder
│   │   │   │   ├── metrics.csv                    # Performance metrics (R², MSE, MAE)
│   │   │   │   ├── predictions_train.csv
│   │   │   │   ├── predictions_val.csv
│   │   │   │   ├── predictions_test.csv
│   │   │   │   ├── coefficients_ols.csv
│   │   │   │   ├── coefficients_ridge.csv
│   │   │   │   ├── coefficients_lasso.csv
│   │   │   │   ├── coefficients_enet.csv
│   │   │   │   ├── coefficients_top20.csv
│   │   │   │   ├── rf_importance.csv
│   │   │   │   └── analysis_metadata.csv
│   │   │   └── ...
│   │   └── Random_genes/
│   │       └── ...
│   │
│   └── NEURAL_NETWORKS/${SAMPLE}/${METHOD}/${ARCHITECTURE}/
│       ├── HVG/
│       │   ├── ${GENE_NAME}/
│       │   │   ├── metrics.csv
│       │   │   ├── predictions_train.csv
│       │   │   ├── predictions_val.csv
│       │   │   ├── predictions_test.csv
│       │   │   └── nn_model.keras
│       │   └── ...
│       └── Random_genes/
│           └── ...
│
├── shap_analysis/${SAMPLE}/${METHOD}/
│   └── [SHAP analysis results]
│
└── figures/${SAMPLE}/${METHOD}/
    ├── model_comparison.pdf
    ├── feature_importance.pdf
    └── results_summary.csv
```

---

## Key Implementation Details

### Feature Storage Format

Features are stored in a single RDS file per gene set (not per-gene):

```r
# Load all features for HVG target genes
features <- readRDS("features_extracted/SAMPLE/pca_lsi/HVG/gene_specific_features.rds")

# Access specific gene data
names(features)  # Returns: c("Nanog", "Sox2", "Pou5f1", ...)

# Each gene contains split-specific data
features$Nanog$train$X      # Feature matrix (cells × features)
features$Nanog$train$y      # Target expression vector
features$Nanog$train$cells  # Cell identifiers
```

### Model Results Organization

Each target gene receives a dedicated results folder:

```r
# Read performance metrics
metrics <- read.csv("models/LINEAR_AND_TREE_BASED/SAMPLE/pca_lsi/HVG/Nanog/metrics.csv")

# Read predictions
predictions <- read.csv("models/LINEAR_AND_TREE_BASED/SAMPLE/pca_lsi/HVG/Nanog/predictions_test.csv")
```

### Method-Specific Paths

The `${METHOD}` placeholder corresponds to the dimensionality reduction approach:

| DIM_REDUCTION_METHOD | Path Suffix |
|---------------------|-------------|
| `linear` | `pca_lsi` |
| `wnn` | `wnn` |
| `scvi_peakvi` | `scvi_peakvi` |
| `multivi` | `multivi` |

---

## Configuration

Define base directories in `config.R`:

```r
# Base directories
BASE_OUTPUT_DIR <- "~/output/scMultiPreDICT_data"
BASE_RESULTS_DIR <- "~/output/scMultiPreDICT_results"

# Dimensionality reduction method
DIM_REDUCTION_METHOD <- "linear"

# Derived paths (computed automatically)
OUTPUT_FEATURES_DIR <- file.path(BASE_OUTPUT_DIR, "features_extracted", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)
OUTPUT_MODELS_LINEAR_DIR <- file.path(BASE_RESULTS_DIR, "models", "LINEAR_AND_TREE_BASED", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)
```

---

## Summary

| Aspect | Implementation |
|--------|---------------|
| Base directories | Two separate (data vs. results) |
| Feature storage | Single file per gene set containing all genes |
| Model outputs | Individual folder per gene with multiple CSV files |
| Path hierarchy | `BASE/TYPE/SAMPLE/METHOD/GENE_SET/GENE/` |
| Metacell naming | `smoothed_train.rds`, `smoothed_val.rds`, `smoothed_test.rds` |
