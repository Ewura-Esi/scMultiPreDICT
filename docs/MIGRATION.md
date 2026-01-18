# Configuration Migration Guide

This document provides a reference for mapping between alternative configuration styles. The scMultiPreDICT pipeline uses flat variable naming for simplicity and consistency with the original validated implementation.

## Variable Name Mappings

### Input/Output Paths

| Alternative (Nested List) | Current (Flat Variable) |
|---------------------------|-------------------------|
| `DATA_PATHS$matrix_dir` | `INPUT_MTX`, `INPUT_FEATURES`, `INPUT_BARCODES` |
| `DATA_PATHS$fragments` | `INPUT_FRAGMENTS` |
| `OUTPUT_PATHS$qc_dir` | `OUTPUT_SEURAT_DIR` |
| `OUTPUT_PATHS$splits_dir` | `OUTPUT_SPLITS_DIR` |
| `OUTPUT_PATHS$metacells_dir` | `OUTPUT_METACELLS_DIR` |
| `OUTPUT_PATHS$features_dir` | `OUTPUT_FEATURES_DIR` |
| `OUTPUT_PATHS$models_dir` | `OUTPUT_MODELS_LINEAR_DIR` |
| `OUTPUT_PATHS$nn_models_dir` | `OUTPUT_MODELS_NN_DIR` |
| `OUTPUT_PATHS$plots_dir` | `OUTPUT_FIGURES_DIR` |

### Quality Control Parameters

| Alternative (Nested List) | Current (Flat Variable) |
|---------------------------|-------------------------|
| `QC_PARAMS$min_features_rna` | `QC_MIN_FEATURES_RNA` |
| `QC_PARAMS$max_features_rna` | `QC_MAX_FEATURES_RNA` |
| `QC_PARAMS$min_count_rna` | `QC_MIN_COUNT_RNA` |
| `QC_PARAMS$max_count_rna` | `QC_MAX_COUNT_RNA` |
| `QC_PARAMS$min_features_atac` | `QC_MIN_FEATURES_ATAC` |
| `QC_PARAMS$max_features_atac` | `QC_MAX_FEATURES_ATAC` |
| `QC_PARAMS$min_count_atac` | `QC_MIN_COUNT_ATAC` |
| `QC_PARAMS$max_count_atac` | `QC_MAX_COUNT_ATAC` |
| `QC_PARAMS$percent_mt_cutoff` | `QC_PERCENT_MT_CUTOFF` |

### Data Splitting Parameters

| Alternative (Nested List) | Current (Flat Variable) |
|---------------------------|-------------------------|
| `SPLIT_PARAMS$train_fraction` | `SPLIT_TRAIN` |
| `SPLIT_PARAMS$val_fraction` | `SPLIT_VAL` |
| `SPLIT_PARAMS$seed` | `SEED_SPLIT` |

### Metacell Creation Parameters

| Alternative (Nested List) | Current (Flat Variable) |
|---------------------------|-------------------------|
| `METACELL_PARAMS$k_neighbors` | `K_NEIGHBORS` |
| `METACELL_PARAMS$seed` | `SEED_METACELL` |
| `METACELL_PARAMS$approaches` | `METACELL_APPROACHES` |
| `METACELL_PARAMS$n_pcs` | `N_PCS` |
| `METACELL_PARAMS$n_lsi` | `N_LSI` |

### Feature Extraction Parameters

| Alternative (Nested List) | Current (Flat Variable) |
|---------------------------|-------------------------|
| `FEATURE_PARAMS$window_size` | `GENE_WINDOW_KB` |
| `FEATURE_PARAMS$atac_strategy` | `ATAC_STRATEGY` |
| `FEATURE_PARAMS$modalities` | `FEATURE_MODALITIES` |

### Model Training Parameters

| Alternative (Nested List) | Current (Flat Variable) |
|---------------------------|-------------------------|
| `MODEL_PARAMS$seed` | `SEED_MODEL` |
| `MODEL_PARAMS$models` | `MODELS_TO_TRAIN` |
| `MODEL_PARAMS$cv_folds` | `CV_FOLDS` |
| `MODEL_PARAMS$rf_n_trees` | `RF_N_TREES` |
| `MODEL_PARAMS$rf_min_node_size` | `RF_MIN_NODE_SIZE` |
| `MODEL_PARAMS$nn_hidden_units` | `NN_HIDDEN_UNITS` |
| `MODEL_PARAMS$nn_n_hidden_layers` | `NN_N_HIDDEN_LAYERS` |
| `MODEL_PARAMS$nn_dropout_rate` | `NN_DROPOUT_RATE` |
| `MODEL_PARAMS$nn_learning_rate` | `NN_LEARNING_RATE` |
| `MODEL_PARAMS$nn_batch_size` | `NN_BATCH_SIZE` |
| `MODEL_PARAMS$nn_epochs` | `NN_MAX_EPOCHS` |
| `MODEL_PARAMS$nn_patience` | `NN_EARLY_STOP_PATIENCE` |

## Configuration Style Comparison

**Nested list style:**
```r
DATA_PATHS <- list(
  matrix_dir = "/path/to/matrix",
  fragments = "/path/to/fragments.tsv.gz"
)

QC_PARAMS <- list(
  min_features_rna = 1000,
  max_features_rna = 10000
)
```

**Flat variable style (current implementation):**
```r
# Input paths
INPUT_MTX <- "/path/to/matrix/matrix.mtx.gz"
INPUT_FEATURES <- "/path/to/matrix/features.tsv.gz"
INPUT_BARCODES <- "/path/to/matrix/barcodes.tsv.gz"
INPUT_FRAGMENTS <- "/path/to/fragments.tsv.gz"

# QC parameters
QC_MIN_FEATURES_RNA <- 1000
QC_MAX_FEATURES_RNA <- 10000
```

## Design Rationale

The flat variable structure was selected for the following reasons:

1. **Simplicity**: Direct variable access (`QC_MIN_FEATURES_RNA`) versus nested access (`QC_PARAMS$min_features_rna`)
2. **Readability**: Variable names are self-documenting
3. **Consistency**: Matches the original validated implementation
4. **Familiarity**: Standard pattern for R configuration files
5. **Reduced typing**: Shorter variable references throughout scripts
