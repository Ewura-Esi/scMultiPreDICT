# ============================================================================
# PIPELINE CONFIGURATION TEMPLATE
# ============================================================================
# This is a TEMPLATE file. To use:
#   1. Copy this file to config.R:  cp config_template.R config.R
#   2. Edit config.R with YOUR dataset-specific paths and parameters
#   3. Run the pipeline: Rscript run_pipeline.R
#
# DO NOT edit this template directly - keep it as a reference.
# ============================================================================

# ============================================================================
# SECTION 1: SAMPLE IDENTIFICATION
# ============================================================================
# Sample name - used for file naming and output organization
SAMPLE_NAME <- "YOUR_SAMPLE_NAME"  # e.g., "E7.5_rep1", "T_Cells", "PBMC_10k"

# Project name (optional, for documentation)
PROJECT_NAME <- "YOUR_PROJECT_NAME"

# ============================================================================
# SECTION 2: INPUT DATA PATHS
# ============================================================================
# Paths to the raw data files (10X Genomics format)
# Use absolute paths or paths relative to your home directory (~/)

# RNA data paths (Matrix Market format)
INPUT_MTX <- "/path/to/your/matrix.mtx.gz"
INPUT_FEATURES <- "/path/to/your/features.tsv.gz"
INPUT_BARCODES <- "/path/to/your/barcodes.tsv.gz"

# ATAC fragments file (must be sorted and indexed with .tbi file)
INPUT_FRAGMENTS <- "/path/to/your/fragments.tsv.gz"

# ============================================================================
# SECTION 3: OUTPUT DIRECTORIES
# ============================================================================
# Base output directory - all outputs will be organized under this
BASE_OUTPUT_DIR <- "~/scMultiPreDICT_output/processed/"

# ============================================================================
# SECTION 3a: DIMENSIONALITY REDUCTION METHOD
# ============================================================================
# Choose which dimensionality reduction approach to use for metacell construction
# Options:
#   "pca_lsi"      - PCA for RNA, LSI for ATAC (default, fast, done in R)
#   "wnn"          - Weighted Nearest Neighbors from Seurat (done in R)
#   "multivi"      - Joint MultiVI autoencoder (requires Python)
#   "scvi_peakvi"  - Separate scVI + PeakVI autoencoders (requires Python)

DIM_REDUCTION_METHOD <- "pca_lsi"

# Short suffix for folder names (auto-derived from DIM_REDUCTION_METHOD)
DIMRED_METHOD_SUFFIX <- switch(
  DIM_REDUCTION_METHOD,
  "pca_lsi" = "pca_lsi",
  "wnn" = "wnn",
  "multivi" = "multivi",
  "scvi_peakvi" = "scvi_peakvi",
  DIM_REDUCTION_METHOD
)

# ----------------------------------------------------------------------------
# Derived output directories (automatically constructed)
# ----------------------------------------------------------------------------
# Steps 01-02: Common preprocessing (NOT method-specific)
OUTPUT_SEURAT_DIR <- file.path(BASE_OUTPUT_DIR, "seurat_obj", SAMPLE_NAME)
OUTPUT_SPLITS_DIR <- file.path(BASE_OUTPUT_DIR, "splits", SAMPLE_NAME)
OUTPUT_TARGET_GENES_DIR <- file.path(BASE_OUTPUT_DIR, "target_genes", SAMPLE_NAME)
OUTPUT_LATENT_DIR <- file.path(BASE_OUTPUT_DIR, "latent_spaces", SAMPLE_NAME)

# Steps 03+: Method-specific outputs
OUTPUT_METACELLS_DIR <- file.path(BASE_OUTPUT_DIR, "metacells", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)
OUTPUT_FEATURES_DIR <- file.path(BASE_OUTPUT_DIR, "features_extracted", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)

# ============================================================================
# SECTION 4: SPECIES/GENOME CONFIGURATION
# ============================================================================
# Supported options: "mouse" or "human"
SPECIES <- "mouse"  # Change to "human" for human data

# Species-specific settings (auto-configured)
if (SPECIES == "mouse") {
  GENOME <- "mm10"
  ENSDB_PACKAGE <- "EnsDb.Mmusculus.v79"
  BSGENOME_PACKAGE <- "BSgenome.Mmusculus.UCSC.mm10"
  MT_PATTERN <- "^mt-"
}

if (SPECIES == "human") {
  GENOME <- "hg38"
  ENSDB_PACKAGE <- "EnsDb.Hsapiens.v86"
  BSGENOME_PACKAGE <- "BSgenome.Hsapiens.UCSC.hg38"
  MT_PATTERN <- "^MT-"
}

# ============================================================================
# SECTION 5: QUALITY CONTROL PARAMETERS
# ============================================================================
# Adjust these thresholds based on your data quality

# RNA QC thresholds
QC_MIN_FEATURES_RNA <- 1000      # Minimum genes detected per cell
QC_MAX_FEATURES_RNA <- 6000      # Maximum genes detected per cell (doublet filter)
QC_MAX_PERCENT_MT <- 35          # Maximum mitochondrial percentage

# ATAC QC thresholds
QC_MIN_COUNT_ATAC <- 1000        # Minimum ATAC fragments per cell
QC_MAX_COUNT_ATAC <- 100000      # Maximum ATAC fragments per cell
QC_MIN_COUNT_RNA <- 1000         # Minimum RNA counts per cell
QC_MAX_COUNT_RNA <- 40000        # Maximum RNA counts per cell
QC_MAX_NUCLEOSOME_SIGNAL <- 1.7  # Maximum nucleosome signal
QC_MIN_TSS_ENRICHMENT <- 3       # Minimum TSS enrichment score

# Gene expression filter: Keep genes expressed in at least this fraction of cells
QC_MIN_GENE_FRACTION <- 0.10

# ============================================================================
# SECTION 6: DATA SPLITTING PARAMETERS
# ============================================================================
# Train/Validation/Test split proportions (must sum to 1.0)
SPLIT_TRAIN <- 0.70
SPLIT_VAL <- 0.20
SPLIT_TEST <- 0.10

# Random seed for reproducibility
SEED_SPLIT <- 42

# ============================================================================
# SECTION 7: METACELL/KNN SMOOTHING PARAMETERS
# ============================================================================
# Dimensionality reduction parameters
PCA_DIMS <- 30                   # Number of PCA dimensions for RNA
LSI_DIMS <- 30                   # Number of LSI dimensions for ATAC
N_VARIABLE_FEATURES <- 3000      # Number of highly variable genes for PCA

# k-NN smoothing parameters
K_NEIGHBORS <- 20                # Number of neighbors for k-NN smoothing

# Random seed for reproducibility
SEED_METACELL <- 2025

# ============================================================================
# SECTION 8: TARGET GENE CONFIGURATION
# ============================================================================
# Target genes define WHICH genes' expression to predict (the response variable Y).
#
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  OPTION 1: Use pre-computed files (REQUIRED for reproducing results)    ║
# ║  OPTION 2: Auto-select (ONLY for analyzing NEW datasets)                ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# ---- OPTION 1: Pre-computed Target Gene Files ----
# For reproducing published results, use the provided files:
#   E7.5_rep1: "data/target_genes/E7.5_rep1/target_genes_hvg_100.txt"
#   E7.5_rep2: "data/target_genes/E7.5_rep2/target_genes_hvg_100.txt"
#   T_Cells:   "data/target_genes/T_Cells/target_genes_hvg_100.txt"
#
# Set to empty string "" to use auto-selection (Option 2)
HVG_GENE_FILE <- ""
RANDOM_GENE_FILE <- ""

# ---- OPTION 2: Auto-selection Parameters (for NEW datasets) ----
N_HVG_GENES <- 100               # Number of HVGs to use as target genes
N_RANDOM_TARGET_GENES <- 100     # Number of random non-HVG genes
TARGET_MIN_DETECTION <- 5        # Min % of cells expressing the gene
TARGET_MAX_DETECTION <- 95       # Max % of cells expressing the gene
TARGET_MIN_MEAN_EXPR <- 0.1      # Min mean expression (log-normalized)
SEED_TARGET_GENES <- 2025        # Random seed for target gene selection

# Path for auto-generated target gene list (created by Step 02b)
TARGET_GENE_FILE <- file.path(OUTPUT_TARGET_GENES_DIR, 
                               sprintf("target_genes_random_%d.txt", N_RANDOM_TARGET_GENES))

# ============================================================================
# SECTION 9: AUTOENCODER LATENT SPACE PATHS (for multivi, scvi_peakvi methods)
# ============================================================================
# These paths are auto-constructed based on OUTPUT_LATENT_DIR
MULTIVI_LATENT_PATH <- file.path(OUTPUT_LATENT_DIR, "latent_multivi_all.csv")
SCVI_LATENT_PATH <- file.path(OUTPUT_LATENT_DIR, "latent_scvi_rna_all.csv")
PEAKVI_LATENT_PATH <- file.path(OUTPUT_LATENT_DIR, "latent_peakvi_atac_all.csv")

# ============================================================================
# SECTION 10: FEATURE EXTRACTION PARAMETERS
# ============================================================================
# Genomic window around TSS for ATAC peak extraction
GENE_WINDOW_KB <- 250            # ±250kb window from TSS

# Minimum peaks required per gene
MIN_PEAKS_PER_GENE <- 1

# Number of HVG expression features to include
N_HVG_FEATURES <- 50

# ============================================================================
# SECTION 11: MODEL TRAINING PARAMETERS
# ============================================================================
# Which gene sets to train models on
MODEL_GENE_SET <- "both"         # Options: "HVG", "Random_genes", "both"

# Random seed for model training
SEED_MODEL <- 123

# --- Random Forest ---
RF_N_TREES <- 500

# --- Deep Learning (Step 06) ---
NN_HIDDEN_UNITS <- 256           # Units in first hidden layer
NN_N_HIDDEN_LAYERS <- 3          # Number of hidden layers (1, 2, or 3)
NN_DROPOUT_RATE <- 0.1           # Dropout rate for regularization
NN_LEARNING_RATE <- 0.001        # Adam optimizer learning rate
NN_BATCH_SIZE <- 256             # Training batch size
NN_MAX_EPOCHS <- 100             # Maximum training epochs
NN_EARLY_STOP_PATIENCE <- 10     # Early stopping patience

# Grid search options (set to FALSE to use fixed params above)
NN_USE_GRID_SEARCH <- TRUE
NN_GRID_UNITS <- c(256)
NN_GRID_DROPOUT <- c(0, 0.1)
NN_GRID_BATCH <- c(256, 512)

# Conda environment for TensorFlow/Keras
CONDA_ENV_NAME <- "your_conda_env"  # Change to your environment name

# ============================================================================
# SECTION 12: COMPUTATIONAL RESOURCES
# ============================================================================
N_CORES <- 4                     # Default cores for local runs
MAX_CORES_TRAINING <- 32         # Maximum cores for model training
USE_DISK_BACKED <- FALSE         # Use disk-backed matrices for large datasets

# Conda environment names
CONDA_ENV_R <- "your_r_env"      # R/Bioconductor environment
CONDA_ENV_PYTHON <- "your_py_env" # Python environment

# SLURM settings
SLURM_PARTITION <- "compute"
SLURM_LOG_DIR <- "~/logs"

# ============================================================================
# SECTION 13: RESULTS OUTPUT DIRECTORIES
# ============================================================================
BASE_RESULTS_DIR <- "~/scMultiPreDICT_output/results/"

# Model output directories
OUTPUT_MODELS_LINEAR_DIR <- file.path(BASE_RESULTS_DIR, "models/LINEAR_AND_TREE_BASED", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)

NN_ARCH_LABEL <- switch(
  as.character(NN_N_HIDDEN_LAYERS),
  "1" = "One_hidden_layer",
  "2" = "Two_hidden_layer",
  "3" = "Three_hidden_layer",
  paste0(NN_N_HIDDEN_LAYERS, "_hidden_layer")
)

OUTPUT_MODELS_NN_DIR <- file.path(BASE_RESULTS_DIR, "models/NEURAL_NETWORKS", SAMPLE_NAME, DIMRED_METHOD_SUFFIX, NN_ARCH_LABEL)
OUTPUT_SHAP_DIR <- file.path(BASE_RESULTS_DIR, "shap_analysis", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)
OUTPUT_FIGURES_DIR <- file.path(BASE_RESULTS_DIR, "figures", SAMPLE_NAME, DIMRED_METHOD_SUFFIX)

# ============================================================================
# SLURM RESOURCE REQUIREMENTS (adjust based on your dataset size)
# ============================================================================
SLURM_RESOURCES <- list(
  step_01 = list(cpus = 8, mem_gb = 64, time_hours = 4, job_name = "qc_preprocess"),
  step_02b = list(cpus = 4, mem_gb = 32, time_hours = 1, job_name = "target_genes"),
  step_020 = list(cpus = 4, mem_gb = 32, time_hours = 1, job_name = "data_split"),
  step_03b = list(cpus = 4, mem_gb = 64, time_hours = 1, job_name = "export_mudata"),
  step_030 = list(cpus = 32, mem_gb = 256, time_hours = 8, job_name = "metacell"),
  step_040 = list(cpus = 32, mem_gb = 256, time_hours = 12, job_name = "feature_extract"),
  step_050 = list(cpus = 32, mem_gb = 128, time_hours = 12, job_name = "model_linear"),
  step_060 = list(cpus = 16, mem_gb = 64, time_hours = 24, job_name = "model_nn")
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Create all output directories
create_output_directories <- function() {
  dirs <- c(
    OUTPUT_SEURAT_DIR,
    OUTPUT_SPLITS_DIR,
    OUTPUT_METACELLS_DIR,
    OUTPUT_FEATURES_DIR,
    OUTPUT_TARGET_GENES_DIR,
    OUTPUT_LATENT_DIR,
    OUTPUT_MODELS_LINEAR_DIR,
    OUTPUT_MODELS_NN_DIR,
    OUTPUT_SHAP_DIR,
    OUTPUT_FIGURES_DIR
  )
  
  for (dir in dirs) {
    dir <- path.expand(dir)
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat(sprintf("Created directory: %s\n", dir))
    }
  }
}

#' Load appropriate annotation packages
load_annotation_packages <- function() {
  if (SPECIES == "mouse") {
    suppressPackageStartupMessages({
      library(EnsDb.Mmusculus.v79)
      library(BSgenome.Mmusculus.UCSC.mm10)
    })
  } else if (SPECIES == "human") {
    suppressPackageStartupMessages({
      library(EnsDb.Hsapiens.v86)
      library(BSgenome.Hsapiens.UCSC.hg38)
    })
  }
}

#' Get EnsDb object based on species
get_ensdb <- function() {
  if (SPECIES == "mouse") {
    return(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
  } else if (SPECIES == "human") {
    return(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
  }
}

#' Print configuration summary
print_config <- function() {
  cat("\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("PIPELINE CONFIGURATION\n")
  cat("=", rep("=", 70), "\n\n", sep = "")
  
  cat("Sample:", SAMPLE_NAME, "\n")
  cat("Species:", SPECIES, "(", GENOME, ")\n")
  cat("Dimensionality Reduction:", DIM_REDUCTION_METHOD, "\n")
  cat("Base output:", BASE_OUTPUT_DIR, "\n\n")
  
  cat("QC Thresholds:\n")
  cat("  RNA features:", QC_MIN_FEATURES_RNA, "-", QC_MAX_FEATURES_RNA, "\n")
  cat("  MT percentage: <", QC_MAX_PERCENT_MT, "%\n")
  cat("  ATAC counts:", QC_MIN_COUNT_ATAC, "-", QC_MAX_COUNT_ATAC, "\n\n")
  
  cat("Data Splits:\n")
  cat("  Train:", SPLIT_TRAIN * 100, "%\n")
  cat("  Validation:", SPLIT_VAL * 100, "%\n")
  cat("  Test:", SPLIT_TEST * 100, "%\n\n")
  
  cat("Metacell Parameters:\n")
  cat("  PCA/LSI dims:", PCA_DIMS, "\n")
  cat("  k-NN neighbors:", K_NEIGHBORS, "\n\n")
}

#' Validate configuration
validate_config <- function() {
  errors <- c()
  
  # Check split proportions
  if (abs(SPLIT_TRAIN + SPLIT_VAL + SPLIT_TEST - 1.0) > 0.001) {
    errors <- c(errors, "Split proportions must sum to 1.0")
  }
  
  # Check species
  if (!SPECIES %in% c("mouse", "human")) {
    errors <- c(errors, "SPECIES must be 'mouse' or 'human'")
  }
  
  # Check input files exist
  if (!file.exists(path.expand(INPUT_MTX))) {
    errors <- c(errors, paste("INPUT_MTX not found:", INPUT_MTX))
  }
  if (!file.exists(path.expand(INPUT_FEATURES))) {
    errors <- c(errors, paste("INPUT_FEATURES not found:", INPUT_FEATURES))
  }
  if (!file.exists(path.expand(INPUT_BARCODES))) {
    errors <- c(errors, paste("INPUT_BARCODES not found:", INPUT_BARCODES))
  }
  if (!file.exists(path.expand(INPUT_FRAGMENTS))) {
    errors <- c(errors, paste("INPUT_FRAGMENTS not found:", INPUT_FRAGMENTS))
  }
  
  if (length(errors) > 0) {
    cat("Configuration errors:\n")
    for (err in errors) {
      cat("  ERROR:", err, "\n")
    }
    stop("Please fix configuration errors before proceeding.")
  }
  
  cat("Configuration validated successfully!\n")
}

# ============================================================================
# AUTO-LOAD MESSAGE
# ============================================================================
cat("\n[CONFIG] Loaded configuration for sample:", SAMPLE_NAME, "\n")
cat("[CONFIG] Species:", SPECIES, "| Genome:", GENOME, "\n")
cat("[CONFIG] Dimensionality Reduction:", DIM_REDUCTION_METHOD, "\n")
cat("[CONFIG] Use print_config() for full configuration\n\n")
