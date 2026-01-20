# Target Gene Lists

## Overview

This directory contains pre-computed target gene lists used in the published analysis. Target genes define **which genes' expression levels to predict** (the response variable Y). The same target genes are used across all three feature modalities to enable direct performance comparison:

- **combined/**: Predicts target gene expression using RNA + ATAC features
- **rna_only/**: Predicts target gene expression using RNA features only
- **atac_only/**: Predicts target gene expression using ATAC features only

## Directory Structure

```
data/target_genes/
├── E7.5_rep1/                          # Mouse ESC dataset (Replicate 1)
│   ├── target_genes_hvg_100.txt        # 100 highly variable genes
│   └── target_genes_random_100.txt     # 100 randomly selected non-HVG genes
├── E7.5_rep2/                          # Mouse ESC dataset (Replicate 2)
│   ├── target_genes_hvg_100.txt
│   └── target_genes_random_100.txt
└── T_Cells/                            # Human PBMC T cells dataset
    ├── target_genes_hvg_100.txt
    └── target_genes_random_100.txt
```

## File Descriptions

- **target_genes_hvg_100.txt**: Top 100 highly variable genes (HVGs) identified from training data. These genes exhibit high cell-to-cell expression variability and are typically of biological interest.

- **target_genes_random_100.txt**: 100 randomly selected non-HVG genes that meet minimum expression thresholds. These serve as a control group for comparison with HVG predictions.

## Reproducing Published Results

> ⚠️ **IMPORTANT**: To reproduce published results exactly, you **MUST** use these pre-computed files. Auto-selection (Step 2b) may yield different genes even with the same random seed due to minor preprocessing variations.

Using pre-computed target gene files makes **Step 2b (`02b_select_target_genes.R`) optional**.

To exactly reproduce the published analysis, specify the target gene files in your `config.R`:

```r
# Example configuration for E7.5_rep1 dataset
HVG_GENE_FILE <- "data/target_genes/E7.5_rep1/target_genes_hvg_100.txt"
RANDOM_GENE_FILE <- "data/target_genes/E7.5_rep1/target_genes_random_100.txt"

# Example configuration for T_Cells dataset
HVG_GENE_FILE <- "data/target_genes/T_Cells/target_genes_hvg_100.txt"
RANDOM_GENE_FILE <- "data/target_genes/T_Cells/target_genes_random_100.txt"
```

## Selecting Target Genes for NEW Datasets

### Option 1: Automated Selection (NEW datasets only)

> **Note**: This option is for analyzing datasets NOT included in the publication.

Execute `02b_select_target_genes.R` to identify target genes from your own data:

```bash
Rscript combined/R/02b_select_target_genes.R
```

This script:
1. Identifies highly variable genes from the training split only (preventing data leakage)
2. Selects the top N HVGs as target genes
3. Randomly samples N non-HVG genes meeting minimum expression criteria

### Option 2: Custom Gene List

Create a plain text file with one gene symbol per line:

```
Gene1
Gene2
Gene3
...
```

Then specify the path in `config.R`:

```r
TARGET_GENE_FILE <- "/path/to/your/custom_genes.txt"
```

## Important Notes

- Gene symbols must match the feature names in your Seurat object exactly
- For mouse data: Use standard gene symbols (e.g., `Nanog`, `Sox2`, `Pou5f1`)
- For human data: Use standard gene symbols (e.g., `NANOG`, `SOX2`, `POU5F1`)
- Identical target genes are used by all three pipelines (combined, rna_only, atac_only) to enable fair comparison
