# Publication Figure Scripts

This directory contains R scripts for generating cross-dataset comparison figures presented in the publication.

## Prerequisites

These scripts require completed analysis results from all three datasets:
- **E7.5_rep1**: Mouse embryonic stem cells (replicate 1)
- **E7.5_rep2**: Mouse embryonic stem cells (replicate 2)
- **T_Cells**: Human peripheral blood mononuclear cell T cells

## Script Descriptions

| Script | Description |
|--------|-------------|
| `Combine_All_Integrated_Analysis_Results.R` | Aggregates model performance results across all datasets into a unified analysis |
| `Multi_Dataset_Combined_Performance_Plots_HVG_SET.R` | Generates performance comparison plots for highly variable gene (HVG) target genes |
| `Multi_Dataset_Combined_Performance_Plots_Non-HVGs.R` | Generates performance comparison plots for randomly selected non-HVG target genes |
| `Multi_Dataset_Combined_Performance_Plots_Single_Modality.R` | Generate Perfromance plots for single modality pipeline (set path to RNA-only or ATAC-only results folder) across datasets |
| `Multi_Dataset_Combined_Alternative_Approach.R` | Evaluates alternative metacell construction approaches (WNN, scVI+PeakVI, MultiVI) |
| `RNA_Vs_Combined_Integrated_Analysis.R` | Comparative analysis of RNA-only versus combined (RNA+ATAC) feature modalities |
| `Normalized_Feature_Counts.R` | Analyzes and visualizes feature count distributions across datasets and modalities |
| `Scatter_Plot_Top_Feature_Importance_for_representative_gene.R` | SHAP-based feature importance visualization for representative target genes |
| `Extract_Gene_List.R` | Utility script for extracting and formatting gene lists from results |

## Usage

Before executing any script, update the input paths at the beginning of each file to point to your results directories:

```r
# Configuration section (modify in each script)
BASE_RESULTS_DIR <- "~/output/scMultiPreDICT_results"
DATASETS <- c("E7.5_rep1", "E7.5_rep2", "T_Cells")
OUTPUT_FIGURES_DIR <- "~/output/paper_figures"
```

Execute scripts individually:

```bash
cd paper_figures
Rscript Multi_Dataset_Combined_Performance_Plots_HVG_SET.R
```

## Output

Generated figures are saved in publication-ready formats (PDF, PNG) suitable for journal submission.
