# Preprocessing Scripts

This directory contains preprocessing scripts for specific datasets that require additional steps before running the main scMultiPreDICT pipeline.

## T_Cells Dataset (Human PBMC)

The T_Cells dataset was derived from the 10x Genomics PBMC 10k multiome dataset. Unlike the mouse embryo datasets (E7.5_rep1, E7.5_rep2) which use whole-tissue data, the T_Cells dataset required:

1. **Full PBMC preprocessing** - QC filtering on both RNA and ATAC modalities
2. **Cell type annotation** - Using SingleR with multiple reference databases
3. **T cell subset extraction** - Selecting cells annotated as T cells
4. **Peak re-calling** - Calling peaks specifically on T cell subset using MACS2

### Script: `00_pbmc_preprocessing_and_tcell_subset.R`

This script performs the complete preprocessing workflow:

```
Full PBMC Multiome Data
         ↓
┌─────────────────────────────────────┐
│  1. Load RNA + ATAC counts          │
│  2. Create Seurat object            │
│  3. QC filtering:                   │
│     - nFeature_RNA: 500-3200        │
│     - percent.mt < 15%              │
│     - nCount_RNA: 500-9000          │
│     - nCount_ATAC: 4000-90000       │
│     - nucleosome_signal < 1.2       │
│     - TSS.enrichment > 4            │
│  4. Normalization + PCA + UMAP      │
│  5. Clustering                      │
│  6. Cell type annotation (SingleR)  │
│     - HumanPrimaryCellAtlas         │
│     - DatabaseImmuneCellExpression  │
│     - MonacoImmuneData              │
└─────────────────────────────────────┘
         ↓
┌─────────────────────────────────────┐
│  7. Subset T cells                  │
│  8. Re-call peaks (MACS2)           │
│  9. Filter low-expression genes     │
│  10. Save T cell Seurat object      │
└─────────────────────────────────────┘
         ↓
   T_Cells Seurat Object
   (Input to Step 02 (Data Splitting))
```

### Dependencies

Additional R packages required:
- `SingleR` - Automated cell type annotation
- `celldex` - Reference datasets for SingleR
- `SingleCellExperiment` - Data structure for SingleR
- `EnsDb.Hsapiens.v86` - Human gene annotations
- `BSgenome.Hsapiens.UCSC.hg38` - Human genome sequence

External tools:
- `MACS2` - Peak calling software

### Usage

1. Update file paths at the top of the script
2. Run the preprocessing script:
   ```bash
   Rscript preprocessing/00_pbmc_preprocessing_and_tcell_subset.R
   ```
3. Use the output Seurat object as input to Step 02 (Data Splitting), skipping Step 01 (QC already done)

### Note

This script is specific to the PBMC T_Cells dataset. For mouse embryo datasets (E7.5_rep1, E7.5_rep2), start directly with Step 01 (Quality Control) of the main pipeline.
