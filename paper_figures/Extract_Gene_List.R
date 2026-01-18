# ==============================================================================
# SETUP
# ==============================================================================
# Install if you don't have it: install.packages("writexl")
library(tidyverse)
library(writexl)

# ==============================================================================
# CONFIGURATION - Update output path
# ==============================================================================
OUTPUT_DIR <- "~/scMultiPreDICT_output/paper_figures/tables/"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. PREPARE THE DATA
# ==============================================================================
# NOTE: This script assumes you have already loaded 'df_delta' from your analysis
# df_delta should contain: Gene, Integrated_Analysis, rna_spearman, spearman, delta_spear
df_export <- df_delta %>%
  mutate(status = case_when( 
    delta_spear >  0.01 ~ "Improved", 
    delta_spear < -0.01 ~ "Worse", 
    TRUE ~ "No change" 
  )) %>%
  
  # 2. SELECT & RENAME COLUMNS
  select(
    Gene = Gene,
    Dataset = Integrated_Analysis,
    
    # The actual spearman correlations 
    RNA_Score = rna_spearman,       
    Integrated_Score = spearman,
    
    Delta = delta_spear,          # The difference
    Status = status               # The label (Improved/Worse/No Change)
  ) %>%
  
  # Sort so the biggest winners are at the top
  arrange(desc(Delta))

# 3. SPLIT INTO LISTS FOR EXCEL TABS
sheets_list <- list(
  "All_Genes"       = df_export,                        # Tab 1: Everything
  "Improved_Genes"  = df_export %>% filter(Status == "Improved"), # Tab 2: Winners
  "Worse_Genes"     = df_export %>% filter(Status == "Worse"),    # Tab 3: Losers
  "No_Change_Genes" = df_export %>% filter(Status == "No change") # Tab 4: Neutral
)

# 4. SAVE TO EXCEL
output_file <- file.path(OUTPUT_DIR, "Supplementary_Table_S1_Gene_Lists.xlsx")
write_xlsx(sheets_list, output_file)

message("Excel file created: ", output_file)