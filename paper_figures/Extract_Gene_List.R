# ==============================================================================
# SETUP
# ==============================================================================
# Install if you don't have it: install.packages("writexl")
library(tidyverse)
library(writexl)

# ==============================================================================
# 1. PREPARE THE DATA
# ==============================================================================

# 1. ADD STATUS COLUMN TO THE MAIN DATAFRAME
# Obtain Gene Status based on whether they improve on performance by adding atac to rna
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
write_xlsx(sheets_list, "~/PROJECTS/2025.PERTURBATION_MODELLING/Combined_Figures_Finale/Tables/Supplementary_Table_S1_Gene_Lists.xlsx")

message("Excel file created with actual scores!")