#!/usr/bin/env Rscript

#' Single-cell RNA-seq Processing Pipeline for COVID-19 PBMC Dataset
#' Dataset: Schulte-Schrepping et al., 2020 (PMID: 32810438)
#' Description: Process and analyze PBMC dynamics in COVID-19 patients
#' GEO: EGAS00001004571

# 1. Setup and Configuration -------------------------------------------
suppressPackageMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(Matrix.utils)
  library(sctransform)
  library(ComplexHeatmap)
  library(circlize)
  library(EpicTools)
  library(SingleR)
  library(scater)
  library(nichenetr)
})

# Set paths
BASE_DIR <- "/datasets"
OUTPUT_DIR <- file.path(BASE_DIR, "PMID32810438")
R_LIB_PATH <- "/R"

# Add custom R library path
.libPaths(c(R_LIB_PATH, .libPaths()))

# 2. Cell Type Subsetting -------------------------------------------
#' Subset Seurat Object by Cell Types
#' @param seurat_obj Seurat object with dmap.labels
#' @return List of Seurat objects for each cell type
subset_cell_types <- function(seurat_obj) {
  cell_type_definitions <- list(
    Megakaryocytes = "Megakaryocytes",
    Monocytes = "Monocytes",
    Erythroid = "Erythroid cells",
    T_cell = c("CD4+ T cells", "CD8+ T cells"),
    B_cell = "B cells",
    Other = c("NK cells", "Eosinophils", "Dendritic cells", 
              "Granulocytes", "Basophils", "HSCs", "MEPs", 
              "CMPs", "NK T cells")
  )
  
  subsets <- map(names(cell_type_definitions), function(type) {
    subset_data <- subset(seurat_obj, 
                         dmap.labels %in% cell_type_definitions[[type]]) %>%
      AddMetaData(metadata = type, col.name = "Cell_type")
    
    # Add standard metadata
    subset_data@meta.data <- subset_data@meta.data %>%
      mutate(
        Disease = "COVID-19",
        Data_ID = "Schulte-Schrepping et al.",
        Disease_group = paste(group_per_sample, 
                            disease_stage, 
                            outcome, 
                            sep = "_")
      )
    
    return(subset_data)
  })
  
  names(subsets) <- names(cell_type_definitions)
  return(subsets)
}

# 3. Cell Count Analysis --------------------------------------------
#' Calculate Cell Counts and Percentages
#' @param cell_subsets List of Seurat objects for each cell type
#' @return DataFrame with cell counts and percentages
analyze_cell_counts <- function(cell_subsets) {
  # Process each cell type
  cell_counts <- map(names(cell_subsets), function(type) {
    FetchData(cell_subsets[[type]], 
             vars = c("Data_ID", "sampleID", 
                     "Disease", "Disease_group", "Cell_type")) %>%
      group_by(Data_NO, Data_ID, sampleID, Disease, Disease_group, Cell_type) %>%
      count() %>%
      spread(Cell_type, n)
  })
  
  # Combine all cell types
  combined_counts <- reduce(cell_counts, full_join, 
                          by = c("Data_ID", "sampleID", 
                                "Disease", "Disease_group")) %>%
    rename(Sample_ID = sampleID)
  
  # Replace NA with 0
  combined_counts[is.na(combined_counts)] <- 0
  
  # Calculate percentages
  count_cols <- combined_counts %>% 
    select(-c(Data_NO:Disease_group))
  
  total_cells <- rowSums(count_cols)
  
  percentages <- map_df(count_cols, ~(. / total_cells) * 100) %>%
    rename_with(~paste0("PCT_", .), everything())
  
  results <- bind_cols(combined_counts, percentages)
  
  return(results)
}

# 4. Main Execution ------------------------------------------------
main <- function() {
  # Create output directory
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  message("Loading Seurat object...")
  seurat_obj <- readRDS(file.path(BASE_DIR, 
    "seurat_COVID19_PBMC_cohort1_10x_jonas_FG_2020-08-15.rds"))
  
  # Create cell type subsets
  message("Creating cell type subsets...")
  cell_subsets <- subset_cell_types(seurat_obj)
  
  # Save cell type subsets
  message("Saving cell type subsets...")
  walk2(cell_subsets, names(cell_subsets), function(subset, name) {
    saveRDS(subset, file.path(OUTPUT_DIR, 
            paste0("subset_", tolower(name), ".rds")))
  })
  
  # Analyze cell counts
  message("Analyzing cell counts...")
  cell_counts <- analyze_cell_counts(cell_subsets)
  
  # Save cell counts
  message("Saving cell count analysis...")
  write_tsv(cell_counts, 
            file.path(OUTPUT_DIR, "cell_type_counts.tsv"))
  
  message("Analysis pipeline completed successfully!")
}

# Run pipeline
if (!interactive()) {
  main()
}
