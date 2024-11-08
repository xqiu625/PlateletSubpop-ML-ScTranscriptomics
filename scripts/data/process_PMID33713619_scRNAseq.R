#!/usr/bin/env Rscript

#' COVID-19 Temporal Immune Response Analysis Pipeline
#' Dataset: Liu et al., 2021 (PMID: 33713619)
#' Description: Time-resolved systems immunology of COVID-19 progression and recovery
#' Data: COVID-19 Cell Atlas (https://www.covid19cellatlas.org/)

# 1. Setup and Configuration -------------------------------------------
suppressPackageMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(SingleCellExperiment)
})

# Set paths
BASE_DIR <- "/datasets"
DATA_DIR <- file.path(BASE_DIR, "PMID33713619")
OUTPUT_DIR <- file.path(BASE_DIR, "PMID33713619")
R_LIB_PATH <- "/R"

.libPaths(c(R_LIB_PATH, .libPaths()))

# 2. Cell Type Definitions ------------------------------------------
CELL_TYPES <- list(
  Platelet = "platelet",
  Monocyte = c("non-classical monocyte", "classical monocyte", "intermediate monocyte"),
  T_cell = c("gamma-delta T cell", "regulatory T cell", 
             "CD4-positive, alpha-beta memory T cell",
             "CD8-positive, alpha-beta memory T cell",
             "naive CD8+ T cell", "naive CD4+ T cell",
             "mucosal invariant T cell (MAIT)",
             "TissueResMemT", "double-positive T cell (DPT)",
             "double negative T cell (DNT)", "TCRVbeta13.1pos"),
  B_cell = c("naive B cell", "memory B cell", "plasmablast"),
  Other = c("NK_CD16hi", "NK_CD56loCD16lo", "plasmacytoid dendritic cell",
            "conventional dendritic cell", "NK_CD56hiCD16lo", "granulocyte")
)

# 3. Data Processing Functions --------------------------------------
#' Process Cell Type Subset
#' @param seurat_obj Seurat object
#' @param cell_type Cell type category
#' @param cell_type_list List of specific cell types
#' @return Processed Seurat object
process_cell_type <- function(seurat_obj, cell_type, cell_type_list) {
  subset_obj <- subset(seurat_obj, cell_type %in% cell_type_list)
  
  # Add metadata
  subset_obj@meta.data <- subset_obj@meta.data %>%
    mutate(
      Cell_type = cell_type,
      Disease = "COVID-19",
      Data_ID = "Liu et al.",
      Disease_group = paste(severity, ever_admitted_to_icu, outcome, sep = "_"),
      donor = paste(donor, timepoint, sep = "_")
    )
  
  return(subset_obj)
}

#' Calculate Cell Counts
#' @param seurat_obj Processed Seurat object
#' @return DataFrame with cell counts
calculate_cell_counts <- function(seurat_obj) {
  FetchData(seurat_obj, 
           vars = c("Data_ID", "donor", 
                   "Disease", "Disease_group", "Cell_type")) %>%
    mutate(donor = trimws(donor)) %>%
    group_by(Data_ID, donor, Disease, Disease_group, Cell_type) %>%
    count(donor) %>%
    spread(Cell_type, n) %>%
    rename(Sample_ID = donor) %>%
    mutate(Sample_ID = as.character(Sample_ID))
}

#' Process All Cell Types
#' @param innate_cells Innate cells Seurat object
#' @param adaptive_cells Adaptive cells Seurat object
#' @return List of processed cell type objects
process_all_cell_types <- function(innate_cells, adaptive_cells) {
  # Process innate cells
  innate_subsets <- list(
    Platelet = process_cell_type(innate_cells, "Platelet", CELL_TYPES$Platelet),
    Monocyte = process_cell_type(innate_cells, "Monocyte", CELL_TYPES$Monocyte),
    Other = process_cell_type(innate_cells, "Other", CELL_TYPES$Other)
  )
  
  # Process adaptive cells
  adaptive_subsets <- list(
    T_cell = process_cell_type(adaptive_cells, "T_cell", CELL_TYPES$T_cell),
    B_cell = process_cell_type(adaptive_cells, "B_cell", CELL_TYPES$B_cell)
  )
  
  c(innate_subsets, adaptive_subsets)
}

#' Calculate Combined Cell Counts
#' @param cell_subsets List of cell type Seurat objects
#' @return DataFrame with combined counts and percentages
calculate_combined_counts <- function(cell_subsets) {
  # Calculate counts for each cell type
  counts <- map(cell_subsets, calculate_cell_counts)
  
  # Combine all counts
  combined_counts <- reduce(counts, full_join, by = c("Sample_ID")) %>%
    mutate(
      Disease_group = coalesce(!!!select(., matches("Disease_group"))),
      Data_ID = coalesce(!!!select(., matches("Data_ID")), "Liu et al."),
      Disease = coalesce(!!!select(., matches("Disease")), "COVID-19")
    ) %>%
    select(Data_ID, Sample_ID, Disease, Disease_group,
           Monocyte, Platelet, T_cell, B_cell, Other_cell)
  
  # Replace NA with 0
  combined_counts[is.na(combined_counts)] <- 0
  
  # Calculate percentages
  count_cols <- combined_counts %>% 
    select(-c(Data_ID:Disease_group))
  total_cells <- rowSums(count_cols)
  
  percentages <- map_df(count_cols, ~(. / total_cells) * 100) %>%
    rename_with(~paste0("PCT_", .), everything())
  
  # Add empty Erythroid columns
  result <- bind_cols(combined_counts, percentages) %>%
    mutate(
      Erythroid = NA_real_,
      PCT_Erythroid = NA_real_,
      Disease_group = str_replace_all(Disease_group, 
                                    c("False" = "Floor", "True" = "ICU"))
    ) %>%
    select(Data_ID, Sample_ID, Disease, Disease_group,
           Monocyte, Platelet, Erythroid, T_cell, B_cell, Other_cell,
           starts_with("PCT_"))
  
  return(result)
}

# 4. Main Execution -----------------------------------------------
main <- function() {
  # Create output directory
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  message("Loading data...")
  innate_cells <- readRDS(file.path(DATA_DIR, "Innate_Cells.rds"))
  adaptive_cells <- readRDS(file.path(DATA_DIR, "Adaptive_Cells.rds"))
  
  # Process cell types
  message("Processing cell types...")
  cell_subsets <- process_all_cell_types(innate_cells, adaptive_cells)
  
  # Save cell type subsets
  message("Saving cell type subsets...")
  walk2(cell_subsets, names(cell_subsets), function(subset, name) {
    saveRDS(subset, file.path(OUTPUT_DIR, 
            paste0("subset_", tolower(name), ".rds")))
  })
  
  # Calculate and save cell counts
  message("Calculating cell counts...")
  cell_counts <- calculate_combined_counts(cell_subsets)
  write_tsv(cell_counts, file.path(OUTPUT_DIR, "cell_type_counts.tsv"))
  
  message("Analysis pipeline completed successfully!")
}

# Run pipeline
if (!interactive()) {
  main()
}
