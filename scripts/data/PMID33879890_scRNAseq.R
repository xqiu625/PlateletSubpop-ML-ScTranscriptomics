#!/usr/bin/env Rscript

#' COVID-19 Multi-omics Immune Response Analysis Pipeline
#' Dataset: Stephenson et al., 2021 (PMID: 33879890)
#' Description: Multi-omics single-cell longitudinal analysis of COVID-19 patients
#' Data: COVID-19 Portal (E-MTAB-10026)


# 1. Setup and Configuration -------------------------------------------
suppressPackageMessages({
  library(Seurat)
  library(SeuratDisk)
  library(tidyverse)
  library(Matrix)
})

# Set paths
BASE_DIR <- "/datasets"
DATA_DIR <- file.path(BASE_DIR, "PMID33879890")
OUTPUT_DIR <- file.path(BASE_DIR, "PMID33879890")
R_LIB_PATH <- "/R"

.libPaths(c(R_LIB_PATH, .libPaths()))

# Define cell type clusters
CELL_TYPES <- list(
  Platelet = "Platelets",
  Monocytes = c("CD14", "CD16", "Mono_prolif"),
  T_cell = c("CD4", "CD8", "Lymph_prolif", "Treg", "gdT", "MAIT"),
  B_cell = c("B_cell", "Plasmablast"),
  Erythroid = "RBC",
  Other = c("DCs", "HSC", "NK_16hi", "NK_56hi", "pDC")
)

# 2. Data Loading Functions -------------------------------------------
#' Convert h5ad to h5seurat
#' @param input_file Input h5ad file path
#' @param output_file Output h5seurat file path
convert_h5ad_to_seurat <- function(input_file, output_file) {
  Convert(input_file, dest = output_file, overwrite = TRUE)
}

#' Load and Process H5AD Data
#' @param file_path Path to h5ad file
#' @return Seurat object
load_h5ad_data <- function(file_path) {
  h5seurat_path <- sub("\\.h5ad$", ".h5seurat", file_path)
  
  # Convert if needed
  if (!file.exists(h5seurat_path)) {
    convert_h5ad_to_seurat(file_path, h5seurat_path)
  }
  
  LoadH5Seurat(h5seurat_path)
}

# 3. Cell Type Processing -------------------------------------------
#' Process Cell Type Subset
#' @param seurat_obj Seurat object
#' @param cell_type Cell type category
#' @param cell_type_list List of specific cell types
#' @return Processed Seurat object
subset_cell_type <- function(seurat_obj, cell_type, cell_type_list) {
  subset_obj <- subset(seurat_obj, initial_clustering %in% cell_type_list)
  
  # Add metadata
  subset_obj@meta.data <- subset_obj@meta.data %>%
    mutate(
      Cell_type = cell_type,
      Disease = "COVID-19",
      Data_ID = "Stephenson et al.",
      Disease_group = paste(Collection_Day,
                          Status_on_day_collection_summary,
                          Outcome, 
                          sep = "_")
    )
  
  return(subset_obj)
}

# 4. Cell Count Analysis --------------------------------------------
#' Calculate Cell Counts for a Subset
#' @param seurat_obj Processed Seurat object
#' @return DataFrame with cell counts
calculate_cell_counts <- function(seurat_obj) {
  FetchData(seurat_obj, 
           vars = c("Data_ID", "sample_id", 
                   "Disease", "Disease_group", "Cell_type")) %>%
    group_by(Data_ID, sample_id, Disease, Disease_group, Cell_type) %>%
    count() %>%
    spread(Cell_type, n) %>%
    rename(Sample_ID = sample_id)
}

#' Process Cell Counts
#' @param cell_subsets List of cell type Seurat objects
#' @return DataFrame with combined counts and percentages
process_cell_counts <- function(cell_subsets) {
  # Calculate counts for each subset
  counts <- map(cell_subsets, calculate_cell_counts)
  
  # Combine all counts
  combined_counts <- reduce(counts, full_join, by = c("Sample_ID")) %>%
    mutate(
      Disease_group = coalesce(!!!select(., matches("Disease_group"))),
      Data_ID = "Stephenson et al.",
      Disease = "COVID-19"
    ) %>%
    select(Data_ID, Sample_ID, Disease, Disease_group,
           Monocyte, Platelet, Erythroid, T_cell, B_cell, Other_cell)
  
  # Replace NA with 0
  combined_counts[is.na(combined_counts)] <- 0
  
  # Calculate percentages
  count_cols <- combined_counts %>% 
    select(-c(Data_ID:Disease_group))
  
  total_cells <- rowSums(count_cols)
  percentages <- map_df(count_cols, ~(. / total_cells) * 100) %>%
    rename_with(~paste0("PCT_", .), everything())
  
  bind_cols(combined_counts, percentages)
}

# 5. Main Execution ------------------------------------------------
main <- function() {
  # Create output directory
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  message("Loading dataset...")
  seurat_data <- load_h5ad_data(file.path(DATA_DIR, 
                               "covid_portal_210320_with_raw.h5ad"))
  
  # Process cell types
  message("Processing cell types...")
  cell_subsets <- map(names(CELL_TYPES), function(type) {
    message(sprintf("Processing %s cells...", type))
    subset_cell_type(seurat_data, type, CELL_TYPES[[type]])
  }) %>%
    set_names(names(CELL_TYPES))
  
  # Save cell type subsets
  message("Saving cell type subsets...")
  walk2(cell_subsets, names(cell_subsets), function(subset, name) {
    saveRDS(subset, file.path(OUTPUT_DIR, 
            paste0("subset_", tolower(name), ".rds")))
  })
  
  # Calculate and save cell counts
  message("Calculating cell counts...")
  cell_counts <- process_cell_counts(cell_subsets)
  write_tsv(cell_counts, file.path(OUTPUT_DIR, "cell_type_counts.tsv"))
  
  message("Analysis pipeline completed successfully!")
}

# Run pipeline
if (!interactive()) {
  main()
}
