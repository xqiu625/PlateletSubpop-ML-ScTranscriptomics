#!/usr/bin/env Rscript

#' COVID-19 Single-cell Atlas Processing Pipeline
#' Dataset: Ren et al., 2021 (PMID: 33657410)
#' Description: Large-scale single-cell atlas of COVID-19 immune features
#' GEO: GSE158055

# 1. Setup and Configuration -------------------------------------------
suppressPackageMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(SingleCellExperiment)
})

# Set paths
BASE_DIR <- "/datasets"
DATA_DIR <- file.path(BASE_DIR, "PMID33657410")
OUTPUT_DIR <- file.path(BASE_DIR, "PMID33657410")
R_LIB_PATH <- "/R"

.libPaths(c(R_LIB_PATH, .libPaths()))

# 2. Data Loading Functions -------------------------------------------
#' Load and Combine Matrix Parts
#' @param base_path Base path for data files
#' @return List of Seurat objects
load_matrix_parts <- function() {
  message("Loading matrix parts...")
  
  # Define paths for both parts
  paths <- list(
    part1 = file.path(DATA_DIR, "filtered_gene_bc_matrices/GSE158055_covid19_part1_folder"),
    part2 = file.path(DATA_DIR, "filtered_gene_bc_matrices/GSE158055_covid19_part2_folder")
  )
  
  # Load both parts
  seurat_objects <- map(paths, function(path) {
    data <- Read10X(path, gene.column = 1)
    CreateSeuratObject(counts = data, project = basename(path))
  })
  
  return(seurat_objects)
}

#' Load and Process Metadata
#' @return List containing sample and cell annotations
load_metadata <- function() {
  message("Loading metadata...")
  
  # Load cell annotations
  cell_annot <- read_csv(file.path(DATA_DIR, "GSE158055_cell_annotation.csv"))
  
  # Load and filter sample metadata
  sample_meta <- read_delim(file.path(DATA_DIR, "GSE158055_sample_metadata.txt")) %>%
    filter(
      str_detect(`characteristics..Sample.type`, "PBMC"),
      str_detect(`characteristics..Sample.time`, "progression|control")
    )
  
  return(list(
    cell_annot = cell_annot,
    sample_meta = sample_meta
  ))
}

# 3. Data Processing ------------------------------------------------
#' Combine and Process Data
#' @param seurat_parts List of Seurat objects
#' @param metadata List of metadata
#' @return Processed Seurat object
process_data <- function(seurat_parts, metadata) {
  message("Processing data...")
  
  # Filter for PBMC samples
  pbmc_samples <- metadata$sample_meta$Sample.name
  pbmc_cells <- which(seurat_parts$part1$sampleID %in% pbmc_samples)
  
  # Subset both parts
  part1_pbmc <- seurat_parts$part1[, pbmc_cells]
  part2_pbmc <- seurat_parts$part2[, pbmc_cells]
  
  # Combine counts
  combined_data <- part1_pbmc
  combined_data@assays$RNA@counts <- combined_data@assays$RNA@counts + 
                                   part2_pbmc@assays$RNA@counts
  combined_data@assays$RNA@data <- combined_data@assays$RNA@data + 
                                  part2_pbmc@assays$RNA@data
  
  # Add metadata
  combined_data$sampleID <- factor(combined_data$sampleID)
  combined_data$celltype <- metadata$cell_annot$celltype[match(colnames(combined_data), 
                                                              rownames(metadata$cell_annot))]
  combined_data$majortype <- metadata$cell_annot$majorType[match(colnames(combined_data), 
                                                                rownames(metadata$cell_annot))]
  combined_data$batch <- metadata$sample_meta$`characteristics...Datasets`[match(combined_data$sampleID,
                                                                               metadata$sample_meta$Sample.name)]
  
  # Quality control
  combined_data <- combined_data %>%
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
    PercentageFeatureSet(pattern = "^RPS", col.name = "percent.rps") %>%
    PercentageFeatureSet(pattern = "^RPL", col.name = "percent.rpl") %>%
    PercentageFeatureSet(pattern = "^RNA\\d8S5", col.name = "percent.rrna") %>%
    NormalizeData(verbose = FALSE)
  
  return(combined_data)
}

# 4. Cell Type Subsetting ------------------------------------------
#' Subset Data by Cell Types
#' @param seurat_obj Processed Seurat object
#' @return List of cell type-specific Seurat objects
subset_cell_types <- function(seurat_obj) {
  message("Subsetting cell types...")
  
  cell_types <- list(
    Platelet = "Mega",
    Monocyte = "Mono",
    T_cell = c("CD4", "CD8"),
    B_cell = "B",
    Other = c("Macro", "DC", "Neu", "Plasma", "NK")
  )
  
  subsets <- map(names(cell_types), function(type) {
    subset_obj <- subset(seurat_obj, majortype %in% cell_types[[type]])
    subset_obj@meta.data <- subset_obj@meta.data %>%
      mutate(
        Cell_type = type,
        Disease = "COVID-19",
        Data_ID = "Ren et al."
      )
    return(subset_obj)
  })
  
  names(subsets) <- names(cell_types)
  return(subsets)
}

# 5. Main Execution -----------------------------------------------
main <- function() {
  # Create output directory
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  seurat_parts <- load_matrix_parts()
  metadata <- load_metadata()
  
  # Process data
  processed_data <- process_data(seurat_parts, metadata)
  
  # Create cell type subsets
  cell_subsets <- subset_cell_types(processed_data)
  
  # Save results
  message("Saving results...")
  saveRDS(processed_data, file.path(OUTPUT_DIR, "complete_dataset.rds"))
  
  walk2(cell_subsets, names(cell_subsets), function(subset, name) {
    saveRDS(subset, file.path(OUTPUT_DIR, 
            paste0("subset_", tolower(name), ".rds")))
  })
  
  message("Analysis pipeline completed successfully!")
}

# Run pipeline
if (!interactive()) {
  main()
}
