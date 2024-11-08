#!/usr/bin/env Rscript

#' Platelet Dataset Integration Pipeline Using Harmony
#' Description: Integrate multiple platelet datasets using Harmony batch correction


# 1. Setup and Configuration -------------------------------------------
suppressPackageMessages({
  library(Seurat)
  library(harmony)
  library(tidyverse)
  library(cowplot)
})

# Set paths
BASE_DIR <- "~/datasets"
INPUT_DIR <- file.path(BASE_DIR, "Platelets")
OUTPUT_DIR <- file.path(BASE_DIR, "harmony_int")
PLOT_DIR <- file.path(BASE_DIR, "Harmony")

# Analysis parameters
ANALYSIS_PARAMS <- list(
  variable_features = 2000,
  pca_dims = 30,
  harmony_dims = 30,
  cluster_resolution = 0.8
)

# 2. Data Loading Functions -------------------------------------------
#' Load and Process All Datasets
#' @return List containing processed matrices and metadata
load_datasets <- function() {
  # Load all RDS files
  message("Loading datasets...")
  datasets <- list.files(INPUT_DIR, pattern = "*.rds", full.names = TRUE) %>%
    map(readRDS)
  
  # Extract count matrices
  matrices <- map(datasets, function(x) x@assays$RNA@counts)
  
  # Find common genes
  common_genes <- reduce(map(matrices, rownames), intersect)
  
  # Subset matrices to common genes
  matrices_filtered <- map(matrices, function(x) x[common_genes,])
  
  # Extract metadata
  metadata <- map_dfr(datasets, function(x) {
    x@meta.data %>%
      select(Cell_ID, Sample_ID, Data_ID)
  })
  
  list(
    matrices = matrices_filtered,
    metadata = metadata,
    common_genes = common_genes
  )
}

# 3. Integration Functions --------------------------------------------
#' Create Integrated Seurat Object
#' @param data_list List containing matrices and metadata
#' @return Integrated Seurat object
create_integrated_object <- function(data_list) {
  message("Creating integrated Seurat object...")
  
  # Combine matrices
  combined_matrix <- do.call(cbind, data_list$matrices)
  
  # Create and process Seurat object
  CreateSeuratObject(counts = combined_matrix) %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(
      selection.method = "vst", 
      nfeatures = ANALYSIS_PARAMS$variable_features, 
      verbose = FALSE
    ) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(
      pc.genes = .$var.genes, 
      npcs = ANALYSIS_PARAMS$pca_dims, 
      verbose = FALSE
    )
}

#' Run Harmony Integration
#' @param seurat_obj Seurat object
#' @return Harmony-integrated Seurat object
run_harmony_integration <- function(seurat_obj) {
  message("Running Harmony integration...")
  
  # Create plot directory if needed
  dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Setup plot
  png(file.path(PLOT_DIR, "Harmony_convergence_Data_ID_Platelets.png"),
      width = 4 * 300, height = 4 * 300,
      units = "px", res = 300, type = 'cairo')
  
  # Run Harmony
  harmony_obj <- seurat_obj %>%
    RunHarmony("Data_ID", plot_convergence = TRUE)
  
  dev.off()
  
  harmony_obj
}

#' Run Downstream Analysis
#' @param harmony_obj Harmony-integrated Seurat object
#' @return Processed Seurat object
run_downstream_analysis <- function(harmony_obj) {
  message("Running downstream analysis...")
  
  harmony_obj %>%
    RunUMAP(
      reduction = "harmony", 
      dims = 1:ANALYSIS_PARAMS$harmony_dims
    ) %>%
    FindNeighbors(
      reduction = "harmony", 
      dims = 1:ANALYSIS_PARAMS$harmony_dims
    ) %>%
    FindClusters(
      resolution = ANALYSIS_PARAMS$cluster_resolution
    )
}

# 4. Main Execution ------------------------------------------------
main <- function() {
  # Create output directory
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  data_list <- load_datasets()
  
  # Create initial Seurat object
  seurat_obj <- create_integrated_object(data_list)
  
  # Add metadata
  seurat_obj@meta.data <- bind_cols(
    seurat_obj@meta.data,
    data_list$metadata
  )
  
  # Run Harmony integration
  harmony_obj <- run_harmony_integration(seurat_obj)
  
  # Save harmony object
  message("Saving Harmony integration results...")
  saveRDS(harmony_obj, 
          file.path(OUTPUT_DIR, "302005_Platelets.rds"))
  
  # Run clustering and UMAP
  clustered_obj <- run_downstream_analysis(harmony_obj)
  
  # Save final object
  message("Saving final clustered object...")
  saveRDS(clustered_obj, 
          file.path(OUTPUT_DIR, "302005_Platelets_cluster.rds"))
  
  message("Integration pipeline completed successfully!")
}

# Run pipeline
if (!interactive()) {
  main()
}
