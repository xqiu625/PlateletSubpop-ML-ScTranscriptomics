#!/usr/bin/env Rscript

#' COVID-19 Single-cell RNA-seq Processing Pipeline
#' Dataset: Combes et al., 2021 (PMID: 33494096)
#' Description: Analysis of COVID-19 patient immune responses
#' GEO: GSE163668

# 1. Setup and Configuration -------------------------------------------
suppressPackageMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(MAST)
  library(scCATCH)
  library(SingleR)
  library(celldex)
  library(scater)
})

# Set paths and parameters
BASE_DIR <- "/datasets"
DATA_DIR <- file.path(BASE_DIR, "PMID33494096")
OUTPUT_DIR <- file.path(BASE_DIR, "PMID33494096")
R_LIB_PATH <- "/R"

# QC parameters
QC_PARAMS <- list(
  nFeature_low = 200,
  nFeature_high = 10000,
  mt_high = 20,
  min_cells = 3
)

# Sample and cell type definitions
SAMPLE_IDS <- paste0("GSM4995", 425:448)

CELL_TYPES <- list(
  Platelet = c('3', '6'),
  Monocytes = c('2', '13', '20', '5', '26'),
  Erythroid = c('18', '16'),
  T_cell = c('1', '9', '15', '24', '28', '10', '19', '32', '33', '4', '36'),
  B_cell = c('21', '12', '25', '31'),
  Other = c('0', '7', '8', '11', '14', '17', '22', '23', '29', '27', 
            '30', '34', '35', '37')
)

# 2. Data Loading Functions -------------------------------------------
#' Process Individual Sample
#' @param sample_id Sample identifier
#' @return Processed Seurat object
process_sample <- function(sample_id) {
  message(sprintf("Processing %s...", sample_id))
  
  # Load data files
  files <- list(
    barcodes = list.files(DATA_DIR, pattern = paste0(sample_id, ".*barcodes.tsv.gz"), full.names = TRUE),
    features = list.files(DATA_DIR, pattern = paste0(sample_id, ".*features.tsv.gz"), full.names = TRUE),
    matrix = list.files(DATA_DIR, pattern = paste0(sample_id, ".*matrix.mtx.gz"), full.names = TRUE)
  )
  
  # Read data
  mat <- readMM(files$matrix)
  features <- read.delim(files$features, header = FALSE)
  barcodes <- read.delim(files$barcodes, header = FALSE)
  
  # Set names
  colnames(mat) <- barcodes$V1
  rownames(mat) <- features$V2
  
  # Create and process Seurat object
  CreateSeuratObject(counts = mat, project = sample_id, 
                    min.cells = QC_PARAMS$min_cells, 
                    min.features = QC_PARAMS$nFeature_low) %>%
    PercentageFeatureSet("^MT-", col.name = "percent.mt") %>%
    subset(nCount_RNA > 500 & 
           nFeature_RNA > QC_PARAMS$nFeature_low & 
           nFeature_RNA < QC_PARAMS$nFeature_high & 
           percent.mt < QC_PARAMS$mt_high) %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", 
                        nfeatures = 4000, 
                        verbose = FALSE)
}

# 3. Integration Functions -------------------------------------------
#' Integrate Multiple Samples
#' @param sample_list List of Seurat objects
#' @return Integrated Seurat object
integrate_samples <- function(sample_list) {
  # SCTransform each dataset
  sample_list <- map(sample_list, function(x) {
    NormalizeData(x) %>%
      SCTransform(vars.to.regress = "percent.mt")
  })
  
  # Integration
  features <- SelectIntegrationFeatures(sample_list, nfeatures = 3000)
  sample_list <- PrepSCTIntegration(sample_list, anchor.features = features)
  
  anchors <- FindIntegrationAnchors(sample_list,
                                   normalization.method = "SCT",
                                   anchor.features = features)
  
  IntegrateData(anchors, normalization.method = "SCT")
}

# 4. Analysis Functions ---------------------------------------------
#' Process Integrated Data
#' @param seurat_obj Integrated Seurat object
#' @return Processed Seurat object
process_integrated_data <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  
  seurat_obj %>%
    FindVariableFeatures(selection.method = "vst", 
                        nfeatures = 4000, 
                        verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:50, verbose = FALSE) %>%
    FindNeighbors(dims = 1:50, verbose = FALSE) %>%
    FindClusters(resolution = 1, verbose = FALSE)
}

# 5. Cell Type Analysis -------------------------------------------
#' Subset Cell Types and Calculate Proportions
#' @param seurat_obj Processed Seurat object
#' @return List of cell type subsets and proportions
analyze_cell_types <- function(seurat_obj) {
  # Create subsets
  subsets <- map(names(CELL_TYPES), function(type) {
    subset_obj <- subset(seurat_obj, idents = CELL_TYPES[[type]])
    subset_obj@meta.data <- subset_obj@meta.data %>%
      mutate(
        Cell_type = type,
        Disease = "COVID-19",
        Data_ID = "Combes et al.",
        Disease_group = paste(covid_status, Qualitative_score,
                            ICU_vs_FLOOR, Outcome, sep = "_")
      )
    return(subset_obj)
  }) %>%
    set_names(names(CELL_TYPES))
  
  # Calculate cell counts
  counts <- map(subsets, function(subset) {
    FetchData(subset, 
             vars = c("Data_ID", "orig.ident", 
                     "Disease", "Disease_group", "Cell_type")) %>%
      group_by(Data_ID, orig.ident, Disease, Disease_group, Cell_type) %>%
      count() %>%
      spread(Cell_type, n)
  }) %>%
    reduce(full_join, by = c("Data_ID", "orig.ident", 
                            "Disease", "Disease_group"))
  
  # Process counts and calculate percentages
  counts[is.na(counts)] <- 0
  counts <- rename(counts, Sample_ID = orig.ident)
  
  total_cells <- rowSums(select(counts, -c(Data_ID:Disease_group)))
  percentages <- map_df(select(counts, -c(Data_ID:Disease_group)), 
                       ~(. / total_cells) * 100) %>%
    rename_with(~paste0("PCT_", .), everything())
  
  list(
    subsets = subsets,
    counts = bind_cols(counts, percentages)
  )
}

# 6. Main Execution ------------------------------------------------
main <- function() {
  # Create output directory
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Process samples
  message("Processing samples...")
  samples <- map(SAMPLE_IDS, process_sample)
  names(samples) <- SAMPLE_IDS
  
  # Integrate samples
  message("Integrating samples...")
  integrated_data <- integrate_samples(samples)
  
  # Process integrated data
  message("Processing integrated data...")
  processed_data <- process_integrated_data(integrated_data)
  
  # Analyze cell types
  message("Analyzing cell types...")
  results <- analyze_cell_types(processed_data)
  
  # Save results
  message("Saving results...")
  saveRDS(processed_data, file.path(OUTPUT_DIR, "complete_analysis.rds"))
  
  walk2(results$subsets, names(results$subsets), function(subset, name) {
    saveRDS(subset, file.path(OUTPUT_DIR, 
            paste0("subset_", tolower(name), ".rds")))
  })
  
  write_tsv(results$counts, file.path(OUTPUT_DIR, "cell_type_counts.tsv"))
  
  message("Analysis pipeline completed successfully!")
}

# Run pipeline
if (!interactive()) {
  main()
}
