#!/usr/bin/env Rscript

#' SLE Single-cell RNA-seq Processing Pipeline
#' Dataset: Mistry et al., 2019 (PMID: 31754025)
#' Description: Analysis of SLE patient immune responses using scRNA-seq
#' GEO: GSE142016

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
DATA_DIR <- file.path(BASE_DIR, "PMID31754025")
OUTPUT_DIR <- file.path(BASE_DIR, "PMID31754025")
R_LIB_PATH <- "/R"

# QC parameters
QC_PARAMS <- list(
  nFeature_low = 200,
  nFeature_high = 10000,
  mt_high = 20,
  min_cells = 3
)

# Sample and cell type definitions
SAMPLE_IDS <- c("SLE1", "SLE2", "SLE3")

CELL_TYPES <- list(
  Platelet = c('24'),
  Monocytes = c('0', '1', '2', '4', '8', '15', '17', '18', '19', '23', '29'),
  T_cell = c('3', '5', '7', '11', '13', '22', '26'),
  B_cell = c('6', '10', '14', '16', '30'),
  Other = c('9', '12', '20', '21', '25', '27', '28', '31')
)

# 2. Data Loading Functions -------------------------------------------
#' Process Individual Sample
#' @param sample_id Sample identifier
#' @return Processed Seurat object
process_sample <- function(sample_id) {
  message(sprintf("Processing %s...", sample_id))
  
  # Find files
  files <- list(
    barcodes = list.files(pattern = paste0('GSM.*', sample_id, '_PBMC_raw_barcodes.tsv.gz')),
    features = list.files(pattern = paste0('GSM.*', sample_id, '_PBMC_raw_genes.tsv.gz')),
    matrix = list.files(pattern = paste0('GSM.*', sample_id, '_PBMC_raw_matrix.mtx.gz'))
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
                        nfeatures = 3000, 
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
    FindClusters(resolution = 1, verbose = FALSE) %>%
    FindAllMarkers(assay = 'RNA', 
                  only.pos = TRUE, 
                  test.use = 'MAST', 
                  verbose = FALSE)
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
        Disease = "SLE",
        Data_ID = "Mistry et al.",
        Disease_group = "SLE"
      )
    return(subset_obj)
  }) %>%
    set_names(names(CELL_TYPES))
  
  # Calculate cell counts
  counts <- map(subsets, function(subset) {
    FetchData(subset, 
             vars = c("Data_ID", "orig.ident", 
                     "Disease", "Disease_group", "Cell_type")) %>%
      mutate(orig.ident = trimws(orig.ident)) %>%
      group_by(Data_ID, orig.ident, Disease, Disease_group, Cell_type) %>%
      count() %>%
      spread(Cell_type, n)
  }) %>%
    reduce(full_join, by = c("Data_ID", "orig.ident", 
                            "Disease", "Disease_group"))
  
  # Process counts and calculate percentages
  counts[is.na(counts)] <- 0
  counts <- rename(counts, Sample_ID = orig.ident) %>%
    mutate(Sample_ID = as.character(Sample_ID))
  
  count_cols <- select(counts, -c(Data_ID:Disease_group))
  total_cells <- rowSums(count_cols)
  
  percentages <- map_df(count_cols, ~(. / total_cells) * 100) %>%
    rename_with(~paste0("PCT_", .), everything())
  
  # Add empty Erythroid columns
  final_counts <- bind_cols(counts, percentages) %>%
    mutate(
      Erythroid = NA_real_,
      PCT_Erythroid = NA_real_
    ) %>%
    select(Data_ID, Sample_ID, Disease, Disease_group,
           Monocyte, Platelet, Erythroid, T_cell, B_cell, Other_cell,
           starts_with("PCT_"))
  
  list(
    subsets = subsets,
    counts = final_counts
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
