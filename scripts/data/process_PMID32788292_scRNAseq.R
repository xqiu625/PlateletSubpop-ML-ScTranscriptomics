#!/usr/bin/env Rscript

#' COVID-19 Single-cell RNA-seq Processing Pipeline
#' Dataset: Arunachalam et al., 2020 (PMID: 32788292)
#' Description: Systems-level immunological assessment of COVID-19 severity
#' GEO: GSE155673

# 1. Setup and Configuration -------------------------------------------
suppressPackageMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(MAST)
  library(scCATCH)
  library(SingleCellExperiment)
  library(celldex)
  library(SingleR)
  library(scater)
})

# Set paths and parameters
BASE_DIR <- "/datasets"
OUTPUT_DIR <- file.path(BASE_DIR, "PMID32788292")
R_LIB_PATH <- "/R"

# QC parameters
QC_PARAMS <- list(
  nFeature_low = 200,
  nFeature_high = 10000,
  mt_high = 20,
  min_cells = 3
)

# 2. Data Loading Functions -------------------------------------------
#' Load and Process Individual Sample
#' @param sample_id Patient ID
#' @return Seurat object
process_sample <- function(sample_id) {
  message(sprintf("Processing sample: %s", sample_id))
  
  # Load data
  barcode.path <- sprintf("GSE155673_%s_barcodes.tsv.gz", sample_id)
  features.path <- "GSE155673_features.tsv.gz"
  matrix.path <- sprintf("GSE155673_%s_matrix.mtx.gz", sample_id)
  
  # Read files
  mat <- readMM(matrix.path)
  features <- read_delim(features.path, delim = "\t", col_names = FALSE)
  barcodes <- read_delim(barcode.path, delim = "\t", col_names = FALSE)
  
  # Set names
  colnames(mat) <- barcodes$X1
  rownames(mat) <- features$X2
  
  # Create Seurat object with QC
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

# 3. Integration and Analysis ----------------------------------------
#' Integrate Multiple Samples
#' @param seurat_list List of Seurat objects
#' @return Integrated Seurat object
integrate_samples <- function(seurat_list) {
  # SCTransform each dataset
  seurat_list <- map(seurat_list, function(x) {
    SCTransform(x, vars.to.regress = c("percent.mt"))
  })
  
  # Integration steps
  features <- SelectIntegrationFeatures(seurat_list, nfeatures = 3000)
  seurat_list <- PrepSCTIntegration(seurat_list, anchor.features = features)
  anchors <- FindIntegrationAnchors(seurat_list,
                                   normalization.method = "SCT",
                                   anchor.features = features)
  
  # Integrate data
  integrated <- IntegrateData(anchors, normalization.method = "SCT")
  
  # Analysis steps
  integrated %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:50, verbose = FALSE) %>%
    FindNeighbors(dims = 1:50, verbose = FALSE) %>%
    FindClusters(resolution = 1, verbose = FALSE)
}

# 4. Cell Type Annotation -------------------------------------------
#' Perform Cell Type Annotation using Multiple Methods
#' @param seurat_obj Seurat object
#' @return Annotated Seurat object
annotate_cell_types <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Find markers
  seurat_obj@misc$markers <- FindAllMarkers(seurat_obj, 
                                           assay = 'RNA',
                                           only.pos = TRUE, 
                                           test.use = 'MAST')
  
  # Multiple annotation methods
  seurat_obj@misc$scCATCH_cellType <- annotate_with_scCATCH(seurat_obj@misc$markers)
  seurat_obj@misc$Canonical_cellType <- annotate_with_canonical_markers(seurat_obj@misc$markers)
  seurat_obj@misc$SingleR_cellType <- annotate_with_singler(seurat_obj)
  
  return(seurat_obj)
}

# 5. Cell Type Subsetting ------------------------------------------
CELL_TYPE_CLUSTERS <- list(
  Platelet = c('0', '6', '12', '21', '40'),
  Monocytes = c('1', '10', '16', '14', '31', '38', '29'),
  Erythroid = '19',
  T_cell = c('2', '3', '5', '36', '26', '33', '8', '18', '32', '9'),
  B_cell = c('11', '13', '25', '37', '7', '24', '39'),
  Other = c('4', '15', '17', '20', '22', '34', '23', '35', '42', 
            '27', '28', '30')
)

subset_cell_types <- function(seurat_obj) {
  map(names(CELL_TYPE_CLUSTERS), function(cell_type) {
    subset_obj <- subset(seurat_obj, idents = CELL_TYPE_CLUSTERS[[cell_type]])
    
    # Add metadata
    subset_obj@meta.data <- subset_obj@meta.data %>%
      mutate(
        Cell_type = cell_type,
        Disease = "COVID-19",
        Data_ID = "Arunachalam et al."
      )
    
    return(subset_obj)
  }) %>% set_names(names(CELL_TYPE_CLUSTERS))
}

# 6. Cell Count Analysis ------------------------------------------
analyze_cell_counts <- function(cell_subsets) {
  # Get counts for each cell type
  cell_counts <- map(names(cell_subsets), function(type) {
    FetchData(cell_subsets[[type]], 
             vars = c("Data_ID", "orig.ident", 
                     "Disease", "Disease_group", "Cell_type")) %>%
      group_by(Data_ID, orig.ident, Disease, Disease_group, Cell_type) %>%
      count() %>%
      spread(Cell_type, n)
  })
  
  # Combine and process counts
  combined_counts <- reduce(cell_counts, full_join, 
                          by = c("Data_ID", "orig.ident", 
                                "Disease", "Disease_group")) %>%
    rename(Sample_ID = orig.ident) %>%
    replace(is.na(.), 0)
  
  # Calculate percentages
  count_cols <- combined_counts %>% 
    select(-c(Data_ID:Disease_group))
  
  total_cells <- rowSums(count_cols)
  
  percentages <- map_df(count_cols, ~(. / total_cells) * 100) %>%
    rename_with(~paste0("PCT_", .), everything())
  
  bind_cols(combined_counts, percentages)
}

# 7. Main Execution -----------------------------------------------
main <- function() {
  # Create output directory
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Process samples
  seurat_list <- map(SAMPLE_METADATA$patient_id, process_sample)
  names(seurat_list) <- SAMPLE_METADATA$patient_id
  
  # Integrate samples
  integrated_data <- integrate_samples(seurat_list)
  
  # Annotate cell types
  annotated_data <- annotate_cell_types(integrated_data)
  
  # Create cell type subsets
  cell_subsets <- subset_cell_types(annotated_data)
  
  # Save results
  saveRDS(annotated_data, file.path(OUTPUT_DIR, "complete_analysis.rds"))
  
  # Save cell type subsets
  walk2(cell_subsets, names(cell_subsets), function(subset, name) {
    saveRDS(subset, file.path(OUTPUT_DIR, 
            paste0("subset_", tolower(name), ".rds")))
  })
  
  # Analyze and save cell counts
  cell_counts <- analyze_cell_counts(cell_subsets)
  write_tsv(cell_counts, file.path(OUTPUT_DIR, "cell_type_counts.tsv"))
  
  message("Analysis pipeline completed successfully!")
}

# Run pipeline
if (!interactive()) {
  main()
}
