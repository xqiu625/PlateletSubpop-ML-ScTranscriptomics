#!/usr/bin/env Rscript

#' Single-cell RNA-seq Processing Pipeline for COVID-19 PBMC Dataset
#' Dataset: Wilk et al., 2020 (PMID: 32514174)
#' Description: Process and analyze single-cell RNA-seq data from COVID-19 patients
#' https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150728

# 1. Setup and Configuration -------------------------------------------
# Load required libraries
suppressPackageMessages({
  library(Seurat)
  library(MAST)
  library(scCATCH)
  library(SingleCellExperiment)
  library(celldex)
  library(SingleR)
  library(scater)
  library(tidyverse)
})

# Set paths
BASE_DIR <- "/bigdata/godziklab/shared/Xinru/302005/datasets"
OUTPUT_DIR <- file.path(BASE_DIR, "PMID32514174")
R_LIB_PATH <- "/bigdata/godziklab/shared/Xinru/R"

# Add custom R library path
.libPaths(c(R_LIB_PATH, .libPaths()))

# 2. Data Loading and Initial Processing -------------------------------
process_initial_data <- function(path) {
  # Read and process matrices
  cm.list <- list.files(pattern = "*.matrices.rds", path = path, full.names = TRUE)
  cm.files <- lapply(cm.list, readRDS)
  names(cm.files) <- sub(path, "", sub("\\_cell.counts.matrices.rds", "", cm.list))
  
  # Process and merge
  cm.pp <- mapply(EpicPreHS, cm.files, orig.ident = names(cm.files), SIMPLIFY = FALSE)
  covid_combined.emat <- mergeCM(cm.pp, type = "emat")
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = covid_combined.emat, min.cells = 10, 
                                  names.field = 1, names.delim = "\\.")
  
  return(seurat_obj)
}

# 3. Quality Control and Preprocessing --------------------------------
perform_qc <- function(seurat_obj) {
  seurat_obj <- seurat_obj %>%
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
    PercentageFeatureSet(pattern = "^RPS", col.name = "percent.rps") %>%
    PercentageFeatureSet(pattern = "^RPL", col.name = "percent.rpl") %>%
    PercentageFeatureSet(pattern = "^RNA\\d8S5", col.name = "percent.rrna")
  
  # SCTransform and dimensionality reduction
  seurat_obj <- seurat_obj %>%
    SCTransform(vars.to.regress = c("percent.mt", "percent.rps", "percent.rpl", 
                                   "percent.rrna", "nCount_RNA", "nFeature_RNA"), 
                verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:50, verbose = FALSE) %>%
    FindNeighbors(dims = 1:50, verbose = FALSE) %>%
    FindClusters(resolution = 1, verbose = FALSE)
  
  return(seurat_obj)
}

# 4. Cell Type Annotation --------------------------------------------
annotate_cell_types <- function(seurat_obj) {
  # Find markers
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj@misc$markers <- FindAllMarkers(seurat_obj, assay = 'RNA', 
                                           only.pos = TRUE, test.use = 'MAST')
  
  # scCATCH annotation
  clu_ann <- scCATCH(object = seurat_obj@misc$markers,
                     species = 'Human',
                     cancer = NULL,
                     tissue = 'Blood')
  seurat_obj@misc$scCATCH_cellType <- clu_ann
  
  # SingleR annotation
  ref.data <- readRDS(file.path(BASE_DIR, "celldex_db/HumanPrimaryCellAtlasData.rds"))
  common <- intersect(rownames(ref.data), rownames(seurat_obj))
  seurat_singler <- DietSeurat(seurat_obj[common,], graphs = "umap")
  sce_obj <- as.SingleCellExperiment(seurat_singler)
  sce_obj <- logNormCounts(sce_obj)
  
  pred.hpca <- SingleR(sce_obj, ref = ref.data[common,], 
                       labels = ref.data$label.main)
  
  # Process SingleR results
  SingleR_label <- data.frame(
    seurat_clusters = seurat_obj$seurat_clusters,
    pruned.labels = pred.hpca$pruned.labels
  ) %>%
    group_by(seurat_clusters) %>%
    count(pruned.labels, sort = TRUE) %>%
    slice_max(n = 1, n)
  
  seurat_obj@misc$SingleR_cellType <- SingleR_label
  
  return(seurat_obj)
}

# 5. Cell Type Subsetting -------------------------------------------
subset_cell_types <- function(seurat_obj) {
  cell_types <- list(
    Platelet = '19',
    Monocytes = c('3','5','6','7','9','11','21','22','28'),
    Erythroid = '15',
    T_cell = c('0', '1', '2', '10', '12', '17'),
    B_cell = c('4','8','14', '16', '18', '20', '27', '29'),
    Other = c('13', '23', '24', '25', '26')
  )
  
  subsets <- lapply(names(cell_types), function(type) {
    subset_obj <- subset(seurat_obj, idents = cell_types[[type]])
    subset_obj@meta.data$Cell_type <- type
    return(subset_obj)
  })
  
  names(subsets) <- names(cell_types)
  return(subsets)
}

# 6. Cell Count Analysis -------------------------------------------
analyze_cell_counts <- function(cell_subsets) {
  cell_counts <- lapply(names(cell_subsets), function(type) {
    FetchData(cell_subsets[[type]], 
             vars = c("Data_NO", "Data_ID", "Donor.full", "Disease", 
                     "Disease_group", "Cell_type")) %>%
      group_by(Data_NO, Data_ID, Donor.full, Disease, Disease_group, Cell_type) %>%
      count() %>%
      spread(Cell_type, n)
  })
  
  # Combine all cell counts
  combined_counts <- reduce(cell_counts, full_join, by = c("Data_NO", "Data_ID", 
                                                          "Donor.full", "Disease", 
                                                          "Disease_group"))
  
  # Calculate percentages
  combined_counts[is.na(combined_counts)] <- 0
  total_cells <- rowSums(combined_counts[, -(1:5)])
  percentages <- combined_counts[, -(1:5)] / total_cells
  names(percentages) <- paste0("PCT_", names(percentages))
  
  final_counts <- bind_cols(combined_counts, percentages)
  return(final_counts)
}

# 7. Main Execution ------------------------------------------------
main <- function() {
  # Create output directory
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Process data
  seurat_obj <- process_initial_data(file.path(BASE_DIR, "302004data01/COVID-19/"))
  seurat_obj <- perform_qc(seurat_obj)
  seurat_obj <- annotate_cell_types(seurat_obj)
  
  # Save processed Seurat object
  saveRDS(seurat_obj, file.path(OUTPUT_DIR, "processed_seurat_object.rds"))
  
  # Generate cell type subsets
  cell_subsets <- subset_cell_types(seurat_obj)
  
  # Save individual cell type objects
  for (type in names(cell_subsets)) {
    saveRDS(cell_subsets[[type]], 
            file.path(OUTPUT_DIR, paste0("subset_", tolower(type), ".rds")))
  }
  
  # Analyze cell counts
  cell_counts <- analyze_cell_counts(cell_subsets)
  write_tsv(cell_counts, file.path(OUTPUT_DIR, "cell_type_counts.tsv"))
  
  message("Analysis pipeline completed successfully!")
}

# Run pipeline
if (!interactive()) {
  main()
}
