#!/usr/bin/env Rscript

#' Description: Process and using Seurat v4 to label cell types
#' 
#' Datasets:
#' - PMID: 32514174 (Wilk et al., 2020)
#' - PMID: 32651212 (Lee et al., 2020)
#' - PMID: 32810438 (Schulte-Schrepping et al., 2020)
#' - PMID: 32788292 (Arunachalam et al., 2020)
#' - PMID: 33657410 (Ren et al., 2021)
#' - PMID: 33713619 (Liu et al., 2021)
#' - PMID: 33879890 (Stephenson et al., 2021)
#' - PMID: 32554932 (Jiang et al., 2020)
#' - PMID: 33296687 (Bernardes et al., 2020)
#' - PMID: 34558746 (Qiu et al., 2021)
#' - PMID: 33494096 (Combes et al., 2021)
#' - PMID: 31754025 (Mistry et al., 2019)

# 1. Setup and Configuration -------------------------------------------
suppressPackageMessages({
  library(Seurat)
  library(SeuratDisk)
  library(tidyverse)
  library(HGNChelper)
  library(SingleR)
  library(celldex)
})

# Set paths
BASE_DIR <- "/datasets"
OUTPUT_DIR <- "/datasets_processed"
R_LIB_PATH <- "/R"

# 2. Processing Functions ---------------------------------------------
#' Load and Preprocess Dataset
#' @param pmid PMID identifier
#' @return Loaded Seurat object
load_dataset <- function(pmid) {
  old_id <- DATASET_MAPPING[[pmid]]
  file_path <- file.path(BASE_DIR, old_id, paste0(old_id, "_covid19.rds"))
  
  # Special case for h5seurat file
  if (pmid == "33879890") {
    file_path <- file.path(BASE_DIR, old_id, "covid_portal_210320_with_raw.h5seurat")
    obj <- LoadH5Seurat(file_path)
    DefaultAssay(obj) <- "raw"
    return(obj)
  }
  
  readRDS(file_path)
}

#' Process Standard Dataset
#' @param obj Seurat object
#' @param gene_name Gene name mapping data
#' @param reference Reference dataset
#' @return Processed Seurat object
process_standard <- function(obj, gene_name, reference) {
  # Extract metadata and counts
  meta <- obj@meta.data
  mat <- obj@assays$RNA@counts
  
  # Process gene names
  dat <- tibble(x = rownames(mat)) %>%
    left_join(gene_name, by = "x") %>%
    filter(!is.na(Suggested.Symbol))
  
  # Update matrix
  mat <- mat[dat$x, ]
  rownames(mat) <- dat$Suggested.Symbol
  
  # Create new object and process
  CreateSeuratObject(counts = mat) %>%
    `@`(meta.data = meta) %>%
    SCTransform(verbose = FALSE) %>%
    process_with_reference(reference)
}

#' Process H5 Dataset
#' @param obj Seurat object
#' @param gene_name Gene name mapping data
#' @param reference Reference dataset
#' @return Processed Seurat object
process_h5 <- function(obj, gene_name, reference) {
  # Extract metadata and counts
  meta <- obj@meta.data
  mat <- obj@assays$raw@counts
  
  # Process gene names
  dat <- tibble(x = rownames(mat)) %>%
    left_join(gene_name, by = "x") %>%
    filter(!is.na(Suggested.Symbol))
  
  # Update matrix
  mat <- mat[dat$x, ]
  rownames(mat) <- dat$Suggested.Symbol
  
  # Create new object and process
  CreateSeuratObject(counts = mat) %>%
    `@`(meta.data = meta) %>%
    SCTransform(verbose = FALSE) %>%
    process_with_reference(reference)
}

#' Process Dataset with Reference
#' @param obj Processed Seurat object
#' @param reference Reference dataset
#' @return Annotated Seurat object
process_with_reference <- function(obj, reference) {
  # Find transfer anchors
  anchors <- FindTransferAnchors(
    reference = reference,
    query = obj,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50,
    recompute.residuals = FALSE
  )
  
  # Map query
  MapQuery(
    anchorset = anchors,
    query = obj,
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
}

#' Process Cell Type Subsets
#' @param obj Seurat object with SingleR labels
#' @param cell_type Cell type to subset
#' @return Processed subset
process_cell_subset <- function(obj, cell_type) {
  subset(obj, pruned.labels == cell_type) %>%
    process_standard(gene_name, reference_pbmc)
}

# 3. Main Execution ------------------------------------------------
main <- function() {
  # Load reference data
  message("Loading reference data...")
  reference_pbmc <- LoadH5Seurat("/seurat_v4/pbmc_multimodal.h5seurat")
  
  # Load gene name mappings
  message("Loading gene name mappings...")
  gene_name <- list.files(pattern = "*.txt", 
                         path = file.path(OUTPUT_DIR, "correct_gene_name"), 
                         full.names = TRUE) %>%
    map_df(read.delim2) %>%
    distinct(x, .keep_all = TRUE)
  
  # Process each dataset
  for (pmid in names(DATASET_MAPPING)) {
    message(sprintf("Processing dataset PMID:%s...", pmid))
    
    # Load and process dataset
    obj <- load_dataset(pmid)
    processed <- if (pmid == "33879890") {
      process_h5(obj, gene_name, reference_pbmc)
    } else {
      process_standard(obj, gene_name, reference_pbmc)
    }
    
    # Save metadata
    write_tsv(processed@meta.data,
              file.path(OUTPUT_DIR, "processed_datasets_v3", 
                       sprintf("PMID%s_annotation_SingleR.tsv", pmid)))
    
    # Process T cell and Monocyte subsets if needed
    if (pmid %in% c("33657410", "33879890")) {
      for (cell_type in c("T_cells", "Monocyte")) {
        subset_processed <- process_cell_subset(obj, cell_type)
        write_tsv(subset_processed@meta.data,
                 file.path(OUTPUT_DIR, "processed_datasets_v3",
                          sprintf("PMID%s_%s_annotation_SeuratV4.tsv", 
                                pmid, gsub("_cells", "", cell_type))))
      }
    }
  }
  
  message("Processing completed successfully!")
}

# Run pipeline
if (!interactive()) {
  main()
}
