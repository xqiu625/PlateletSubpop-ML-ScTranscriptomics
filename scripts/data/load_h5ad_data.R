#!/usr/bin/env Rscript

#' Load H5AD Data and Convert to Seurat Object
#' Description: Read AnnData h5ad file and convert to Seurat for ML pipeline

# 1. Setup -------------------------------------------------------------
suppressPackageMessages({
  library(Seurat)
  library(SeuratDisk)
  library(tidyverse)
})

#' Convert H5AD to Seurat Object
#' @param h5ad_path Path to .h5ad file
#' @param output_path Optional path to save .rds file
#' @return Seurat object
convert_h5ad_to_seurat <- function(h5ad_path, output_path = NULL) {

  message(sprintf("Loading h5ad file: %s", h5ad_path))

  # Method 1: Using SeuratDisk (preferred)
  # First convert h5ad to h5seurat
  h5seurat_path <- gsub("\\.h5ad$", ".h5seurat", h5ad_path)

  if (!file.exists(h5seurat_path)) {
    message("Converting h5ad to h5seurat format...")
    Convert(h5ad_path, dest = "h5seurat", overwrite = TRUE)
  }

  # Load h5seurat
  message("Loading h5seurat file...")
  seurat_obj <- LoadH5Seurat(h5seurat_path)

  # Print info
  message(sprintf("Loaded Seurat object:"))
  message(sprintf("  - Cells: %d", ncol(seurat_obj)))
  message(sprintf("  - Genes: %d", nrow(seurat_obj)))
  message(sprintf("  - Assays: %s", paste(Assays(seurat_obj), collapse = ", ")))

  # Check metadata columns
  message(sprintf("  - Metadata columns: %s",
                  paste(head(colnames(seurat_obj@meta.data), 10), collapse = ", ")))

  # Save if output path provided
  if (!is.null(output_path)) {
    message(sprintf("Saving to: %s", output_path))
    saveRDS(seurat_obj, output_path)
  }

  return(seurat_obj)
}

#' Quick Data Summary
#' @param seurat_obj Seurat object
print_data_summary <- function(seurat_obj) {
  cat("\n=== Data Summary ===\n")
  cat(sprintf("Cells: %d\n", ncol(seurat_obj)))
  cat(sprintf("Genes: %d\n", nrow(seurat_obj)))

  cat("\n=== Metadata Columns ===\n")
  print(colnames(seurat_obj@meta.data))

  # Check for common columns
  common_cols <- c("Outcome", "Severity", "Disease", "Cluster",
                   "seurat_clusters", "Cell_type")
  available <- common_cols[common_cols %in% colnames(seurat_obj@meta.data)]

  if (length(available) > 0) {
    cat("\n=== Available Classification Columns ===\n")
    for (col in available) {
      cat(sprintf("\n%s:\n", col))
      print(table(seurat_obj[[col, drop = TRUE]]))
    }
  }
}

# 2. Main Execution ----------------------------------------------------
main <- function() {
  # Set paths for HPCC
  BASE_DIR <- "/bigdata/godziklab/shared/Xinru/302005"
  H5AD_FILE <- file.path(BASE_DIR, "302005_platelet_harmony_integrated.h5ad")
  OUTPUT_FILE <- file.path(BASE_DIR, "302005_platelet_harmony_integrated.rds")

  # Check if file exists
  if (!file.exists(H5AD_FILE)) {
    stop(sprintf("File not found: %s", H5AD_FILE))
  }

  # Convert and load
  seurat_obj <- convert_h5ad_to_seurat(H5AD_FILE, OUTPUT_FILE)

  # Print summary
  print_data_summary(seurat_obj)

  message("\nData loading completed!")
  message(sprintf("Seurat object saved to: %s", OUTPUT_FILE))
}

# Run
if (!interactive()) {
  main()
}
