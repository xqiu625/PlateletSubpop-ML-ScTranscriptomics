#!/usr/bin/env Rscript

#' Gene Name Processing Pipeline
#' 
#' This script processes and standardizes gene names from a Seurat object using HGNChelper.
#' It handles special cases like AB- prefixes and dot notations, and resolves multiple symbols.

# 1. Setup and Configuration -------------------------------------------

# Load required libraries
suppressPackageMessages({
  library(Seurat)      # For single-cell analysis
  library(HGNChelper)  # For gene name standardization
  library(tidyverse)   # For data manipulation
  library(stringr)     # For string operations
})

# Set up paths
INPUT_DIR <- "/bigdata/godziklab/shared/Xinru/302005/datasets"
OUTPUT_DIR <- file.path(INPUT_DIR, "processed_datasets")
R_LIB_PATH <- "/bigdata/godziklab/shared/Xinru/R"

# 2. Helper Functions -------------------------------------------------

#' Process Gene Names with Special Prefixes
#' @param gene_data DataFrame with gene names
#' @param prefix Pattern to match (e.g., "^AB-")
#' @return DataFrame with processed gene names
process_prefix_genes <- function(gene_data, prefix) {
  # Extract genes with prefix
  prefix_genes <- gene_data %>% 
    filter(grepl(prefix, x)) %>%
    mutate(y = gsub(prefix, "", x))
  
  # Check symbols without prefix
  checked_symbols <- checkGeneSymbols(prefix_genes$y, species = "human") %>%
    as_tibble() %>%
    filter(!is.na(Suggested.Symbol)) %>%
    rename(y = x)
  
  # Merge and format results
  prefix_genes %>%
    inner_join(checked_symbols, by = "y") %>%
    select(x, Suggested.Symbol) %>%
    distinct(x, .keep_all = TRUE)
}

#' Process Gene Names with Dots
#' @param gene_data DataFrame with gene names
#' @return DataFrame with processed gene names
process_dot_genes <- function(gene_data) {
  # Extract genes with dots
  dot_genes <- gene_data %>%
    filter(grepl("\\..", x)) %>%
    mutate(y = gsub("\\..*", "", x))
  
  # Check symbols without dots
  checked_symbols <- checkGeneSymbols(dot_genes$y, species = "human") %>%
    as_tibble() %>%
    filter(!is.na(Suggested.Symbol)) %>%
    rename(y = x)
  
  # Merge and format results
  dot_genes %>%
    inner_join(checked_symbols, by = "y") %>%
    select(x, Suggested.Symbol.y) %>%
    distinct(x, .keep_all = TRUE) %>%
    rename(Suggested.Symbol = Suggested.Symbol.y)
}

#' Clean Multiple Symbols
#' @param gene_data DataFrame with gene names
#' @return DataFrame with cleaned gene names
clean_multiple_symbols <- function(gene_data) {
  # Split data into single and multiple symbols
  multiple_symbols <- gene_data %>% 
    filter(grepl("///", Suggested.Symbol)) %>%
    mutate(
      Suggested.Symbol = str_split_fixed(Suggested.Symbol, "///", 2)[,1],
      Suggested.Symbol = trimws(Suggested.Symbol)
    )
  
  single_symbols <- gene_data %>% 
    filter(!grepl("///", Suggested.Symbol))
  
  # Combine results
  bind_rows(multiple_symbols, single_symbols)
}

# 3. Main Processing Function -----------------------------------------

#' Process Gene Names from Seurat Object
#' @param seurat_obj Seurat object containing gene expression data
#' @return DataFrame with processed gene names
process_gene_names <- function(seurat_obj) {
  # Set RNA as default assay
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Extract gene list
  gene_list <- data.frame(Gene = rownames(seurat_obj@assays$RNA@counts))
  
  # Initial HGNC check
  initial_check <- checkGeneSymbols(gene_list$Gene, species = "human") %>%
    as_tibble()
  
  # Process standard gene names
  standard_genes <- initial_check %>%
    filter(!is.na(Suggested.Symbol)) %>%
    select(x, Suggested.Symbol)
  
  # Process special cases
  ab_genes <- process_prefix_genes(initial_check, "^AB-")
  dot_genes <- process_dot_genes(initial_check)
  
  # Combine all results
  combined_genes <- bind_rows(
    standard_genes,
    ab_genes,
    dot_genes
  )
  
  # Clean multiple symbols
  clean_multiple_symbols(combined_genes)
}

# 4. Main Execution --------------------------------------------------

main <- function() {
  # Create output directory if needed
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  message("Reading Seurat object...")
  seurat_data <- readRDS(file.path(INPUT_DIR, 
    "platelets_seurat_integrated_SingleRfiltered_clustered_findallmarkers.rds"))
  
  # Process gene names
  message("Processing gene names...")
  corrected_genes <- process_gene_names(seurat_data)
  
  # Save results
  message("Saving processed gene names...")
  output_file <- file.path(OUTPUT_DIR, "Corrected_gene_name.rds")
  saveRDS(corrected_genes, file = output_file)
  
  message("Gene name processing completed successfully!")
  message("Results saved to: ", output_file)
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}
