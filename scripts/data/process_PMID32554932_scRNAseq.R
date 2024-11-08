#!/usr/bin/env Rscript

#' Single-cell RNA-seq Processing Pipeline for ARDS and Sepsis Dataset
#' Dataset: Jiang et al., 2020 (PMID: 32554932)
#' Description: Analysis of ARDS and sepsis patient immune responses
#' GEO: GSE151263
#' 
#' @author [Your Name]
#' @date 2024-11-07

# 1. Setup and Configuration -------------------------------------------
suppressPackageMessages({
  library(Seurat)
  library(tidyverse)
  library(HGNChelper)
  library(MAST)
  library(scCATCH)
  library(SingleR)
  library(celldex)
  library(scater)
})

# Set paths
BASE_DIR <- "/datasets"
DATA_DIR <- file.path(BASE_DIR, "PMID32554932")
OUTPUT_DIR <- file.path(BASE_DIR, "PMID32554932")
R_LIB_PATH <- "/R"

.libPaths(c(R_LIB_PATH, .libPaths()))

# Sample definitions
SAMPLES <- c("ARDS1", "ARDS2", "ARDS3", "Sepsis1", "Sepsis2", "Sepsis3", "Sepsis4")

# 2. Data Processing Functions ----------------------------------------
#' Process Gene Names and Create Seurat Object
#' @param file_path Path to data file
#' @param sample_name Sample identifier
#' @return Processed Seurat object
process_sample <- function(file_path, sample_name) {
  message(sprintf("Processing %s...", sample_name))
  
  # Read data
  df <- read.delim2(file_path) %>%
    mutate(X = trimws(X)) %>%
    filter(X != '<NA>', !duplicated(X))
  
  # Process gene names
  gene_check <- checkGeneSymbols(df$X, species = "human")
  processed_symbols <- gene_check %>%
    filter(!is.na(Suggested.Symbol), Approved == 'FALSE') %>%
    mutate(
      Symbol = str_split_fixed(Suggested.Symbol, "///", 2)[,1] %>% trimws(),
      X = x
    ) %>%
    select(X, Symbol) %>%
    distinct(X, .keep_all = TRUE) %>%
    filter(!Symbol %in% c(df$X, .$Symbol[duplicated(.$Symbol)]))
  
  # Update gene names
  df <- left_join(df, processed_symbols, by = "X") %>%
    mutate(X = ifelse(!is.na(Symbol), Symbol, X))
  
  rownames(df) <- df$X
  df <- df %>% select(-X, -Symbol)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = df, project = sample_name, 
                                  min.cells = 3, min.features = 200)
  
  # Add QC metrics
  seurat_obj <- seurat_obj %>%
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
    PercentageFeatureSet(pattern = "^RPS", col.name = "percent.rps") %>%
    PercentageFeatureSet(pattern = "^RPL", col.name = "percent.rpl") %>%
    PercentageFeatureSet(pattern = "^RNA\\d8S5", col.name = "percent.rrna")
  
  return(seurat_obj)
}

#' Integrate Multiple Samples
#' @param sample_list List of processed Seurat objects
#' @return Integrated Seurat object
integrate_samples <- function(sample_list) {
  # SCTransform each dataset
  sample_list <- map(sample_list, function(x) {
    NormalizeData(x) %>%
      SCTransform(vars.to.regress = c("percent.mt"))
  })
  
  # Integration
  features <- SelectIntegrationFeatures(sample_list, nfeatures = 3000)
  sample_list <- PrepSCTIntegration(sample_list, anchor.features = features)
  
  anchors <- FindIntegrationAnchors(sample_list,
                                   normalization.method = "SCT",
                                   anchor.features = features)
  
  integrated <- IntegrateData(anchors, normalization.method = "SCT")
  
  return(integrated)
}

#' Process Integrated Data
#' @param seurat_obj Integrated Seurat object
#' @return Processed Seurat object
process_integrated_data <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  
  seurat_obj <- seurat_obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", 
                        nfeatures = 3000, 
                        verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:50, verbose = FALSE) %>%
    FindNeighbors(dims = 1:50, verbose = FALSE) %>%
    FindClusters(resolution = 1, verbose = FALSE)
  
  # Find markers
  seurat_obj@misc$markers <- FindAllMarkers(seurat_obj, 
                                           assay = 'RNA',
                                           only.pos = TRUE, 
                                           test.use = 'MAST')
  
  return(seurat_obj)
}

# 3. Cell Type Annotation -------------------------------------------
#' Annotate Cell Types Using Multiple Methods
#' @param seurat_obj Processed Seurat object
#' @return Annotated Seurat object
annotate_cell_types <- function(seurat_obj) {
  # scCATCH annotation
  seurat_obj@misc$scCATCH_cellType <- scCATCH(
    object = seurat_obj@misc$markers,
    species = 'Human',
    cancer = NULL,
    tissue = 'Blood'
  )
  
  # SingleR annotation
  ref_data <- readRDS(file.path(BASE_DIR, "celldex_db/HumanPrimaryCellAtlasData.rds"))
  
  # Prepare data for SingleR
  common_genes <- intersect(rownames(ref_data), rownames(seurat_obj))
  seurat_subset <- DietSeurat(seurat_obj[common_genes,], graphs = "umap")
  sce_obj <- as.SingleCellExperiment(seurat_subset)
  sce_obj <- logNormCounts(sce_obj)
  
  # Run SingleR
  pred.hpca <- SingleR(sce_obj, ref = ref_data[common_genes,],
                       labels = ref_data$label.main)
  
  # Process results
  singler_results <- data.frame(
    seurat_clusters = seurat_obj$seurat_clusters,
    pruned.labels = pred.hpca$pruned.labels
  ) %>%
    group_by(seurat_clusters) %>%
    count(pruned.labels, sort = TRUE) %>%
    slice_max(n = 1, n)
  
  names(singler_results)[1] <- "cluster"
  seurat_obj@misc$SingleR_cellType <- singler_results
  
  return(seurat_obj)
}

# 4. Cell Type Subsetting ------------------------------------------
CELL_CLUSTERS <- list(
  Platelet = '18',
  Monocytes = c('0', '2', '4', '5', '15', '17', '20', '24', '25', '27', '14'),
  T_cell = c('1', '3', '7', '8', '11', '13', '22', '12', '10', '16'),
  B_cell = c('6', '19', '21'),
  Other = c('9', '23', '26')
)

#' Subset Cell Types and Calculate Counts
#' @param seurat_obj Annotated Seurat object
#' @return List of cell type subsets and counts
subset_and_count_cells <- function(seurat_obj) {
  # Create subsets
  subsets <- map(names(CELL_CLUSTERS), function(type) {
    subset_obj <- subset(seurat_obj, idents = CELL_CLUSTERS[[type]])
    subset_obj@meta.data <- subset_obj@meta.data %>%
      mutate(
        Cell_type = type,
        Disease = "Sepsis",
        Data_ID = "Jiang et al.",
        Disease_group = paste("mechanically ventilated", "icu",
                            orig.ident, sep = "_")
      )
    return(subset_obj)
  }) %>%
    set_names(names(CELL_CLUSTERS))
  
  # Calculate counts
  counts <- map(subsets, function(subset) {
    FetchData(subset, vars = c("Data_ID", "orig.ident", 
                              "Disease", "Disease_group", "Cell_type")) %>%
      group_by(Data_ID, orig.ident, Disease, Disease_group, Cell_type) %>%
      count() %>%
      spread(Cell_type, n)
  }) %>%
    reduce(full_join, by = c("Data_ID", "orig.ident", 
                            "Disease", "Disease_group"))
  
  # Process counts
  counts[is.na(counts)] <- 0
  counts <- rename(counts, Sample_ID = orig.ident)
  
  # Calculate percentages
  total_cells <- rowSums(select(counts, -c(Data_ID:Disease_group)))
  percentages <- map_df(select(counts, -c(Data_ID:Disease_group)), 
                       ~(. / total_cells) * 100) %>%
    rename_with(~paste0("PCT_", .), everything())
  
  # Combine results
  final_counts <- bind_cols(
    counts,
    percentages,
    tibble(
      Erythroid = NA_real_,
      PCT_Erythroid = NA_real_
    )
  )
  
  return(list(
    subsets = subsets,
    counts = final_counts
  ))
}

# 5. Main Execution -----------------------------------------------
main <- function() {
  # Create output directory
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Process individual samples
  message("Processing samples...")
  sample_data <- map(SAMPLES, function(sample) {
    file_path <- file.path(DATA_DIR, paste0("GSM4569780_", sample, "_processed_UMI.txt"))
    process_sample(file_path, sample)
  })
  names(sample_data) <- SAMPLES
  
  # Integrate samples
  message("Integrating samples...")
  integrated_data <- integrate_samples(sample_data)
  
  # Process integrated data
  message("Processing integrated data...")
  processed_data <- process_integrated_data(integrated_data)
  
  # Annotate cell types
  message("Annotating cell types...")
  annotated_data <- annotate_cell_types(processed_data)
  
  # Subset and count cells
  message("Creating cell type subsets and calculating counts...")
  results <- subset_and_count_cells(annotated_data)
  
  # Save results
  message("Saving results...")
  saveRDS(annotated_data, file.path(OUTPUT_DIR, "complete_analysis.rds"))
  
  walk2(results$subsets, names(results$subsets), function(subset, name) {
    saveRDS(subset, file.path(OUTPUT_DIR, paste0("subset_", tolower(name), ".rds")))
  })
  
  write_tsv(results$counts, file.path(OUTPUT_DIR, "cell_type_counts.tsv"))
  
  message("Analysis pipeline completed successfully!")
}

# Run pipeline
if (!interactive()) {
  main()
}
