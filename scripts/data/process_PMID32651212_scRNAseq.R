#!/usr/bin/env Rscript

#' Single-cell RNA-seq Processing Pipeline for COVID-19 and Influenza Dataset
#' Dataset: Lee et al., 2020 (PMID: 32651212)
#' Description: Process and analyze single-cell RNA-seq data comparing COVID-19 and influenza patients
#' GSE149689: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689

# 1. Setup and Configuration -------------------------------------------
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
OUTPUT_DIR <- file.path(BASE_DIR, "PMID32651212")
R_LIB_PATH <- "/bigdata/godziklab/shared/Xinru/R"
REFERENCE_DIR <- file.path(BASE_DIR, "celldex_db")

# Add custom R library path
.libPaths(c(R_LIB_PATH, .libPaths()))

# 2. Quality Control and Processing -----------------------------------
process_seurat_object <- function(seurat_obj) {
  seurat_obj <- seurat_obj %>%
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
    PercentageFeatureSet(pattern = "^RPS", col.name = "percent.rps") %>%
    PercentageFeatureSet(pattern = "^RPL", col.name = "percent.rpl") %>%
    PercentageFeatureSet(pattern = "^RNA\\d8S5", col.name = "percent.rrna")
  
  # Normalization and dimensionality reduction
  seurat_obj <- seurat_obj %>%
    SCTransform(vars.to.regress = c("percent.mt", "percent.rps", "percent.rpl", 
                                   "percent.rrna", "nCount_RNA", "nFeature_RNA"), 
                verbose = FALSE) %>%
    NormalizeData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:50, verbose = FALSE) %>%
    FindNeighbors(dims = 1:50, verbose = FALSE) %>%
    FindClusters(resolution = 1, verbose = FALSE)
  
  return(seurat_obj)
}

# 3. Cell Type Annotation -------------------------------------------
annotate_cell_types <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Find markers
  seurat_obj@misc$markers <- FindAllMarkers(seurat_obj, assay = 'RNA', 
                                           only.pos = TRUE, test.use = 'MAST')
  
  # Canonical marker annotation
  canonical_markers <- annotate_canonical_markers(seurat_obj@misc$markers)
  seurat_obj@misc$Canonical_cellType <- canonical_markers
  
  # scCATCH annotation
  clu_ann <- scCATCH(object = seurat_obj@misc$markers,
                     species = 'Human',
                     cancer = NULL,
                     tissue = 'Blood')
  seurat_obj@misc$scCATCH_cellType <- clu_ann
  
  # SingleR annotation
  singler_results <- perform_singler_annotation(seurat_obj)
  seurat_obj@misc$SingleR_cellType <- singler_results
  
  return(seurat_obj)
}

# Helper function for canonical marker annotation
annotate_canonical_markers <- function(marker_data) {
  markers <- marker_data %>% 
    select(cluster, gene) %>%
    group_by(cluster) %>%
    summarize(genes = paste(gene, collapse = "\", \"")) %>%
    mutate(genes = paste0("\"", genes, "\""))
  
  # Define cell type markers
  cell_types <- tibble(
    marker_pairs = list(
      c("CD14", "LYZ") ~ "CD14+ Mono",
      "MS4A1" ~ "B",
      c("IL7R", "CCR7") ~ "CD4+ T",
      c("IL7R", "CD27") ~ "CD4+ T",
      "CD8A" ~ "CD8+ T",
      c("FCER1A", "CST3") ~ "DC",
      c("CD123", "GZMB") ~ "DC",
      c("GYPB", "AHSP") ~ "Eryth precursor",
      c("FCGR3A", "MS4A7") ~ "FCGR3A+ Mono",
      c("GNLY", "NKG7") ~ "NK",
      "PPBP" ~ "Platelet"
    )
  ) %>%
    unnest(marker_pairs)
  
  # Apply cell type annotations
  markers <- markers %>%
    mutate(cell_type = apply_cell_type_rules(genes, cell_types))
  
  return(select(markers, cluster, cell_type))
}

# 4. Disease Group Assignment ----------------------------------------
assign_disease_groups <- function(metadata) {
  disease_groups <- tribble(
    ~pattern, ~group,
    "nCoV [1,3,4,7,8,9]", "severe COVID-19",
    "nCoV [2,5,6,10,11]", "mild COVID-19",
    "Flu [1-5]", "severe influenza",
    "Normal_[1-4]", "healthy donor"
  )
  
  metadata %>%
    mutate(Disease_group = case_when(
      str_detect(cells, "nCoV [1,3,4,7,8,9]") ~ "severe COVID-19",
      str_detect(cells, "nCoV [2,5,6,10,11]") ~ "mild COVID-19",
      str_detect(cells, "Flu [1-5]") ~ "severe influenza",
      str_detect(cells, "Normal_[1-4]") ~ "healthy donor",
      TRUE ~ NA_character_
    ))
}

# 5. Cell Type Subsetting ------------------------------------------
define_cell_subsets <- function(seurat_obj) {
  cell_types <- list(
    Platelet = '17',
    Monocytes = c('0', '7', '10', '11', '13', '14', '15', '19', '23', '25', '26', '28', '29'),
    Erythroid = '20',
    T_cell = c('5', '6', '2', '1', '16', '21', '24', '27'),
    B_cell = c('4', '9', '22', '32'),
    Other = c('3', '8', '12', '18', '30', '31', '33')
  )
  
  subsets <- map(names(cell_types), function(type) {
    subset(seurat_obj, idents = cell_types[[type]]) %>%
      AddMetaData(metadata = type, col.name = "Cell_type")
  })
  
  names(subsets) <- names(cell_types)
  return(subsets)
}

# 6. Cell Count Analysis ------------------------------------------
analyze_cell_counts <- function(cell_subsets) {
  # Process each cell type
  cell_counts <- map(names(cell_subsets), function(type) {
    FetchData(cell_subsets[[type]], 
             vars = c("Data_NO", "Data_ID", "sample", "Disease", 
                     "Disease_group", "Cell_type")) %>%
      group_by(Data_NO, Data_ID, sample, Disease, Disease_group, Cell_type) %>%
      count() %>%
      spread(Cell_type, n)
  })
  
  # Combine counts and calculate percentages
  counts_combined <- reduce(cell_counts, full_join, 
                          by = c("Data_NO", "Data_ID", "sample", 
                                "Disease", "Disease_group")) %>%
    replace(is.na(.), 0)
  
  total_cells <- select(counts_combined, -c(Data_NO:Disease_group)) %>%
    rowSums()
  
  percentages <- counts_combined %>%
    select(-c(Data_NO:Disease_group)) %>%
    map_df(~(. / total_cells) * 100) %>%
    rename_with(~paste0("PCT_", .), everything())
  
  bind_cols(counts_combined, percentages)
}

# 7. Visualization -----------------------------------------------
create_umap_plot <- function(seurat_obj, output_file) {
  png(file = output_file,
      width = 12 * 300,
      height = 12 * 300,
      units = "px",
      res = 300,
      type = 'cairo')
  
  print(DimPlot(seurat_obj,
                reduction = "umap",
                label = TRUE,
                label.size = 12,
                pt.size = 0.01) + 
          NoLegend())
  
  dev.off()
}

# 8. Main Execution ---------------------------------------------
main <- function() {
  # Create output directory
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  message("Loading Seurat object...")
  seurat_obj <- readRDS(file.path(BASE_DIR, "filtered_PMID32651212.rds"))
  
  # Process data
  message("Processing data...")
  seurat_obj <- process_seurat_object(seurat_obj)
  
  # Annotate cell types
  message("Annotating cell types...")
  seurat_obj <- annotate_cell_types(seurat_obj)
  
  # Assign disease groups
  message("Assigning disease groups...")
  seurat_obj@meta.data <- assign_disease_groups(seurat_obj@meta.data)
  
  # Create cell type subsets
  message("Creating cell type subsets...")
  cell_subsets <- define_cell_subsets(seurat_obj)
  
  # Save results
  message("Saving results...")
  saveRDS(seurat_obj, file.path(OUTPUT_DIR, "processed_seurat_object.rds"))
  
  # Save cell type subsets
  walk2(cell_subsets, names(cell_subsets), function(subset, name) {
    saveRDS(subset, file.path(OUTPUT_DIR, paste0(tolower(name), "_subset.rds")))
  })
  
  # Analyze cell counts
  message("Analyzing cell counts...")
  cell_counts <- analyze_cell_counts(cell_subsets)
  write_tsv(cell_counts, file.path(OUTPUT_DIR, "cell_type_counts.tsv"))
  
  # Create UMAP visualization
  message("Creating UMAP plot...")
  create_umap_plot(seurat_obj, file.path(OUTPUT_DIR, "umap_plot.png"))
  
  message("Analysis pipeline completed successfully!")
}

# Run pipeline
if (!interactive()) {
  main()
}
