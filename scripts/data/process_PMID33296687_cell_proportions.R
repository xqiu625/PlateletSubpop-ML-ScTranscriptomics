#!/usr/bin/env Rscript

#' Cell Proportion Analysis Pipeline for COVID-19 Dataset
#' Dataset: Bernardes et al., 2020 (PMID: 33296687)
#' Description: Analysis of cell type proportions in COVID-19 patients
#' Data: FastGenomics Portal (https://beta.fastgenomics.org/datasets/detail-dataset-0b28ff50dc64451896628dfdd5198fd1)

# 1. Setup and Configuration -------------------------------------------
suppressPackageMessages({
  library(tidyverse)
})

# Set paths
BASE_DIR <- "/datasets"
INPUT_FILE <- file.path(BASE_DIR, "PMID33296687/")
OUTPUT_DIR <- file.path(BASE_DIR, "PMID33296687")
R_LIB_PATH <- "/R"

.libPaths(c(R_LIB_PATH, .libPaths()))

# Define metadata
STUDY_METADATA <- list(
  Disease = "COVID-19",
  Data_ID = "Bernardes et al."
)

# Define cell type groupings
CELL_TYPE_GROUPS <- list(
  Monocyte = "Monocytes",
  Platelet = "Megakaryocytes",
  T_cell = c("CD4..T.cells", "CD8..T.cells", "Proliferative.Lymphocytes"),
  B_cell = c("Bcells", "Plasmablasts"),
  Other_cell = c("Cell.precursors..HSCs..MEPs..CMPs..GMPs.",
                 "Dendritic.cells", "Granulocytes", "NK.cells")
)

# 2. Data Processing Functions ----------------------------------------
#' Process Cell Proportions
#' @param data Input data frame with cell proportions
#' @return Processed data frame with cell counts and percentages
process_cell_proportions <- function(data) {
  # Add metadata
  processed_data <- data %>%
    rename(Sample_ID = Sample) %>%
    mutate(
      Disease = STUDY_METADATA$Disease,
      Data_ID = STUDY_METADATA$Data_ID,
      Disease_group = paste(Pseudotime_name, Time.point, sep = "_")
    )
  
  # Calculate cell type groups
  cell_counts <- processed_data %>%
    transmute(
      Data_ID = Data_ID,
      Sample_ID = Sample_ID,
      Disease = Disease,
      Disease_group = Disease_group,
      Monocyte = Monocytes,
      Platelet = Megakaryocytes,
      T_cell = reduce(select(processed_data, !!!CELL_TYPE_GROUPS$T_cell), `+`),
      B_cell = reduce(select(processed_data, !!!CELL_TYPE_GROUPS$B_cell), `+`),
      Other_cell = reduce(select(processed_data, !!!CELL_TYPE_GROUPS$Other_cell), `+`)
    )
  
  # Handle missing values
  cell_counts[is.na(cell_counts)] <- 0
  
  # Calculate percentages
  count_cols <- cell_counts %>%
    select(-c(Data_ID:Disease_group))
  
  total_cells <- rowSums(count_cols)
  
  percentages <- map_df(count_cols, ~(. / total_cells) * 100) %>%
    rename_with(~paste0("PCT_", .), everything())
  
  # Combine results and add empty Erythroid columns
  final_data <- bind_cols(cell_counts, percentages) %>%
    mutate(
      Erythroid = NA_real_,
      PCT_Erythroid = NA_real_
    ) %>%
    select(
      Data_ID, Sample_ID, Disease, Disease_group,
      Monocyte, Platelet, Erythroid, T_cell, B_cell, Other_cell,
      starts_with("PCT_")
    )
  
  return(final_data)
}

# 3. Main Execution ------------------------------------------------
main <- function() {
  # Create output directory
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Read input data
  message("Reading cell proportion data...")
  raw_data <- read_csv2(INPUT_FILE)
  
  # Process data
  message("Processing cell proportions...")
  processed_data <- process_cell_proportions(raw_data)
  
  # Save results
  message("Saving results...")
  output_file <- file.path(OUTPUT_DIR, "cell_type_proportions.tsv")
  write_tsv(processed_data, output_file)
  
  message(sprintf("Results saved to: %s", output_file))
  message("Analysis completed successfully!")
}

# Run pipeline
if (!interactive()) {
  main()
}
