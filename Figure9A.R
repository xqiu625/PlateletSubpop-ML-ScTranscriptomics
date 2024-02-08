# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)
library(pheatmap)
library(RColorBrewer)

# Set working directory to Interaction Score files
setwd('Interaction_score/')

# Read interaction scores and annotations
file_names <- list.files(pattern = ".*SUM")
interaction_scores <- do.call("rbind", lapply(file_names, read.delim2))

# Read annotations
annotations <- read.delim2('/bigdata/godziklab/shared/Xinru/302005/datasets/per_cell_type_annotation.tsv') %>%
  mutate(
    Disease2 = ifelse(Disease %in% c('Flu', 'Lung_infect', 'Non_covid'), 'SSH', Disease),
    Severity = ifelse(Severity == "", "SLE", Severity)
  ) %>%
  filter(pbmc == "Y", Severity != "HC_LPS") %>%
  select(Sample_ID, Severity, Outcome, Disease2) %>%
  distinct(Sample_ID, .keep_all = TRUE)

# Merge scores with annotations
interaction_scores <- left_join(interaction_scores, annotations, by = "Sample_ID") %>%
  mutate(
    Score = as.numeric(Score),
    Outcome = ifelse(Severity == "HC", "HC", Outcome)
  )

# Function to generate heatmap for a given grouping variable
generate_heatmap <- function(data, group_var, file_name) {
  data_prep <- data %>%
    group_by(CellType_Pair, !!sym(group_var)) %>%
    summarise(Mean = mean(Score, na.rm = TRUE)) %>%
    mutate(zscore = scale(Mean)) %>%
    select(CellType_Pair, !!sym(group_var), zscore) %>%
    spread(!!sym(group_var), zscore) %>%
    select(-CellType_Pair) %>%
    as.data.frame()
  
  row.names(data_prep) <- data$CellType_Pair
  matrix_data <- as.matrix(data_prep)
  
  pheatmap(matrix_data,
           border_color = NA,
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(200),
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize = 10,
           main = ""
  )
  
  ggsave(file_name, width = 4, height = 3.2, units = "in", dpi = 300)
}

# Generate heatmaps
setwd('/Interaction_score_figure')
generate_heatmap(interaction_scores, "Severity", "receptorLigand_Severity.png")
generate_heatmap(interaction_scores, "Disease2", "receptorLigand_Disease.png")
generate_heatmap(interaction_scores, "Outcome", "receptorLigand_Outcome.png")
