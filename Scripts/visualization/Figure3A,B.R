# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(ggpubr)
library(stringr)


# Load the dataset
dat0 <- readRDS("Platelets_cluster_Anno.rds")
DefaultAssay(dat0) <- "RNA"


# Function to perform operations and generate plots for a given group.by parameter
generate_plots <- function(df, group_by_var, feature, file_name) {
  df_avg <- AverageExpression(
    df,
    group.by = group_by_var,
    slot = "data",
    return.seurat = TRUE
  )
  
  df_avg_frame <- as.data.frame(df_avg@assays[["RNA"]]@data)
  cluster_averages_df <- subset(df_avg_frame, row.names(df_avg_frame) == feature)
  print(t(cluster_averages_df))
  
  dpi <- 300
  
  # Save violin plot
  png(file = paste0(file_name, ".png"), width = dpi * 4, height = dpi * 4, units = "px", res = dpi, type = 'cairo')
  VlnPlot(
    df,
    features = feature,
    sort = 'increasing',
    group.by = group_by_var,
    adjust = 1,
    y.max = NULL,
    same.y.lims = FALSE,
    log = FALSE,
    ncol = NULL,
    combine = TRUE,
    pt.size = 0
  )
  dev.off()
}

generate_plots(dat0, "Severity", "ITGA2B", "ITGA2B_Severity")
generate_plots(dat0, "Severity", "actin cytoskeleton reorganization", "actin_cytoskeleton_reorganization_Severity")
