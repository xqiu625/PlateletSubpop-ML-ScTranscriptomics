library(Seurat)
library(ggplot2)

# Assuming dat0 is your Seurat object
dat0 <- readRDS("PlateletsOnly_cluster_Anno.rds")
DefaultAssay(dat0) <- "RNA"

# Assuming marker_Survival and marker_Fatal are predefined
# Ensure Survival_list and Fatal_list are correctly formatted as lists of feature vectors
Survival_list <- c("AIF1","FOS","CD74","JUN","JUNB","HLA-DRA","MNDA","RPL39","RPS21","RPS18","EEF1A1","RPS28","RPL34","S100A8","S100A11","S100A12")
Fatal_list <- c("HBA2","HBB","HPSE","SLC25A37","TMCC2")

# Add module scores for "Survival" and "Fatal" markers
dat0 <- AddModuleScore(dat0, features = Survival_list, name = "Survival", search = TRUE)
dat0 <- AddModuleScore(dat0, features = Fatal_list, name = "Fatal", search = TRUE)

# Function to generate ridge plots for "Survival" and "Fatal" scores
generate_ridge_plot <- function(seurat_obj, score_name, group_by, title, file_name) {
  png(file_name, width = dpi * 6, height = dpi * 4, units = "px", res = dpi, type = 'cairo')
  VlnPlot(seurat_obj, features = paste0(score_name, "1"), 
          sort = TRUE, pt.size = 0, fill.by = group_by) + 
    ggtitle(title) + 
    ylab(group_by) +
    theme(legend.position = "none", 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) +
    coord_flip()
  dev.off()
}

# Set working directory and DPI
dpi <- 300

# Generate ridge plots
generate_ridge_plot(dat0, "Fatal", "seurat_clusters", "Fatal", "Fatal_clusters_Ridge.png")
generate_ridge_plot(dat0, "Survival", "seurat_clusters", "Survival", "Survival_clusters_Ridge.png")
generate_ridge_plot(dat0, "Fatal", "Severity", "Fatal", "Fatal_Severity_Ridge.png")
generate_ridge_plot(dat0, "Survival", "Severity", "Survival", "Survival_Severity_Ridge.png")
