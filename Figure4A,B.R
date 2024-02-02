library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)

# Read the Seurat object
dat0 <- readRDS("PlateletsOnly_cluster_Anno.rds")

# Update metadata within the Seurat object
dat0 <- dat0 %>%
  mutate(seurat_clusters = paste0("C", seurat_clusters),
         seurat_clusters = factor(seurat_clusters,
                                  levels = c("C0", "C1", "C2", "C3", "C4", "C5", "C6",
                                             "C7", "C8", "C9", "C10", "C11", "C12")))

# Define DPI for plots
dpi <- 300

# Generate UMAP plot grouped by 'seurat_clusters'
png(file = "UMAP_DataSource_Platelets.png", width = dpi * 12,height = dpi * 9,units = "px",res = dpi,type = 'cairo')
DimPlot(dat0, group.by = "Data_Source",  label.size = 12) + ggtitle("")
dev.off()

# Generate UMAP plot grouped by 'seurat_clusters'
png(file = "UMAP_cluster_Platelets.png", width = dpi * 12, height = dpi * 9, units = "px", res = dpi, type = 'cairo')
DimPlot(dat0, group.by = "seurat_clusters", label.size = 12) + NoLegend()
dev.off()
