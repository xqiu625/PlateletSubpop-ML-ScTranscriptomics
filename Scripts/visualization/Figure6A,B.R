library(monocle3)
library(ggplot2)


# Load the CellDataSet
cds <- readRDS("pla_cds_C0_C1_C4_C6.rds")
dpi = 300
# Generating Monocle 3 Trajectory plot colored by Seurat clusters
png(file = "Monocle3_Traject_by_cluster.png", width = dpi * 6, height = dpi * 4, units = "px", res = dpi, type = 'cairo')
plot_cells(cds,
           color_cells_by = "seurat_clusters",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = TRUE,
           label_principal_points = TRUE,
           graph_label_size = 1.5) +
  scale_color_manual(values = c('#377eb8', '#4daf4a', '#e41a1c', '#984ea3'))
dev.off()

# Preparing CellDataSet for pseudotime plot
cds <- order_cells(cds, 
                   reduction_method = "UMAP",
                   root_pr_nodes = "Y_144")

# Generating Pseudotime plot
png(file = "pseudo_by_cluster.png", width = dpi * 6, height = dpi * 4, units = "px", res = dpi, type = 'cairo')
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5) 
dev.off()
