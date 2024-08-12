# Load necessary libraries
library(Seurat)
library(ggplot2)
library(cowplot) # For theme_cowplot
library(reshape2) # For melt function


# Load the Seurat object
identity <- readRDS('PlateletsOnly_cluster_Anno.rds')

# Preparing the data for plotting
# Transposing and selecting features for simplicity
features <- c('MYL9', 'TMEM140', 'DEPP1', 'JUNB', 'HPSE', 'JCHAIN', 'CRBN', 'S100A8', 'RPLP0', 'HBB', 'CD74', 'TPT1', 'MALAT1')
pbmc <- as.data.frame(t(identity@assays[["RNA"]]@data[features, ]))
pbmc$Cell <- rownames(pbmc)

# Update Seurat cluster identities in metadata
identity@meta.data$seurat_clusters <- factor(paste0("C", identity@meta.data$seurat_clusters),
                                             levels = paste0("C", 0:12))

# Adding cluster identities to the data frame
pbmc$Idents <- identity@meta.data$seurat_clusters[rownames(pbmc)]

# Melting the data frame for ggplot
pbmc_melt <- reshape2::melt(pbmc, id.vars = c("Cell", "Idents"), variable.name = "Feat", value.name = "Expr")

# Generating the plot
a <- ggplot(pbmc_melt, aes(x = Idents, y = Expr, fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 20) +
  theme(legend.position = "none",
        panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 20, face = "bold"),
        strip.text.y.left = element_text(angle = 0)) +
  labs(title = "", x = "Cluster", y = "Expression Level")

# Saving the plot
ggsave("cluster_stacked_violin_marker.png", plot = a, width = 12, height = 24, dpi = 300, units = "in", type = 'cairo')
