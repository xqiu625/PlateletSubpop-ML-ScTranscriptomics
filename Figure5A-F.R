library(dplyr)
library(Seurat)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(limma)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(rstatix)

# Load gene sets
genesets <- getGmt("h.all.v7.5.1.symbols.gmt")

# Load data and prepare
dat0 <- readRDS("Platelets_cluster_Anno.rds")
DefaultAssay(dat0) <- "RNA"

# Extract expression data
df.data <- GetAssayData(object = dat0, slot = "data")

# Perform GSVA
gsvascore <- gsva(data.matrix(df.data), genesets, parallel.sz = 2)
saveRDS(gsvascore, 'gsvascore.RDS')


# Prepare data frame for heatmap annotation
df.group <- data.frame(umi = rownames(dat0@meta.data),
                       cluster = dat0$seurat_clusters,
                       stringsAsFactors = FALSE)

# Create a heatmap annotation based on cluster
ha.t <- HeatmapAnnotation(Cluster = df.group$cluster)

# Generate heatmap
p <- Heatmap(as.matrix(gsvascore),
             show_column_names = FALSE,
             cluster_rows = TRUE,
             cluster_columns = TRUE,
             top_annotation = ha.t,
             column_split = df.group$cluster,
             row_names_gp = gpar(fontsize = 8),
             row_names_max_width = max_text_width(rownames(gsvascore), gp = gpar(fontsize = 8)))

# Save the heatmap to a PNG file
dpi = 300
png(file = 'seurat_cluster_gsvascore.png',
    width = dpi * 30, height = dpi * 6, units = "px", res = dpi, type = 'cairo')
print(p)
dev.off()


# Define function for plotting
plot_hallmark_scores <- function(hallmark_name, df_gsvascore_combined) {
  df_filtered <- df_gsvascore_combined %>% filter(hallmark == hallmark_name)
  
  stat.test <- wilcox_test(score ~ cluster, data = df_filtered) %>%
    add_xy_position(x = "cluster")
  
  p <- ggplot(df_filtered, aes(x = cluster, y = score, fill = cluster)) +
    geom_bar(stat = "summary", fun = "mean") +
    geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2) +
    scale_fill_manual(values = rainbow(n_distinct(df_filtered$cluster))) +
    theme_minimal() +
    labs(title = hallmark_name, x = "Cluster", y = "GSVA Score") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), legend.position = "none")

  
  ggsave(paste0(hallmark_name, "_cluster_hallmark.png"), plot = p, width = 12, height = 4, dpi = 300)
}

# Loop through each hallmark for plotting
hallmarks <- unique(df_gsvascore_combined$hallmark)
for (hallmark_name in hallmarks) {
  plot_hallmark_scores(hallmark_name, df_gsvascore_combined)
}