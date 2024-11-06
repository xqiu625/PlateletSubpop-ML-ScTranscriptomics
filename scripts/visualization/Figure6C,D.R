# Load necessary libraries
library(dplyr)
library(Seurat)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Set global options and paths
dpi <- 300

# Function to perform differential expression analysis and generate plots
analyze_clusters <- function(dat0, cluster1, cluster2, output_prefix) {
  # Subset data for the clusters of interest
  subset_clusters <- subset(dat0, seurat_clusters %in% c(cluster1, cluster2))
  Idents(subset_clusters) <- 'seurat_clusters'
  
  # Find markers
  markers <- FindAllMarkers(object = subset_clusters, assay = 'RNA', only.pos = TRUE, test.use = 'MAST')
  saveRDS(markers, paste0(output_prefix, '_Markers.rds'))
  write.table(markers, file = paste0(output_prefix, "_vs_", cluster1, ".txt"), sep = "\t", quote = F)
  
  # Prepare data for EnhancedVolcano
  markers <- markers %>%
    mutate(avg_log2FC = as.numeric(avg_log2FC),
           p_val_adj = as.numeric(p_val_adj),
           p_val = as.numeric(p_val),
           avg_log2FC = ifelse(cluster == cluster1, -avg_log2FC, avg_log2FC),
           p_val_adj = ifelse(p_val_adj == 0, .Machine$double.eps, p_val_adj)) %>%
    select(gene, p_val_adj, avg_log2FC) %>%
    filter(abs(avg_log2FC) > 1)
  
  # Volcano plot
  png(file = paste0(output_prefix, "_vs_", cluster1, ".png"), width = dpi * 12, height = dpi * 6, units = "px", res = dpi, type = 'cairo')
  EnhancedVolcano(markers,
                  lab = markers$gene,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  xlim = c(min(markers$avg_log2FC) - 0.5, max(markers$avg_log2FC) + 0.5),
                  xlab = bquote(~Log[2]~ 'fold change'),
                  ylab = bquote(~-Log[10]~ 'adjusted P value'),
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  cutoffLineType = 'twodash',
                  cutoffLineWidth = 0.8,
                  pointSize = 4.0,
                  labSize = 5.0,
                  labCol = 'black',
                  labFace = 'bold',
                  boxedLabels = TRUE,
                  legendPosition = 'none',
                  drawConnectors = TRUE,
                  widthConnectors = 1.0,
                  colConnectors = 'black')
  dev.off()
  
  # Perform GO enrichment analysis for upregulated genes in cluster2 compared to cluster1
  gene_ids <- bitr(markers$gene, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)$ENTREZID
  enrichedGO <- enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  enrichedGO <- setReadable(enrichedGO, 'org.Hs.eg.db', 'ENTREZID')
  
  # Dot plot for GO enrichment
  if (!is.null(enrichedGO) && nrow(enrichedGO) > 0) {
    p <- dotplot(enrichedGO, showCategory = 10)
    png(file = paste0(output_prefix, "_vs_", cluster1, "_GO_enrichment.png"), width = dpi * 6, height = dpi * 6, units = "px", res = dpi, type = 'cairo')
    print(p + theme(axis.text.y = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 16)))
    dev.off()
  }
}

# Assuming necessary libraries have been loaded
library(dplyr)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Function to perform analysis and visualization for a specific comparison
analyze_and_plot <- function(data, cluster_target, cluster_reference, direction, output_dir) {
  # Define the condition for filtering based on direction
  direction_condition <- if (direction == "UP") {
    filter(cluster == cluster_target & p_val_adj < 0.05 & avg_log2FC > 1)
  } else {
    filter(cluster == cluster_reference & p_val_adj < 0.05 & avg_log2FC < -1)
  }
  
  # Perform filtering
  edo <- data %>%
    direction_condition %>%
    filter(pct.1 != 0 & pct.2 != 0)
  
  # Proceed if there are genes to analyze
  if (nrow(edo) > 0) {
    gene_ <- edo$gene
    B_fun.gene.1 <- bitr(gene_, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
    B_fun.gene <- B_fun.gene.1$ENTREZID
    
    # EnrichGO analysis
    B_fun_bp <- enrichGO(gene = B_fun.gene, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
    B_fun_bp <- setReadable(B_fun_bp, 'org.Hs.eg.db', 'ENTREZID')
    
    # Dot plot visualization
    if (nrow(B_fun_bp) > 0) {
      p3 <- dotplot(B_fun_bp, showCategory = 10) +
        theme(axis.text.y = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 16)) +
        ggtitle(paste("Cluster", cluster_target, "vs", cluster_reference, direction))
      
      output_file <- file.path(output_dir, paste("C", cluster_target, "vsC", cluster_reference, "_", direction, ".png", sep = ""))
      dpi <- 300
      png(file = output_file, width = dpi * 6, height = dpi * 6, units = "px", res = dpi, type = 'cairo')
      print(p3)
      dev.off()
    }
  }
}

# Load data
data <- readRDS("Platelets_C4vsC0.rds")
output_dir <- "."

# Perform analysis and visualization for each comparison
analyze_and_plot(data, "4", "0", "UP", output_dir)
analyze_and_plot(data, "0", "4", "DOWN", output_dir)

# Repeat for C1 vs C0 comparison
data <- readRDS("Platelets_C1vsC0.rds")
analyze_and_plot(data, "1", "0", "UP", output_dir)
analyze_and_plot(data, "0", "1", "DOWN", output_dir)

