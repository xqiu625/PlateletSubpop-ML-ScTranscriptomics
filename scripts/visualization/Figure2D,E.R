library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
library(tibble)
library(viridis)
library(forcats)
library(stringr)

dpi <- 300

perform_enrichment <- function(cluster_id) {
  degs_ <- markers %>%
    filter(gene %in% markers_gene, cluster == cluster_id)
  
  gene_ <- degs_$gene
  B_fun.gene.1 <- bitr(gene_, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  B_fun.gene <- B_fun.gene.1$ENTREZID
  
  B_fun_bp <- enrichGO(
    gene = B_fun.gene,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  ) %>%
    setReadable('org.Hs.eg.db', 'ENTREZID')
  
  # File naming based on cluster ID
  file_name <- paste0("ML_feature_", cluster_id, "_go.png")
  
  png(file = file_name, width = dpi * 12, height = dpi * 6, units = "px", res = dpi, type = 'cairo')
  
  print(
    barplot(B_fun_bp, showCategory = 20, drop = TRUE) +
      geom_col(width = 0.8) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
      theme(
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)
      )
  )
  
  dev.off()
}

# Perform enrichment for "S" and "NS" clusters
perform_enrichment("S")
perform_enrichment("NS")
