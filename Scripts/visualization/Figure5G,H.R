library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(org.Hs.eg.db)
library(stringr)
library(ggnewscale)
library(DOSE)
library(dplyr)
library(ggupset)
library(AnnotationHub)
library(MeSHDbi)
library(meshes)
library(enrichplot)

# AnnotationHub setup
ah <- AnnotationHub(localHub = TRUE)
hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
file_hsa <- hsa[[1]]
db <- MeSHDbi::MeSHDb(file_hsa)

# Load and prepare data
data <- readRDS('Platelets_cluster_marker.rds')
data$cluster <- paste0("C", data$cluster)

# Analyzing cluster C10
group_ <- "C4"
edo <- data %>%
  filter(cluster == group_, p_val_adj < 0.05, avg_log2FC > 0.4, pct.1 != 0, pct.2 != 0)

# Mapping genes to ENTREZ IDs and filtering
gene_ <- edo$gene
B_fun.gene.1 <- bitr(gene_, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
B_fun.gene <- B_fun.gene.1$ENTREZID
names(B_fun.gene.1)[1] <- "gene"
edo2 <- merge(edo, B_fun.gene.1, by = "gene")
edo3 <- edo2[order(-edo2$avg_log2FC),]
edo4 <- edo3 %>% select(ENTREZID, avg_log2FC)
edo6 <- setNames(as.numeric(edo4$avg_log2FC), edo4$ENTREZID)

# Enrichment analysis
enrichMeSH <- enrichMeSH(names(edo6), MeSHDb = db, category = 'C', database = 'gendoo')
enrichMeSH <- setReadable(enrichMeSH, 'org.Hs.eg.db', 'ENTREZID')

# Concept network plot
p3 <- cnetplot(enrichMeSH, foldChange = edo6, colorEdge = TRUE)

# Save plot
output_file <- paste(group_, '_enrichMeSH.png', sep = '')
dpi <- 300
png(file = output_file, width = dpi * 12, height = dpi * 6, units = "px", res = dpi, type = 'cairo')
print(p3 + ggtitle(group_) +
        theme(axis.text.y = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 16)))
dev.off()





for (i in 1:length(clusters)) {
  group_ <- clusters[i]
  
  # Filter for significant DEGs in the cluster
  edo <- data %>%
    filter(cluster == group_ & p_val_adj < 0.05 & avg_log2FC > 0.4 & pct.1 != 0 & pct.2 != 0)
  
  # Map gene symbols to ENTREZ IDs
  gene_ <- edo$gene
  B_fun.gene.1 <- bitr(gene_, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  
  # Filter out mappings that resulted in NAs
  B_fun.gene.1 <- B_fun.gene.1[!is.na(B_fun.gene.1$ENTREZID),]
  
  # Perform GO enrichment analysis
  if (nrow(B_fun.gene.1) > 0) {
    B_fun_bp <- enrichGO(gene = B_fun.gene.1$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
    
    # Merge DEGs data with ENTREZ IDs
    edo2 <- merge(edo, B_fun.gene.1, by.x = "gene", by.y = "SYMBOL")
    edo3 <- edo2[order(-edo2$avg_log2FC),]
    edo4 <- dplyr::select(edo3, ENTREZID, avg_log2FC)
    edo6 <- setNames(as.numeric(edo4$avg_log2FC), edo4$ENTREZID)
    
    # Calculate pairwise term similarity and plot
    if (length(B_fun_bp) > 0) {
      B_fun_bp <- pairwise_termsim(B_fun_bp)
      p3 <- treeplot(B_fun_bp, hclust_method = "average")
      
      # Save plot
      output_path <- paste('clusterProfiler', group_, '_GO.png', sep = '/')
      dpi <- 300
      png(file = output_path, width = dpi * 16, height = dpi * 6, units = "px", res = dpi, type = 'cairo')
      print(p3 + ggtitle(group_) +
              theme(axis.text.y = element_text(size = 10),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 16)))
      dev.off()
    }
  }
}
