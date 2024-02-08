# Load necessary libraries
library(dplyr)
library(SeuratDisk)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(RColorBrewer)

# Set the working directory to where DEG files are located
setwd("~/DEG/Platelets")

# Function to process DEG files
process_DEG_file <- function(file_name) {
  # Construct the file path and read the DEG data
  degs <- read.table(file_name, sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
  
  # Loop through each cluster in the DEG data
  unique(degs$cluster) %>%
    walk(function(cluster) {
      degs_clustered <- degs %>%
        filter(cluster == cluster & p_val_adj < 0.05 & pct.1 != 0 & pct.2 != 0 & avg_log2FC > 0.25)
      
      if(nrow(degs_clustered) > 1) {
        perform_enrichment_analysis(degs_clustered, gsub(".txt","",basename(file_name)), cluster)
      }
    })
}

# Function to perform enrichment analysis and save results
perform_enrichment_analysis <- function(degs, name, cluster) {
  gene_ids <- bitr(degs$gene, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)$ENTREZID
  
  # GO enrichment analysis
  B_fun_bp <- enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  
  # KEGG pathway analysis
  B_fun_kegg <- enrichKEGG(gene = gene_ids, organism = 'hsa', pvalueCutoff = 0.05)
  
  # Save results and plots for GO and KEGG
  save_results_and_plots(B_fun_bp, paste0(name, "_", cluster, '_go-bp'), "")
  save_results_and_plots(B_fun_kegg, paste0(name, "_", cluster, '_kegg'), "")
}

# Function to save enrichment results and generate plots
save_results_and_plots <- function(enrichment_results, file_prefix, output_dir) {
  if(dim(enrichment_results)[1] > 0) {
    write.table(enrichment_results, file = paste0(output_dir, "/", file_prefix, '.txt'), sep='\t', quote = FALSE, row.names = F)
    
    dpi <- 300
    png(file = paste0(output_dir, "/", file_prefix, '.png'), width = dpi * 12, height = dpi * 6, units = "px", res = dpi, type = 'cairo')
    barplot(enrichment_results, showCategory = 20, drop = TRUE) +
      theme(axis.text.y = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16))
    dev.off()
  }
}

# Process each file
list.files(pattern = "*.txt") %>% map(~process_DEG_file(.))


