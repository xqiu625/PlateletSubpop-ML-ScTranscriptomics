
# Load required libraries
library(SeuratDisk)
library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(GO.db)
library(org.Hs.eg.db)
library(KEGGREST)
library(HGNChelper)
library(tidyr)
library(ggpubr)
library(rstatix)

# Load dataset
df <- readRDS("Platelets.rds")
DefaultAssay(dat0) <- "RNA"

# Function to calculate module scores
add_module_score <- function(df, pathway_id, module_name) {
  gene_data <- as.data.frame(keggGet(pathway_id)[[1]]$GENE)
  names(gene_data) <- c("gene")
  gene_data <- gene_data %>% dplyr::filter(grepl(";", gene))
  gene_data$gene <- str_split_fixed(gene_data$gene, ";", 2)[,1]
  gene_data_checked <- checkGeneSymbols(gene_data$gene, species="human") %>%
    dplyr::filter(!is.na(Suggested.Symbol)) %>%
    dplyr::distinct(Suggested.Symbol, .keep_all = TRUE)
  gene_list <- list(as.list(gene_data_checked$Suggested.Symbol))
  df <- AddModuleScore(df, features = gene_list, name = module_name, search = TRUE)
  return(df)
}

# Add module scores
df <- add_module_score(df, "hsa05022", "neurodegeneration")
df <- add_module_score(df, "hsa05171", "COVID-19")

# Save metadata
saveRDS(df@meta.data, "Module_meta.rds")

# Read processed data
df_viz <- readRDS("Module_meta.rds")
modules <- c("COVID-19", "neurodegeneration")

# Function to generate bar plots
plot_module <- function(df, variable, module_name, palette, output_dir) {
  df1 <- df %>% dplyr::select(!!variable, !!module_name) %>% dplyr::filter(!is.na(!!sym(variable)))
  names(df1) <- c("Group", "Module")
  df1$Module <- as.numeric(df1$Module)
  stat.test <- df1 %>% wilcox_test(Module ~ Group)
  stat.test <- stat.test %>% add_xy_position(fun = "max", x = "Group") %>% filter(p.adj.signif != "ns")
  
  png(file = file.path(output_dir, paste0(module_name, "_", variable, "_module.png")), width = 1200, height = 1200, res = 300)
  print(
    ggbarplot(df1, x = "Group", y = "Module", fill = "Group", palette = palette, add = c("mean_se")) +
      theme_minimal() +
      ggtitle(module_name) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  )
  dev.off()
}

# Plot settings
output_dir <- "/Module_score/"
palette_category <- c('#fee5d9','#fcbba1','#fc9272','#fb6a4a','#de2d26','#a50f15', '#bdbdbd')
palette_disease <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c')
palette_death <- c('#FFBE7A','#8ECFC9','#82B0D2','#FA7F6F')

# Generate plots
for (module in modules) {
  plot_module(df_viz, "Category", module, palette_category, output_dir)
  plot_module(df_viz, "Disease2", module, palette_disease, output_dir)
  plot_module(df_viz, "Death", module, palette_death, output_dir)
}
