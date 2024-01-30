library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

dat01 <- read.delim2('per_cell_type_annotation.tsv')
dat <- dat01 %>% 
  filter(label != "")
dim(dat)  
n_cells_01 <- dat %>%
  dplyr::group_by(Sample_ID, label) %>%
  dplyr::count(label) %>%
  tidyr::spread(label, n)

# Filter out samples
pcts3 <- read.delim2('ncells_perCellType.txt')
dat0 <- readRDS("PlateletsOnly_cluster_Anno.rds")
dat01 <- dat0@meta.data
dat02 <- dat01 %>% 
  select(Sample_ID, Data_NO, Disease, Outcome, Severity) %>% 
  distinct(Sample_ID, .keep_all = T)
meta <- merge(pcts3, dat02, by = "Sample_ID")

meta$Severity <- ifelse(meta$Disease == "SLE", "SLE", meta$Severity)
meta$Outcome <- ifelse(meta$Severity == "HC", "HC", meta$Outcome)

# Cell type proportion
meta <- meta %>% select("PCT_B_cell", "PCT_Cell_precursors", "PCT_DC",
                        "PCT_Monocyte", "PCT_Neutrophils", "PCT_NK_cell",
                        "PCT_Platelets", "PCT_T_cells", Disease, Outcome, Severity)

# Renaming columns
colnames(meta) <- sub("PCT_", "", colnames(meta))

cell_types <- c("B cell", "Cell precursors", "DC", "Monocyte",
                "Neutrophil", "NK cell", "Platelets", "T cell")

# Loop for analyzing each cell type
for (i in seq_along(cell_types)) {
  # Custom analysis for each cell type
  meta_category <- meta %>%
    select(cell_types[i], Severity) %>%
    filter(Severity != "")
  names(meta_category) <- c("Percentage", "Severity")
  meta_category$Percentage <- as.numeric(meta_category$Percentage)
  stat.test <- meta_category %>% wilcox_test(Percentage ~ Severity)
  stat.test <- stat.test %>% filter(p.adj.signif != "ns")
  meta_category$Severity <- factor(meta_category$Severity,
                                   levels = c("HC", "CV", "ML", "MD", "SV", "FT", "SLE"))
  setwd("<path_to_save_plots>")  # Replace with a generic path or a specific path as needed
  dpi = 300
  png(file = paste0(cell_types[i], "_Severity.png"), 
      width = dpi * 4, height = dpi * 4, units = "px", res = dpi, type = 'cairo')
  print(ggbarplot(meta_category,
                  x = "Severity",
                  y = "Percentage",
                  fill = "Severity",
                  palette = c('#fee5d9','#fcbba1','#fc9272','#fb6a4a','#de2d26','#a50f15', '#bdbdbd'),
                  add = c("mean_se")) +
          stat_pvalue_manual(stat.test,
                             y.position = 50,
                             step.increase = 0.05,
                             tip.length = 0.01,
                             label = "p.adj.signif") +
          theme_minimal(base_size = 11) +
          theme(legend.position = "none") +
          ggtitle(cell_types[i])
  )
  dev.off()
}
