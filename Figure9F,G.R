# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggpubr)
library(rstatix)
library(ggplot2)



# Load dataset
df_viz <- readRDS('Module_meta_4.rds')

# Define modules to analyze
modules <- c("T.cell.differentiation", "B.cell.proliferation")

# Function to generate plots for each module
generate_module_plots <- function(module_name) {
  message("Processing module: ", module_name)
  
  df_filtered <- df_viz %>%
    select(Outcome, all_of(module_name)) %>%
    filter(!is.na(Outcome)) %>%
    mutate(Score = as.numeric(.data[[module_name]]),
           Outcome = factor(Outcome, levels = c("Unknown", "HC", "S", "NS")))
  
  # Perform statistical test
  stat.test <- df_filtered %>%
    wilcox_test(Score ~ Outcome) %>%
    add_xy_position(fun = "max", x = "Outcome") %>%
    filter(p.adj.signif != "ns")
  
  # Generate plot
  plot_path <- sprintf("/%s_Outcome_module.png", module_name)
  png(plot_path, width = 300 * 4, height = 300 * 4, units = "px", res = 300, type = 'cairo')
  ggbarplot(df_filtered,
            x = "Outcome", y = "Score", fill = "Outcome",
            palette = c('#FFBE7A', '#8ECFC9', '#82B0D2', '#FA7F6F'),
            add = "mean_se") +
    stat_pvalue_manual(stat.test, y.position = 0.03, step.increase = 0.01, tip.length = 0.005, label = "p.adj.signif") +
    theme_minimal(base_size = 11) +
    ggtitle(module_name) +
    theme(legend.position = "none")
  dev.off()
}

# Apply function to each module
lapply(modules, generate_module_plots)
