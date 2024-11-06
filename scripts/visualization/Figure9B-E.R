library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)

file_names <- list.files(pattern = ".*T_cell_interaction_score.txt")
data <- do.call("rbind",lapply(file_names,FUN=function(file_names){read.delim2(file_names)}))


plot_interaction_scores <- function(data, group_by_var, file_name) {
  # Prepare data
  data_prepared <- data %>%
    filter(CellType_Pair == 'Platelets - T_cells') %>%
    filter(LR_Pair %in% c('CXCL12 - CXCR4', 'CCL2 - CCR2', 'CCL3 - CCR5')) %>%
    group_by(LR_Pair, !!sym(group_by_var)) %>%
    mutate(Mean = mean(Score, na.rm = TRUE)) %>%
    distinct(LR_Pair, !!sym(group_by_var), .keep_all = TRUE) %>%
    mutate(zscore = (Mean - mean(Mean)) / sd(Mean)) %>%
    select(LR_Pair, !!sym(group_by_var), zscore) %>%
    mutate(!!sym(group_by_var) := factor(!!sym(group_by_var), levels = c("HC", "CV", "ML", "MD", "SV", "FT", "SLE", "Unknown", "S", "NS", "SSH", "COVID-19", "Sepsis")))
  
  # Plot
  p <- ggplot(data_prepared, aes_string(x = group_by_var, y = "LR_Pair")) +
    geom_point(aes(size = zscore), color = '#bcbddc') +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(size = 14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title = element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  
  # Save plot
  ggsave(file_name, plot = p, width = 4, height = 2, units = "in", dpi = 300)
}

file_names <- list.files(pattern = ".*B_cell_interaction_score.txt")
data <- do.call("rbind",lapply(file_names,FUN=function(file_names){read.delim2(file_names)}))

plot_interaction_scores <- function(data, grouping_var, file_name) {
  data %>%
    filter(CellType_Pair == 'Platelets - B_cells') %>%
    filter(grepl("CD40$",LR_Pair)) %>%
    group_by(LR_Pair, {{ grouping_var }}) %>%
    mutate(Mean = mean(Score, na.rm = TRUE)) %>%
    distinct(LR_Pair, {{ grouping_var }}, .keep_all = TRUE) %>%
    mutate(zscore = (Mean - mean(Mean)) / sd(Mean)) %>%
    ggplot(aes(x = {{ grouping_var }}, y = LR_Pair, size = zscore)) +
    geom_point(color = '#bcbddc') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title = element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid")) +
    ggsave(file_name, width = 4, height = 2, dpi = 300)
}

