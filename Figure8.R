library(dplyr)
library(ggpubr)
library(rstatix)
library(ggplot2)


modules <- c("response.to.type.I.interferon", "response.to.interferon.beta",
             "coagulation", "OXPHOS", "Glycolysis", "MHC.ClassII",
             "response.to.interferon.gamma", "translational.initiation")

generate_plot <- function(module_name) {
  df <- df_viz %>%
    select(Severity, !!sym(module_name)) %>%
    filter(!is.na(Severity)) %>%
    rename(Score = !!sym(module_name)) %>%
    mutate(Score = as.numeric(Score), Severity = Severity)
  
  stat.test <- df %>% 
    wilcox_test(Score ~ Severity) %>%
    add_xy_position(fun = "max", x = "Severity") %>%
    filter(p.adj.signif != "ns")
  
  p <- ggbarplot(df,
                 x = "Severity", y = "Score", fill = "Severity",
                 palette = c('#fee5d9','#fcbba','#fc9272','#fb6a4a','#de2d26','#a50f5', '#bdbdbd'),
                 add = c("mean_se")) +
    stat_pvalue_manual(stat.test, y.position = 0.05, step.increase = 0.0, tip.length = 0.005, label = "p.adj.signif") +
    theme_minimal(base_size = ) +
    ggtitle(paste(module_name, "Severity Module")) +
    theme(legend.position = "none")
  
  ggsave(paste0(module_name, "_Severity_module.png"), plot = p, device = "png", path = "~/", dpi = 300, width = 4, height = 4)
}

# Generate plots for selected modules
lapply(modules, generate_plot)
