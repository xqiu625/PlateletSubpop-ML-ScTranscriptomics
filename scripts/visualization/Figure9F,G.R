# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggpubr)
library(rstatix)
library(ggplot2)

# Load dataset
df <- readRDS('Module_meta_4.rds')

df$Outcome <- factor(df$Outcome,
                     levels = c("Unknown", "HC", "S", "FT"))
# Define modules to analyze
modules <- c("T.cell.differentiation", "B.cell.proliferation")

for (i in 1: length(modules)){
  module_name = modules[i]
  print(module_name)
  df1 <-df %>%
    dplyr::select(Outcome,modules[i]) %>%
    dplyr::filter(!is.na(Outcome))
  df1[,2]<- as.numeric(df1[,2])
  names(df1) <- c("Outcome", "Score")
  stat.test <- df1 %>% wilcox_test(Score ~ Outcome)
  stat.test <- stat.test %>%
    add_xy_position(fun = "max",x = "Outcome") %>%
    filter(p.adj.signif != "ns")
  dpi = 300
  png(file = paste0(module_name,"_Outcome_module.png"),
      width = dpi * 4,height = dpi * 4,
      units = "px",res = dpi,type = 'cairo')
  print(ggbarplot(df1,
                  x = "Outcome",
                  y =  "Score",
                  fill = "Outcome",
                  palette = c('#FFD16F', '#A1CEED', '#456681', '#ad3c53'),
                  add = c("mean_se")) +
          theme(plot.title = element_text(size=10)) +
          
          # stat_pvalue_manual(stat.test,
          #                    y.position = 0.5,
          #                    step.increase = 0.05,
          #                    tip.length = 0.01,
          #                    label = "p.adj.signif") +
          theme_minimal() +
          ggtitle(module_name) +
          theme(legend.position = "none") +
          theme_minimal(base_size = 11)
  )
  dev.off()
}
