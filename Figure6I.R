library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(pROC)
library(cowplot)

# Load data
df <- read.delim2('per_cell_type_annotation.tsv')
dat0 <- readRDS("PlateletsOnly_cluster_Anno.rds")

# Preprocessing for ROC analysis
df <- distinct(df, Sample_ID, Outcome, Severity)
dat01 <- dat0@meta.data %>%
  mutate(seurat_clusters = paste0("C", seurat_clusters)) %>%
  group_by(Sample_ID, seurat_clusters) %>%
  count(name = "n") %>%
  spread(key = seurat_clusters, value = n, fill = 0) %>%
  left_join(dat03, by = "Sample_ID")

# Calculate percentage of each cluster
total_col <- rowSums(dat01[,-(1:3)]) # Exclude first three columns: Sample_ID, Outcome, Severity
pcts <- as.data.frame(sweep(dat01[,-(1:3)], 1, total_col, FUN = "/") * 100)
pcts2 <- bind_cols(dat01[,1:3], pcts) # Rebind with initial columns

# Filter for non-unknown and non-HC Outcome outcomes
meta3 <- pcts2 %>%
  select(Outcome, starts_with("C")) %>%
  filter(Outcome != "Unknown", Outcome != "HC") %>%
  mutate(Outcome = ifelse(Outcome == "NS", 1, 0)) # Assuming NS=Non-survivor as positive class

# ROC Curve Analysis
setwd("~/Platelet_project")
pdf("ROC_platelet-cluster_ratio.pdf", height = 8, width = 8)
mycol <- c("#007b00", "#24e0b8", "#ffcc51", "#ff8b76", "#ff3031", "#ff6600", "#f5de2e", "#d00e80", "#f673f5", "#8fce00", "#ff6600", "#f5de2e", "#d00e80")

# Initialize data frame to store AUC values
auc_out <- data.frame(Name = character(), AUC = numeric(), `AUC CI` = character(), Color = character())

# Loop through clusters to plot ROC and calculate AUC
for (i in 2:ncol(meta3)) {
  roc_res <- roc(meta3$Outcome, meta3[[i]])
  auc_res <- auc(roc_res)
  ci_res <- ci.auc(roc_res)
  
  # Plot ROC Curve
  plot(roc_res, col = mycol[i-1], lwd = 2, legacy.axes = TRUE, main = "", print.thres = FALSE)
  if(i == 2) { # Add legend from the second iteration
    legend("bottomright", legend = names(meta3)[i], col = mycol[i-1], lwd = 2, bty = "n")
  }
  
  # Store AUC results
  auc_out <- rbind(auc_out, c(Name = names(meta3)[i], AUC = as.numeric(auc_res), `AUC CI` = paste(round(ci_res[1], 3), round(ci_res[3], 3), sep = "-"), Color = mycol[i-1]))
}

dev.off()

# Assuming the legend needs to be outside the loop to include all clusters' AUCs
