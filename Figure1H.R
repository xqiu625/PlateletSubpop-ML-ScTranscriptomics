library(cowplot)
library(tidyr)
library(ggplot2)
library(stringr)
library(rstatix)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

dat01 <- read.delim2('per_cell_type_annotation.tsv')
dat <- dat01 %>%
  filter(label != "")
dim(dat)
n_cells_01 <- dat %>%
  dplyr::group_by(Sample_ID, label) %>%
  dplyr::count(label) %>%
  tidyr::spread(label, n)

pcts3 <- read.delim2('<path_to_your_data>/ncells_perCellType.txt')
dat0 <- readRDS("<path_to_your_data>/PlateletsOnly_cluster_Anno.rds")
dat01 <- dat0@meta.data
dat02 <- dat01 %>%
  select(Sample_ID, Data_NO, Disease, Outcome, Severity) %>%
  distinct(Sample_ID, .keep_all = T)
meta <- merge(pcts3, dat02, by = "Sample_ID")
meta$Severity <- ifelse(meta$Disease == "SLE", "SLE", meta$Severity)
meta$Outcome <- ifelse(meta$Severity == "HC", "HC", meta$Outcome)
row.names(meta) <- meta$Sample_ID

# Convert percentages to numeric
meta <- meta %>%
  mutate(across(starts_with("PCT_"), as.numeric)) %>%
  mutate(platelet_Tcell_ratio = PCT_Platelets / PCT_T_cells)

meta3 <- meta %>%
  dplyr::select(Outcome, platelet_Tcell_ratio, starts_with("PCT_")) %>%
  dplyr::filter(Outcome != "Unknown", Outcome != "HC")

names(meta3)[1] <- "Outcome"
names(meta3)[2:ncol(meta3)] <- str_replace(names(meta3)[2:ncol(meta3)], "PCT_", "")

library("pROC")
mycol <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')


pdf("ROC_platelet-tcell_ratio.pdf", height = 8, width = 8)
auc.out <- c()
df <- meta3

# Initial ROC plot
x <- plot.roc(df[,1], df[,2], ylim = c(0,1), xlim = c(1,0),
              ci = TRUE,
              main = "",
              print.thres = FALSE,
              col = mycol[2],
              lwd = 2,
              legacy.axes = T) # Standard axis: "1-specificity"

ci.lower <- round(as.numeric(x$ci[1]), 3)
ci.upper <- round(as.numeric(x$ci[3]), 3)
auc.ci <- c(colnames(meta3)[2], round(as.numeric(x$auc), 3), paste(ci.lower, ci.upper, sep = "-"))
auc.out <- rbind(auc.out, auc.ci)

# Loop to plot additional ROC curves
for (i in 3:ncol(df)) {
  x <- plot.roc(df[,1], df[,i],
                add = T,
                ci = TRUE,
                print.thres = FALSE,
                col = mycol[i],
                lwd = 2,
                legacy.axes = T)
  
  ci.lower <- round(as.numeric(x$ci[1]), 3)
  ci.upper <- round(as.numeric(x$ci[3]), 3)
  
  auc.ci <- c(colnames(df)[i], round(as.numeric(x$auc), 3), paste(ci.lower, ci.upper, sep = "-"))
  auc.out <- rbind(auc.out, auc.ci)
}

# Comparing ROC curves
p.out <- c()
for (i in 2:(ncol(df) - 1)) {
  for (j in (i + 1):ncol(df)) {
    p <- roc.test(df[,1], df[,i], df[,j], method = "bootstrap")
    p.tmp <- c(colnames(df)[i], colnames(df)[j], p$p.value)
    p.out <- rbind(p.out, p.tmp)
  }
}

# Saving p-value comparisons to a data frame
p.out <- as.data.frame(p.out)
colnames(p.out) <- c("ROC1", "ROC2", "p.value")

# Saving AUC and AUC CI to a data frame
auc.out <- as.data.frame(auc.out)
colnames(auc.out) <- c("Name", "AUC", "AUC CI")
auc.out$Color <- mycol
auc.out$AUC <- as.numeric(auc.out$AUC)
auc.out <- auc.out[order(-auc.out$AUC),]

# Drawing the legend
legend.name <- paste(colnames(df)[2:length(df)], "AUC", auc.out$AUC, sep = " ")
legend("bottomright",
       legend = legend.name,
       col = auc.out$Color,
       lwd = 2,
       bty = "n")

# Close the PDF device
dev.off()
