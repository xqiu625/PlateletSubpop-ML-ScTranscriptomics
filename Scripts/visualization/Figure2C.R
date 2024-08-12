library(dplyr)
library(EnhancedVolcano)

# Set parameters
dpi <- 300

# Create Volcano Plot
png(file = "NS_S_volcano_ML_features.png", width = dpi * 12, height = dpi * 6, units = "px", res = dpi, type = 'cairo')

EnhancedVolcano(markers,
                lab = markers$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = markers_gene,
                xlim = c(min(markers$avg_log2FC) - 0.5, max(markers$avg_log2FC) + 0.5),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'adjusted P value'),
                title = "",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 4.0,
                labSize = 2.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                legendPosition = 'none',
                maxoverlapsConnectors = Inf)

dev.off()
