library(ggplot2)

# Sample data
data <- data.frame(
  model = rep(c("DNN", "XGB"), each = 4),
  metric = rep(c("accuracy", "f1", "precision", "recall"), 2),
  value = c(0.8945, 0.9326, 0.9358, 0.9294, 0.9049, 0.7519, 0.6556, 0.8815)
)

# File path to save the plot
file_path <- "mlmodel_metric_compare.png"

# Open a PNG device
png(filename = file_path, width = 8, height = 6, units = 'in', res = 300)

# Create the plot
ggplot(data = data, aes(x = metric, y = value, fill = model)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  theme_minimal() +
  scale_fill_manual(values = c('#91AD9E','#8A95A9')) +
  geom_text(aes(label = value), vjust = 1.6, color = "white",
            position = position_dodge(0.9), size = 3.5)

# Close the PNG device
dev.off()
