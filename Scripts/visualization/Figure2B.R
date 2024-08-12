library(dplyr)
library(VennDiagram)
library(grid)


# Read data
dnn <- read.csv("feature_importance.csv")
xgb <- read.csv("xbg_importance_matrix.csv")
logit <- read.csv("logit_features.csv")

# Calculate the 95th percentile for Importance and Gain
dnn_percentile_95th <- quantile(dnn$Importance, probs = 0.95)
xgb_percentile_95th <- quantile(xgb$Gain, probs = 0.95)

# Filter data based on the 95th percentile
dnn2 <- dnn %>% filter(Importance >= dnn_percentile_95th)
xgb2 <- xgb %>% filter(Gain >= xgb_percentile_95th)

# Create lists for Venn diagram
features_dnn <- dnn2$Feature
features_xgb <- xgb2$Feature

list_for_venn <- list(DNN = features_dnn, XGB = features_xgb)

# Generate and plot the Venn diagram
venn_plot <- venn.diagram(
  x = list_for_venn,
  filename = NULL
)

# Plotting the Venn diagram directly
grid.newpage()
grid.draw(venn_plot)

# Saving the Venn diagram to a file
png(filename = "venn_diagram.png", width = 800, height = 800)
grid.draw(venn_plot)
dev.off()
