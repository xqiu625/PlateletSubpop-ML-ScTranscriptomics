library(ggplot2)
library(dplyr)

# Function to set working directory and plot DPI
plot_dpi <- 300
plot_width <- 6
plot_height <- 4

# Function to read data, preprocess, and generate plot
generate_plot <- function(file_name, outcome_col, legend_title, colors, output_file) {
  # Read and preprocess the data
  df <- read.csv(file_name)
  
  # Update factor levels based on outcome_col
  if (outcome_col == "Severity") {
    levels_map <- c("HC" = "Healthy Control", "CV" = "Convalescence", "ML" = "Mild",
                    "MD" = "Moderate", "SV" = "Severe", "FT" = "Fatal", "SLE" = "SLE")
  } else { # Assumes outcome_col == "Outcomes"
    levels_map <- c("HC" = "Healthy Control", "NS" = "Non-survivor", "S" = "Survivor",
                    "Unknown" = "Unknown")
  }
  
  df[[outcome_col]] <- factor(df[[outcome_col]], levels = names(levels_map))
  df[[outcome_col]] <- levels_map[df[[outcome_col]]]
  
  # Factorize clusters
  df$Clusters <- factor(df$Clusters, levels = c(paste0("C", 0:12)))
  
  # Generate plot
  png(file = output_file, width = plot_width, height = plot_height, units = 'in', res = plot_dpi, type = 'cairo')
  ggplot(df, aes(fill = .data[[outcome_col]], y = Frequency, x = Clusters)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = colors) +
    theme(axis.text = element_text(size = 16, face = "bold")) +
    xlab("Cluster") + ylab("Percentage") +
    guides(fill = guide_legend(title = legend_title)) +
    theme_minimal()
  dev.off()
}

# Generate plots
generate_plot(
  "Supp_Table2.csv", 
  "Severity", 
  "Severity", 
  c('#8dd3c7','#ffffb3','#bebada','#fdb462','#80b1d3','#fb8072','#b3de69'), 
  "Category ~ cluster_density-plot_v2.png"
)

generate_plot(
  "Supp_Table3.csv", 
  "Outcomes", 
  "Outcome", 
  c('#984ea3', '#377eb8','#4daf4a','#e41a1c'), 
  "Death ~ cluster_density-plot_v2.png"
)
