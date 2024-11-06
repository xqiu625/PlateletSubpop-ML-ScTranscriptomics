library(dplyr)
library(VennDiagram)
library(readr)
library(tidyr)
library(ggplot2)
library(Cairo)


# Function to process DEG files and write filtered results
process_DEG_files <- function(file_path) {
  name <- gsub(".txt", "", basename(file_path))
  degs <- read_delim(file_path, delim = "\t", col_types = cols(), .default = col_double())
  unique(degs$cluster) %>%
    walk(~{
      degs %>%
        filter(cluster == .x, p_val_adj < 0.05, pct.1 != 0, pct.2 != 0, avg_log2FC > 0.25) %>%
        write_delim(file = paste0('VENN_data/', name, "|", .x, ".txt"), delim = '\t')
    })
}

# Process each DEG file
list.files(pattern = "*.txt") %>% map(~process_DEG_files(.))

# Venn Diagram Function
create_venn_diagram <- function(file_names, output_file_name, colors_palette = "YlGn") {
  dataset <- map(file_names, read_delim, delim = "\t", col_types = cols(), .default = col_double()) %>%
    map(., `[[`, "gene") %>% set_names(gsub(".txt", "", file_names))
  
  # Calculate union and intersections for Venn diagram
  venn_data <- calculate_venn(dataset)
  
  # Venn diagram creation
  CairoPDF(file = output_file_name, width = 9, height = 6)
  plot_venn(venn_data, colors = brewer.pal(n = min(7, length(dataset)), name = colors_palette))
  dev.off()
}

# Function placeholders for 'calculate_venn' and 'plot_venn'
# Implement these functions based on your specific needs and data structure

# Example usage for "vs HC" comparisons
setwd('VENN_data/')
fnames_vs_HC <- c("SLE_HC|HC.txt","COVID-19_HC|HC.txt","Sepsis_HC|HC.txt","SSH|HC.txt")

create_venn_diagram(fnames_vs_HC, "/Platelet_vs_HC-HC_venn.pdf", "YlGn")

# Example usage for "vs Disease" comparisons
fnames_vs_Disease <- c("SLE_HC|SLE.txt", "COVID-19_HC|COVID-19.txt", "Sepsis_HC|Sepsis.txt", "SSH_HC|SSH.txt")

create_venn_diagram(fnames_vs_Disease, "/Platelet_vs_HC-Disease_venn.pdf", "YlOrBr")
