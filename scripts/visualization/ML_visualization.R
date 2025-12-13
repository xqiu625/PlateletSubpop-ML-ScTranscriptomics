#!/usr/bin/env Rscript

#' Machine Learning Visualization Functions
#' Description: Generate publication-quality figures for ML model evaluation
#'
#' Includes:
#' - ROC curves (single and comparison)
#' - Confusion matrices
#' - Feature importance plots
#' - Model performance comparison
#' - Cross-validation results

# 1. Setup -------------------------------------------------------------
suppressPackageMessages({
  library(ggplot2)
  library(pROC)
  library(ROCR)
  library(cowplot)
  library(gridExtra)
  library(RColorBrewer)
  library(viridis)
  library(tidyverse)
  library(scales)
})

# Color palettes
MODEL_COLORS <- c(
  "XGBoost" = "#E41A1C",
  "RandomForest" = "#4DAF4A",
  "SVM_Radial" = "#377EB8",
  "SVM_Linear" = "#984EA3",
  "ElasticNet" = "#FF7F00",
  "GBM" = "#A65628",
  "DNN" = "#F781BF"
)

OUTCOME_COLORS <- c(
  "S" = "#4DAF4A",     # Survivor - green
  "NS" = "#E41A1C"     # Non-survivor - red
)

# 2. ROC Curve Functions -----------------------------------------------
#' Plot Single ROC Curve
#' @param roc_obj ROC object from pROC
#' @param model_name Name of the model
#' @param color Color for the curve
#' @return ggplot object
plot_single_roc <- function(roc_obj, model_name = "Model", color = "#E41A1C") {
  auc_val <- as.numeric(auc(roc_obj))

  # Extract coordinates
  roc_df <- data.frame(
    Sensitivity = roc_obj$sensitivities,
    Specificity = roc_obj$specificities
  )

  ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity)) +
    geom_path(color = color, size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    annotate("text", x = 0.7, y = 0.2, size = 5,
             label = sprintf("AUC = %.3f", auc_val)) +
    labs(
      title = sprintf("%s ROC Curve", model_name),
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)"
    ) +
    coord_equal() +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

#' Plot Multiple ROC Curves for Comparison
#' @param results List of model evaluation results
#' @param title Plot title
#' @return ggplot object
plot_roc_comparison <- function(results, title = "ROC Curve Comparison") {

  # Combine all ROC data
  roc_data <- lapply(names(results), function(name) {
    roc_obj <- results[[name]]$roc
    data.frame(
      Model = name,
      Sensitivity = roc_obj$sensitivities,
      Specificity = roc_obj$specificities,
      AUC = results[[name]]$auc
    )
  }) %>% bind_rows()

  # Create legend labels with AUC
  auc_labels <- sapply(names(results), function(name) {
    sprintf("%s (AUC = %.3f)", name, results[[name]]$auc)
  })

  p <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity,
                            color = Model, group = Model)) +
    geom_path(size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = MODEL_COLORS, labels = auc_labels) +
    labs(
      title = title,
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)",
      color = "Model"
    ) +
    coord_equal() +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = c(0.7, 0.25),
      legend.background = element_rect(fill = "white", color = "gray80"),
      legend.key.size = unit(0.8, "lines"),
      panel.grid.minor = element_blank()
    )

  return(p)
}

#' Plot ROC with Confidence Interval (Bootstrap)
#' @param roc_obj ROC object from pROC
#' @param model_name Model name
#' @param n_boot Number of bootstrap iterations
#' @return ggplot object
plot_roc_with_ci <- function(roc_obj, model_name = "Model", n_boot = 2000) {
  # Calculate CI
  ci_obj <- ci.se(roc_obj, boot.n = n_boot, progress = "none")

  # Create data frame
  ci_df <- data.frame(
    Specificity = as.numeric(rownames(ci.se(roc_obj))),
    Lower = ci_obj[, 1],
    Sensitivity = ci_obj[, 2],
    Upper = ci_obj[, 3]
  )

  auc_ci <- ci.auc(roc_obj, boot.n = n_boot)

  ggplot(ci_df, aes(x = 1 - Specificity)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "#E41A1C") +
    geom_line(aes(y = Sensitivity), color = "#E41A1C", size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    annotate("text", x = 0.65, y = 0.2, size = 4,
             label = sprintf("AUC = %.3f\n(95%% CI: %.3f - %.3f)",
                           auc_ci[2], auc_ci[1], auc_ci[3])) +
    labs(
      title = sprintf("%s ROC Curve with 95%% CI", model_name),
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)"
    ) +
    coord_equal() +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# 3. Confusion Matrix Functions ----------------------------------------
#' Plot Confusion Matrix Heatmap
#' @param cm Confusion matrix from caret
#' @param model_name Model name
#' @param normalize Whether to show percentages
#' @return ggplot object
plot_confusion_matrix <- function(cm, model_name = "Model", normalize = FALSE) {

  cm_table <- as.data.frame(cm$table)

  if (normalize) {
    # Calculate percentages by reference
    cm_table <- cm_table %>%
      group_by(Reference) %>%
      mutate(Percentage = Freq / sum(Freq) * 100) %>%
      ungroup()

    label_var <- "Percentage"
    label_format <- function(x) sprintf("%.1f%%", x)
    fill_var <- "Percentage"
    legend_title <- "Percentage"
  } else {
    cm_table$Label <- cm_table$Freq
    label_var <- "Freq"
    label_format <- function(x) as.character(x)
    fill_var <- "Freq"
    legend_title <- "Count"
  }

  # Extract metrics
  accuracy <- cm$overall["Accuracy"]
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]

  p <- ggplot(cm_table, aes(x = Reference, y = Prediction, fill = .data[[fill_var]])) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = label_format(.data[[fill_var]])),
              color = "white", size = 6, fontface = "bold") +
    scale_fill_gradient(low = "#3288BD", high = "#D53E4F",
                        name = legend_title) +
    scale_x_discrete(position = "top") +
    labs(
      title = sprintf("%s Confusion Matrix", model_name),
      subtitle = sprintf("Accuracy: %.1f%% | Sensitivity: %.1f%% | Specificity: %.1f%%",
                        accuracy * 100, sensitivity * 100, specificity * 100),
      x = "Actual",
      y = "Predicted"
    ) +
    coord_equal() +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text = element_text(size = 12),
      legend.position = "right",
      panel.grid = element_blank()
    )

  return(p)
}

#' Plot Multiple Confusion Matrices
#' @param results List of model evaluation results
#' @param ncol Number of columns
#' @return Combined ggplot
plot_confusion_matrices_grid <- function(results, ncol = 3) {
  plots <- lapply(names(results), function(name) {
    plot_confusion_matrix(results[[name]]$confusion_matrix, name)
  })

  plot_grid(plotlist = plots, ncol = ncol)
}

# 4. Feature Importance Functions --------------------------------------
#' Plot Feature Importance Bar Chart
#' @param importance_df Data frame with feature importance
#' @param model_name Model name to plot
#' @param top_n Number of top features
#' @return ggplot object
plot_feature_importance <- function(importance_df, model_name, top_n = 20) {

  plot_data <- importance_df %>%
    filter(Model == model_name) %>%
    arrange(desc(Importance)) %>%
    head(top_n)

  color <- MODEL_COLORS[model_name]
  if (is.na(color)) color <- "#377EB8"

  ggplot(plot_data, aes(x = reorder(Gene, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = color, alpha = 0.8) +
    coord_flip() +
    labs(
      title = sprintf("%s Feature Importance (Top %d)", model_name, top_n),
      x = "Gene",
      y = "Importance Score"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 9)
    )
}

#' Plot Consensus Feature Importance
#' @param consensus_df Data frame with consensus features
#' @param top_n Number of features to show
#' @return ggplot object
plot_consensus_features <- function(consensus_df, top_n = 20) {

  plot_data <- head(consensus_df, top_n)

  ggplot(plot_data, aes(x = reorder(Gene, n_models + avg_importance/100),
                        y = avg_importance, fill = factor(n_models))) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_viridis_d(name = "# Models\nSelected", option = "plasma") +
    labs(
      title = "Consensus Top Features Across Models",
      subtitle = "Color indicates number of models that selected the feature",
      x = "Gene",
      y = "Average Importance Score"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.text.y = element_text(size = 9)
    )
}

#' Plot Feature Importance Heatmap
#' @param importance_df Data frame with feature importance
#' @param top_n Number of top features per model
#' @return ggplot object
plot_importance_heatmap <- function(importance_df, top_n = 15) {

  # Get top features from each model
  top_genes <- importance_df %>%
    group_by(Model) %>%
    top_n(top_n, Importance) %>%
    pull(Gene) %>%
    unique()

  # Filter and reshape
  plot_data <- importance_df %>%
    filter(Gene %in% top_genes) %>%
    select(Gene, Model, Importance) %>%
    pivot_wider(names_from = Model, values_from = Importance, values_fill = 0)

  # Convert to matrix for ordering
  mat <- as.matrix(plot_data[, -1])
  rownames(mat) <- plot_data$Gene

  # Order by mean importance
  gene_order <- names(sort(rowMeans(mat), decreasing = TRUE))

  # Reshape for ggplot
  plot_long <- importance_df %>%
    filter(Gene %in% top_genes) %>%
    mutate(Gene = factor(Gene, levels = rev(gene_order)))

  ggplot(plot_long, aes(x = Model, y = Gene, fill = Importance)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "magma") +
    labs(
      title = "Feature Importance Across Models",
      x = "Model",
      y = "Gene",
      fill = "Importance"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8)
    )
}

# 5. Performance Comparison Functions ----------------------------------
#' Plot Model Performance Comparison
#' @param summary_df Summary data frame with metrics
#' @return ggplot object
plot_performance_comparison <- function(summary_df) {

  # Reshape for plotting
  plot_data <- summary_df %>%
    pivot_longer(cols = c(Accuracy, Sensitivity, Specificity, AUC, F1_Score),
                 names_to = "Metric", values_to = "Value")

  ggplot(plot_data, aes(x = Model, y = Value, fill = Model)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~Metric, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = MODEL_COLORS) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = "Model Performance Comparison",
      x = NULL,
      y = "Score"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.text = element_text(face = "bold")
    )
}

#' Plot Radar Chart for Model Comparison
#' @param summary_df Summary data frame
#' @return ggplot object
plot_radar_comparison <- function(summary_df) {
  # Prepare data
  metrics <- c("Accuracy", "Sensitivity", "Specificity", "AUC", "F1_Score")

  plot_data <- summary_df %>%
    select(Model, all_of(metrics)) %>%
    pivot_longer(cols = -Model, names_to = "Metric", values_to = "Value") %>%
    group_by(Model) %>%
    mutate(
      angle = (as.numeric(factor(Metric)) - 1) * 2 * pi / length(metrics),
      x = Value * cos(angle),
      y = Value * sin(angle)
    )

  # Add closing point for each model
  plot_data <- plot_data %>%
    group_by(Model) %>%
    group_modify(~ bind_rows(.x, .x[1, ]))

  ggplot(plot_data, aes(x = x, y = y, color = Model, group = Model)) +
    geom_polygon(aes(fill = Model), alpha = 0.1) +
    geom_path(size = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = MODEL_COLORS) +
    scale_fill_manual(values = MODEL_COLORS) +
    coord_equal() +
    labs(title = "Model Performance Radar Chart") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

# 6. CV Results Functions ----------------------------------------------
#' Plot Cross-Validation Results
#' @param cv_results Resamples object from caret
#' @return ggplot object
plot_cv_results <- function(cv_results) {
  bwplot(cv_results, metric = "ROC",
         main = "Cross-Validation ROC Comparison",
         par.settings = list(
           box.rectangle = list(col = "black"),
           box.umbrella = list(col = "black"),
           plot.symbol = list(col = "black", pch = 19)
         ))
}

#' Plot CV Density
#' @param cv_results Resamples object
#' @return ggplot object
plot_cv_density <- function(cv_results) {
  densityplot(cv_results, metric = "ROC",
              main = "Cross-Validation ROC Distribution",
              auto.key = TRUE)
}

# 7. Master Visualization Function -------------------------------------
#' Generate All Visualizations
#' @param ml_results Results from ML pipeline
#' @param output_dir Directory for saving plots
#' @param format Image format ("png" or "pdf")
#' @param width Plot width in inches
#' @param height Plot height in inches
generate_all_plots <- function(ml_results,
                               output_dir = "results/figures",
                               format = "png",
                               width = 10,
                               height = 8) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("Generating visualizations...")

  # 1. ROC Comparison
  message("  - ROC curves...")
  p_roc <- plot_roc_comparison(ml_results$results)
  ggsave(file.path(output_dir, paste0("roc_comparison.", format)),
         p_roc, width = width, height = height, dpi = 300)

  # 2. Confusion Matrices
  message("  - Confusion matrices...")
  p_cm <- plot_confusion_matrices_grid(ml_results$results)
  ggsave(file.path(output_dir, paste0("confusion_matrices.", format)),
         p_cm, width = width * 1.5, height = height, dpi = 300)

  # 3. Feature Importance (for each model)
  message("  - Feature importance plots...")
  for (model_name in names(ml_results$models)) {
    if (model_name %in% unique(ml_results$feature_importance$Model)) {
      p_imp <- plot_feature_importance(ml_results$feature_importance, model_name)
      ggsave(file.path(output_dir,
                       paste0("feature_importance_", model_name, ".", format)),
             p_imp, width = width, height = height, dpi = 300)
    }
  }

  # 4. Consensus Features
  message("  - Consensus features...")
  p_consensus <- plot_consensus_features(ml_results$consensus_features)
  ggsave(file.path(output_dir, paste0("consensus_features.", format)),
         p_consensus, width = width, height = height, dpi = 300)

  # 5. Importance Heatmap
  message("  - Importance heatmap...")
  p_heatmap <- plot_importance_heatmap(ml_results$feature_importance)
  ggsave(file.path(output_dir, paste0("importance_heatmap.", format)),
         p_heatmap, width = width, height = height * 1.2, dpi = 300)

  # 6. Performance Comparison
  message("  - Performance comparison...")
  p_perf <- plot_performance_comparison(ml_results$summary)
  ggsave(file.path(output_dir, paste0("performance_comparison.", format)),
         p_perf, width = width * 1.2, height = height * 0.6, dpi = 300)

  message(sprintf("All plots saved to: %s", output_dir))
}

# 8. Example Usage -----------------------------------------------------
if (!interactive()) {
  message("ML Visualization Functions Loaded")
  message("\nUsage:")
  message("  # Load results")
  message("  ml_results <- readRDS('results/ml_comparison/ml_pipeline_results.rds')")
  message("")
  message("  # Generate all plots")
  message("  generate_all_plots(ml_results)")
  message("")
  message("  # Or generate individual plots:")
  message("  plot_roc_comparison(ml_results$results)")
  message("  plot_confusion_matrix(ml_results$results$XGBoost$confusion_matrix, 'XGBoost')")
  message("  plot_feature_importance(ml_results$feature_importance, 'RandomForest')")
}
