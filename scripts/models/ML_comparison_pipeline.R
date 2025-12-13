#!/usr/bin/env Rscript

#' Comprehensive Machine Learning Pipeline for Platelet Classification
#' Description: Compare multiple classifiers with cross-validation and feature importance
#'
#' Models included:
#' - XGBoost
#' - Random Forest
#' - Support Vector Machine (SVM)
#' - Logistic Regression
#' - Gradient Boosting Machine (GBM)

# 1. Setup and Configuration -------------------------------------------
suppressPackageMessages({
  library(Seurat)
  library(tidyverse)
  library(caret)
  library(xgboost)
  library(randomForest)
  library(e1071)           # SVM
  library(glmnet)          # Logistic Regression with regularization
  library(gbm)             # Gradient Boosting
  library(pROC)
  library(ROCR)
  library(doParallel)      # Parallel processing
})

# Configuration
CONFIG <- list(
  random_seed = 42,
  train_ratio = 0.8,
  cv_folds = 5,
  cv_repeats = 3,
  top_features = 200,
  n_cores = parallel::detectCores() - 1
)

# 2. Data Preparation Functions ----------------------------------------
#' Prepare Data for Machine Learning
#' @param seurat_obj Seurat object
#' @param outcome_col Column name for outcome variable
#' @param outcome_levels Vector of outcome levels to include
#' @param n_features Number of top variable features to select
#' @return List containing features, labels, and gene names
prepare_ml_data <- function(seurat_obj,
                            outcome_col = "Outcome",
                            outcome_levels = c("NS", "S"),
                            n_features = CONFIG$top_features) {

  message("Preparing data for ML...")

  # Subset to specified outcomes
  cells_to_keep <- seurat_obj[[outcome_col, drop = TRUE]] %in% outcome_levels
  data_subset <- seurat_obj[, cells_to_keep]

  # Extract normalized expression matrix
  expr_matrix <- t(GetAssayData(data_subset, slot = "data"))

  # Get labels
  labels <- factor(data_subset[[outcome_col, drop = TRUE]],
                   levels = outcome_levels)

  # Feature selection: top variable genes
  gene_vars <- apply(expr_matrix, 2, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(n_features, ncol(expr_matrix))])
  expr_matrix <- expr_matrix[, top_genes]

  message(sprintf("Data prepared: %d samples, %d features",
                  nrow(expr_matrix), ncol(expr_matrix)))
  message(sprintf("Class distribution: %s",
                  paste(names(table(labels)), table(labels), sep = "=", collapse = ", ")))

  list(
    features = as.matrix(expr_matrix),
    labels = labels,
    gene_names = top_genes
  )
}

#' Split Data with Stratification
#' @param features Feature matrix
#' @param labels Factor of labels
#' @param train_ratio Proportion for training
#' @return List with train/test splits
split_data <- function(features, labels, train_ratio = CONFIG$train_ratio) {
  set.seed(CONFIG$random_seed)

  train_idx <- createDataPartition(labels, p = train_ratio, list = FALSE)

  list(
    train_x = features[train_idx, ],
    train_y = labels[train_idx],
    test_x = features[-train_idx, ],
    test_y = labels[-train_idx]
  )
}

# 3. Cross-Validation Setup --------------------------------------------
#' Setup Cross-Validation Control
#' @param folds Number of CV folds
#' @param repeats Number of CV repeats
#' @return trainControl object
setup_cv_control <- function(folds = CONFIG$cv_folds,
                             repeats = CONFIG$cv_repeats) {
  trainControl(
    method = "repeatedcv",
    number = folds,
    repeats = repeats,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final",
    verboseIter = FALSE,
    allowParallel = TRUE
  )
}

# 4. Model Training Functions ------------------------------------------
#' Train All Models with Cross-Validation
#' @param train_data Data frame with features and Outcome column
#' @param cv_control Cross-validation control object
#' @return List of trained models
train_all_models <- function(train_data, cv_control) {

  # Setup parallel processing
  cl <- makeCluster(CONFIG$n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))

  models <- list()

  # 1. XGBoost
  message("Training XGBoost...")
  models$XGBoost <- train(
    Outcome ~ ., data = train_data,
    method = "xgbTree",
    trControl = cv_control,
    metric = "ROC",
    tuneGrid = expand.grid(
      nrounds = c(50, 100),
      max_depth = c(3, 6),
      eta = c(0.1, 0.3),
      gamma = 0,
      colsample_bytree = 0.8,
      min_child_weight = 1,
      subsample = 0.8
    ),
    verbose = FALSE
  )

  # 2. Random Forest
  message("Training Random Forest...")
  models$RandomForest <- train(
    Outcome ~ ., data = train_data,
    method = "rf",
    trControl = cv_control,
    metric = "ROC",
    tuneGrid = expand.grid(mtry = c(10, 20, 50)),
    ntree = 500,
    importance = TRUE
  )

  # 3. SVM (Radial)
  message("Training SVM (Radial)...")
  models$SVM_Radial <- train(
    Outcome ~ ., data = train_data,
    method = "svmRadial",
    trControl = cv_control,
    metric = "ROC",
    tuneLength = 5,
    preProcess = c("center", "scale")
  )

  # 4. SVM (Linear)
  message("Training SVM (Linear)...")
  models$SVM_Linear <- train(
    Outcome ~ ., data = train_data,
    method = "svmLinear",
    trControl = cv_control,
    metric = "ROC",
    tuneGrid = expand.grid(C = c(0.1, 1, 10)),
    preProcess = c("center", "scale")
  )

  # 5. Logistic Regression with Elastic Net
  message("Training Elastic Net Logistic Regression...")
  models$ElasticNet <- train(
    Outcome ~ ., data = train_data,
    method = "glmnet",
    trControl = cv_control,
    metric = "ROC",
    tuneGrid = expand.grid(
      alpha = c(0, 0.5, 1),  # 0=Ridge, 1=Lasso, 0.5=Elastic Net
      lambda = 10^seq(-4, -1, length = 10)
    ),
    preProcess = c("center", "scale")
  )

  # 6. Gradient Boosting Machine
  message("Training GBM...")
  models$GBM <- train(
    Outcome ~ ., data = train_data,
    method = "gbm",
    trControl = cv_control,
    metric = "ROC",
    tuneGrid = expand.grid(
      n.trees = c(100, 200),
      interaction.depth = c(3, 5),
      shrinkage = c(0.01, 0.1),
      n.minobsinnode = 10
    ),
    verbose = FALSE
  )

  message("All models trained successfully!")
  return(models)
}

# 5. Evaluation Functions ----------------------------------------------
#' Evaluate Model on Test Set
#' @param model Trained caret model
#' @param test_x Test features
#' @param test_y Test labels
#' @param model_name Name of the model
#' @return List with evaluation metrics
evaluate_model <- function(model, test_x, test_y, model_name) {
  # Predictions
  pred_class <- predict(model, test_x)
  pred_prob <- predict(model, test_x, type = "prob")[, 2]

  # Confusion Matrix
  cm <- confusionMatrix(pred_class, test_y, positive = levels(test_y)[2])

  # ROC and AUC
  roc_obj <- roc(test_y, pred_prob, levels = levels(test_y))
  auc_val <- as.numeric(auc(roc_obj))

  # Calculate additional metrics
  precision <- cm$byClass["Pos Pred Value"]
  recall <- cm$byClass["Sensitivity"]
  f1 <- 2 * (precision * recall) / (precision + recall)

  list(
    model_name = model_name,
    confusion_matrix = cm,
    roc = roc_obj,
    auc = auc_val,
    accuracy = cm$overall["Accuracy"],
    sensitivity = cm$byClass["Sensitivity"],
    specificity = cm$byClass["Specificity"],
    precision = precision,
    f1_score = f1,
    predictions = pred_class,
    probabilities = pred_prob
  )
}

#' Compare All Models
#' @param models List of trained models
#' @param test_x Test features
#' @param test_y Test labels
#' @return List with all evaluation results
evaluate_all_models <- function(models, test_x, test_y) {
  results <- lapply(names(models), function(name) {
    evaluate_model(models[[name]], test_x, test_y, name)
  })
  names(results) <- names(models)
  return(results)
}

#' Create Summary Table
#' @param results List of evaluation results
#' @return Data frame with metrics
create_summary_table <- function(results) {
  metrics_df <- do.call(rbind, lapply(results, function(r) {
    data.frame(
      Model = r$model_name,
      Accuracy = round(r$accuracy, 4),
      Sensitivity = round(r$sensitivity, 4),
      Specificity = round(r$specificity, 4),
      Precision = round(r$precision, 4),
      F1_Score = round(r$f1_score, 4),
      AUC = round(r$auc, 4)
    )
  }))
  rownames(metrics_df) <- NULL
  metrics_df[order(-metrics_df$AUC), ]
}

# 6. Feature Importance Functions --------------------------------------
#' Extract Feature Importance from Models
#' @param models List of trained models
#' @param gene_names Vector of gene names
#' @return Data frame with feature importance
extract_feature_importance <- function(models, gene_names) {
  importance_list <- list()

  # XGBoost importance
  if ("XGBoost" %in% names(models)) {
    xgb_imp <- varImp(models$XGBoost)$importance
    xgb_imp$Gene <- rownames(xgb_imp)
    xgb_imp$Model <- "XGBoost"
    xgb_imp$Importance <- xgb_imp$Overall
    importance_list$XGBoost <- xgb_imp[, c("Gene", "Importance", "Model")]
  }

  # Random Forest importance
  if ("RandomForest" %in% names(models)) {
    rf_imp <- varImp(models$RandomForest)$importance
    rf_imp$Gene <- rownames(rf_imp)
    rf_imp$Model <- "RandomForest"
    rf_imp$Importance <- rf_imp$Overall
    importance_list$RandomForest <- rf_imp[, c("Gene", "Importance", "Model")]
  }

  # GBM importance
  if ("GBM" %in% names(models)) {
    gbm_imp <- varImp(models$GBM)$importance
    gbm_imp$Gene <- rownames(gbm_imp)
    gbm_imp$Model <- "GBM"
    gbm_imp$Importance <- gbm_imp$Overall
    importance_list$GBM <- gbm_imp[, c("Gene", "Importance", "Model")]
  }

  # Elastic Net importance (coefficients)
  if ("ElasticNet" %in% names(models)) {
    en_imp <- varImp(models$ElasticNet)$importance
    en_imp$Gene <- rownames(en_imp)
    en_imp$Model <- "ElasticNet"
    en_imp$Importance <- en_imp$Overall
    importance_list$ElasticNet <- en_imp[, c("Gene", "Importance", "Model")]
  }

  # Combine all
  importance_df <- do.call(rbind, importance_list)
  rownames(importance_df) <- NULL

  # Rank within each model
  importance_df <- importance_df %>%
    group_by(Model) %>%
    mutate(Rank = rank(-Importance)) %>%
    ungroup()

  return(importance_df)
}

#' Get Consensus Top Features
#' @param importance_df Feature importance data frame
#' @param top_n Number of top features per model
#' @return Data frame with consensus features
get_consensus_features <- function(importance_df, top_n = 20) {
  # Get top features from each model
  top_features <- importance_df %>%
    group_by(Model) %>%
    top_n(top_n, Importance) %>%
    ungroup()

  # Count how many models selected each gene
  consensus <- top_features %>%
    group_by(Gene) %>%
    summarise(
      n_models = n(),
      avg_importance = mean(Importance),
      avg_rank = mean(Rank),
      .groups = "drop"
    ) %>%
    arrange(-n_models, -avg_importance)

  return(consensus)
}

# 7. Main Pipeline -----------------------------------------------------
#' Run Complete ML Pipeline
#' @param seurat_obj Seurat object with platelet data
#' @param output_dir Directory for saving results
#' @return List with all results
run_ml_pipeline <- function(seurat_obj, output_dir = "results/ml_comparison") {

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("=" %>>% rep(50) %>>% paste(collapse = ""))
  message("Starting ML Pipeline")
  message("=" %>>% rep(50) %>>% paste(collapse = ""))

  # Step 1: Prepare data
  ml_data <- prepare_ml_data(seurat_obj)

  # Step 2: Split data
  data_split <- split_data(ml_data$features, ml_data$labels)

  # Prepare training data frame
  train_data <- data.frame(data_split$train_x)
  train_data$Outcome <- data_split$train_y

  # Step 3: Setup CV
  cv_control <- setup_cv_control()

  # Step 4: Train models
  message("\n--- Training Models ---")
  models <- train_all_models(train_data, cv_control)

  # Step 5: Evaluate models
  message("\n--- Evaluating Models ---")
  results <- evaluate_all_models(models, data_split$test_x, data_split$test_y)

  # Step 6: Create summary
  summary_table <- create_summary_table(results)
  message("\n=== Model Performance Summary ===")
  print(summary_table)

  # Step 7: Extract feature importance
  message("\n--- Extracting Feature Importance ---")
  importance_df <- extract_feature_importance(models, ml_data$gene_names)
  consensus_features <- get_consensus_features(importance_df)

  message("\nTop 10 Consensus Features:")
  print(head(consensus_features, 10))

  # Step 8: Save results
  message("\n--- Saving Results ---")

  # Save summary table
  write.csv(summary_table,
            file.path(output_dir, "model_comparison_summary.csv"),
            row.names = FALSE)

  # Save feature importance
  write.csv(importance_df,
            file.path(output_dir, "feature_importance_all_models.csv"),
            row.names = FALSE)

  # Save consensus features
  write.csv(consensus_features,
            file.path(output_dir, "consensus_top_features.csv"),
            row.names = FALSE)

  # Save models
  saveRDS(models, file.path(output_dir, "trained_models.rds"))

  # Save full results
  saveRDS(list(
    models = models,
    results = results,
    summary = summary_table,
    feature_importance = importance_df,
    consensus_features = consensus_features,
    data_split = data_split,
    gene_names = ml_data$gene_names
  ), file.path(output_dir, "ml_pipeline_results.rds"))

  message(sprintf("\nResults saved to: %s", output_dir))
  message("=" %>>% rep(50) %>>% paste(collapse = ""))
  message("ML Pipeline Completed Successfully!")
  message("=" %>>% rep(50) %>>% paste(collapse = ""))

  return(list(
    models = models,
    results = results,
    summary = summary_table,
    feature_importance = importance_df,
    consensus_features = consensus_features
  ))
}

# 8. Example Usage -----------------------------------------------------
if (!interactive()) {
  # Load your Seurat object
  # seurat_obj <- readRDS("path/to/Platelets_cluster_Anno.rds")
  # DefaultAssay(seurat_obj) <- "RNA"

  # Run pipeline
  # results <- run_ml_pipeline(seurat_obj)

  message("To run the pipeline, use:")
  message("  seurat_obj <- readRDS('your_data.rds')")
  message("  results <- run_ml_pipeline(seurat_obj)")
}

# Utility function for pipe operator
`%>>%` <- function(x, f) f(x)
