#!/usr/bin/env Rscript

#' Platelet Classification Pipeline Using XGBoost
#' Description: Predict cell outcomes using  gene expression data

# 1. Setup and Configuration -------------------------------------------
suppressPackageMessages({
  library(xgboost)
  library(Seurat)
  library(caret)
  library(tidyverse)
  library(ROCR)
})

# Model parameters
MODEL_PARAMS <- list(
  # Data splitting
  train_ratio = 0.8,
  random_seed = 123,
  
  # XGBoost parameters
  xgb_params = list(
    booster = "gbtree",
    objective = "multi:softprob",
    eta = 0.3,
    gamma = 0,
    max_depth = 6,
    min_child_weight = 1,
    subsample = 0.8,
    colsample_bytree = 0.8,
    scale_pos_weight = 1,
    max_delta_step = 0,
    num_parallel_tree = 1,
    lambda = 1,
    alpha = 0
  ),
  
  # Training parameters
  nrounds = 100,
  early_stopping_rounds = 10,
  print_every_n = 10
)

# 2. Data Processing Functions -----------------------------------------
#' Prepare Data for XGBoost
#' @param seurat_obj Seurat object
#' @param outcome_levels Vector of outcome levels
#' @return List containing processed features and labels
prepare_data <- function(seurat_obj, outcome_levels = c("NS", "S")) {
  # Input validation
  if (!all(outcome_levels %in% unique(seurat_obj$Outcome))) {
    stop("Some outcome levels not found in data")
  }
  
  # Subset data
  data_subset <- subset(seurat_obj, Outcome %in% outcome_levels)
  
  # Extract features and labels
  features <- t(GetAssayData(data_subset, slot = "data"))
  labels <- data_subset$Outcome
  
  # Encode labels
  labels_factor <- factor(labels, levels = outcome_levels)
  labels_encoded <- as.numeric(labels_factor) - 1
  
  # Check for missing values
  if (any(is.na(features))) {
    warning("Missing values found in features matrix")
  }
  
  list(
    features = features,
    labels = labels_encoded,
    label_levels = levels(labels_factor)
  )
}

#' Split Data into Training and Test Sets
#' @param features Feature matrix
#' @param labels Encoded labels
#' @param split_ratio Training data proportion
#' @return List containing training and test data matrices
split_data <- function(features, labels, split_ratio = MODEL_PARAMS$train_ratio) {
  set.seed(MODEL_PARAMS$random_seed)
  
  # Create stratified split
  train_idx <- createDataPartition(labels, p = split_ratio, list = FALSE)
  
  # Split data
  train <- list(
    x = features[train_idx, ],
    y = labels[train_idx]
  )
  
  test <- list(
    x = features[-train_idx, ],
    y = labels[-train_idx]
  )
  
  # Create XGBoost matrices
  list(
    train = xgb.DMatrix(data = as.matrix(train$x), label = train$y),
    test = xgb.DMatrix(data = as.matrix(test$x), label = test$y),
    train_raw = train,
    test_raw = test
  )
}

# 3. Model Functions -------------------------------------------------
#' Train XGBoost Model
#' @param data_matrices List containing train and test XGBoost matrices
#' @param num_classes Number of outcome classes
#' @return Trained XGBoost model
train_model <- function(data_matrices, num_classes) {
  # Set number of classes
  params <- MODEL_PARAMS$xgb_params
  params$num_class <- num_classes
  
  # Setup watchlist for early stopping
  watchlist <- list(
    train = data_matrices$train,
    test = data_matrices$test
  )
  
  # Train model
  xgb.train(
    params = params,
    data = data_matrices$train,
    nrounds = MODEL_PARAMS$nrounds,
    watchlist = watchlist,
    early_stopping_rounds = MODEL_PARAMS$early_stopping_rounds,
    print_every_n = MODEL_PARAMS$print_every_n,
    maximize = FALSE
  )
}

#' Evaluate Model Performance
#' @param model Trained XGBoost model
#' @param data_matrices Data matrices including test data
#' @param label_levels Original label levels
#' @return List containing evaluation metrics
evaluate_model <- function(model, data_matrices, label_levels) {
  # Generate predictions
  pred_probs <- predict(model, data_matrices$test)
  pred_probs_matrix <- matrix(pred_probs, ncol = length(label_levels), byrow = TRUE)
  
  # Convert to class predictions
  pred_classes <- max.col(pred_probs_matrix) - 1
  
  # Create factors for confusion matrix
  pred_factor <- factor(pred_classes, levels = 0:(length(label_levels)-1))
  actual_factor <- factor(data_matrices$test_raw$y, levels = 0:(length(label_levels)-1))
  
  # Calculate metrics
  conf_matrix <- confusionMatrix(pred_factor, actual_factor)
  
  # For binary classification, calculate ROC and AUC
  if (length(label_levels) == 2) {
    pred <- prediction(pred_probs_matrix[,2], data_matrices$test_raw$y)
    perf <- performance(pred, "tpr", "fpr")
    auc <- performance(pred, "auc")@y.values[[1]]
  } else {
    perf <- NULL
    auc <- NULL
  }
  
  # Feature importance
  importance <- xgb.importance(
    feature_names = colnames(data_matrices$test_raw$x),
    model = model
  )
  
  list(
    confusion_matrix = conf_matrix,
    roc_curve = perf,
    auc = auc,
    feature_importance = importance,
    predictions = list(
      probabilities = pred_probs_matrix,
      classes = pred_classes
    )
  )
}

# 4. Main Execution ------------------------------------------------
main <- function() {
  # Load data
  message("Loading data...")
  seurat_obj <- readRDS("Platelets_cluster_Anno.rds")
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Prepare data
  message("Preparing data...")
  data_list <- prepare_data(seurat_obj)
  data_matrices <- split_data(data_list$features, data_list$labels)
  
  # Train model
  message("Training model...")
  model <- train_model(data_matrices, length(data_list$label_levels))
  
  # Evaluate model
  message("Evaluating model...")
  results <- evaluate_model(model, data_matrices, data_list$label_levels)
  
  # Print results
  print(results$confusion_matrix)
  if (!is.null(results$auc)) {
    message(sprintf("AUC: %.3f", results$auc))
  }
  
  # Plot feature importance
  message("Plotting feature importance...")
  png("feature_importance.png", width = 10, height = 8, units = "in", res = 300)
  xgb.plot.importance(results$feature_importance)
  dev.off()
  
  # Save results
  message("Saving results...")
  saveRDS(list(
    model = model,
    results = results,
    params = MODEL_PARAMS
  ), "platelet_xgboost_results.rds")
  
  message("Classification pipeline completed successfully!")
}

# Run pipeline
if (!interactive()) {
  main()
}
