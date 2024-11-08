#!/usr/bin/env Rscript

#' Platelet Classification Pipeline Using Deep Learning
#' Description: Classify cell outcomes based on gene expression


# 1. Setup and Configuration -------------------------------------------
suppressPackageMessages({
  library(reticulate)
  library(Seurat)
  library(keras)
  library(tensorflow)
  library(caret)
  library(xgboost)
  library(tidyverse)
})

# Verify TensorFlow availability
if (!py_module_available("tensorflow")) {
  stop("TensorFlow is not available. Please install tensorflow-cpu or tensorflow-gpu.")
}

# Model parameters
MODEL_PARAMS <- list(
  hidden_units = 128,
  dropout_rates = c(0.4, 0.3),
  learning_rate = 0.001,
  epochs = 20,
  batch_size = 50,
  validation_split = 0.2,
  train_split = 0.8,
  random_seed = 123
)

# 2. Data Processing Functions -----------------------------------------
#' Prepare Data for Classification
#' @param seurat_obj Seurat object
#' @param outcome_levels Vector of outcome levels
#' @return List containing processed features and labels
prepare_data <- function(seurat_obj, outcome_levels = c("NS", "S")) {
  # Subset data
  data_subset <- subset(seurat_obj, Outcome %in% outcome_levels)
  
  # Extract features and labels
  features <- t(GetAssayData(data_subset, slot = "data"))
  labels <- data_subset$Outcome
  
  # Encode labels
  labels_factor <- factor(labels, levels = outcome_levels)
  labels_encoded <- as.numeric(labels_factor) - 1
  
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
#' @return List containing training and test sets
split_data <- function(features, labels, split_ratio = MODEL_PARAMS$train_split) {
  set.seed(MODEL_PARAMS$random_seed)
  
  train_idx <- createDataPartition(labels, p = split_ratio, list = FALSE)
  
  list(
    train = list(
      x = features[train_idx, ],
      y = labels[train_idx]
    ),
    test = list(
      x = features[-train_idx, ],
      y = labels[-train_idx]
    )
  )
}

# 3. Model Functions -------------------------------------------------
#' Create Neural Network Model
#' @param input_shape Input feature dimension
#' @param n_classes Number of output classes
#' @return Keras model
create_model <- function(input_shape, n_classes) {
  model <- keras_model_sequential()
  
  # First layer block
  model %>%
    layer_dense(
      units = MODEL_PARAMS$hidden_units, 
      activation = 'relu',
      input_shape = input_shape
    ) %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = MODEL_PARAMS$dropout_rates[1])
  
  # Second layer block
  model %>%
    layer_dense(
      units = MODEL_PARAMS$hidden_units, 
      activation = 'relu'
    ) %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = MODEL_PARAMS$dropout_rates[2])
  
  # Output layer
  model %>%
    layer_dense(
      units = n_classes, 
      activation = 'softmax'
    )
  
  # Compile model
  model %>% compile(
    loss = 'sparse_categorical_crossentropy',
    optimizer = optimizer_rmsprop(
      learning_rate = MODEL_PARAMS$learning_rate
    ),
    metrics = c('accuracy', 'AUC')
  )
  
  return(model)
}

#' Train Model with Callbacks
#' @param model Keras model
#' @param train_data Training data list
#' @param validation_split Validation split ratio
#' @return Training history
train_model <- function(model, train_data, 
                       validation_split = MODEL_PARAMS$validation_split) {
  # Define callbacks
  callbacks <- list(
    callback_early_stopping(
      monitor = "val_loss",
      patience = 5,
      restore_best_weights = TRUE
    ),
    callback_reduce_lr_on_plateau(
      monitor = "val_loss",
      factor = 0.5,
      patience = 3,
      min_lr = 1e-6
    )
  )
  
  # Train model
  history <- model %>% fit(
    x = train_data$x,
    y = train_data$y,
    epochs = MODEL_PARAMS$epochs,
    batch_size = MODEL_PARAMS$batch_size,
    validation_split = validation_split,
    callbacks = callbacks,
    verbose = 1
  )
  
  return(history)
}

#' Evaluate Model Performance
#' @param model Trained model
#' @param test_data Test data list
#' @param label_levels Original label levels
#' @return List containing evaluation metrics
evaluate_model <- function(model, test_data, label_levels) {
  # Generate predictions
  predictions <- model %>% predict(test_data$x)
  predicted_classes <- apply(predictions, 1, which.max) - 1
  
  # Convert to factors for confusion matrix
  actual <- factor(test_data$y, levels = 0:(length(label_levels)-1))
  predicted <- factor(predicted_classes, levels = levels(actual))
  
  # Calculate metrics
  conf_matrix <- confusionMatrix(predicted, actual)
  
  # Calculate ROC curve if binary classification
  if (length(label_levels) == 2) {
    roc_obj <- roc(test_data$y, predictions[,2])
    auc_value <- auc(roc_obj)
  } else {
    roc_obj <- NULL
    auc_value <- NULL
  }
  
  list(
    confusion_matrix = conf_matrix,
    roc = roc_obj,
    auc = auc_value,
    predictions = predictions
  )
}

# 4. Main Execution ------------------------------------------------
main <- function() {
  # Load data
  message("Loading data...")
  seurat_obj <- readRDS("PlateletsOnly_cluster_Anno.rds")
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Prepare data
  message("Preparing data...")
  data_list <- prepare_data(seurat_obj)
  data_split <- split_data(data_list$features, data_list$labels)
  
  # Create and train model
  message("Creating and training model...")
  model <- create_model(
    input_shape = ncol(data_split$train$x),
    n_classes = length(unique(data_list$labels))
  )
  
  history <- train_model(model, data_split$train)
  
  # Evaluate model
  message("Evaluating model...")
  results <- evaluate_model(
    model, 
    data_split$test,
    data_list$label_levels
  )
  
  # Print results
  print(results$confusion_matrix)
  if (!is.null(results$auc)) {
    message(sprintf("AUC: %.3f", results$auc))
  }
  
  # Save model and results
  message("Saving results...")
  save_model_hdf5(model, "platelet_classifier_model.h5")
  saveRDS(results, "platelet_classifier_results.rds")
  
  message("Classification pipeline completed successfully!")
}

# Run pipeline
if (!interactive()) {
  main()
}
