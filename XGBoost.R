library(reticulate)
if (!py_module_available("tensorflow")) {
  stop("Tensorflow is not available in the specified conda environment.")
}

library(Seurat)
library(keras)
library(tensorflow)
library(caret)
library(xgboost)

# Load and preprocess data
dat0 <- readRDS("Platelets_cluster_Anno.rds")
DefaultAssay(dat0) <- "RNA"
dat1 <- subset(dat0, Outcome %in% c("S", "NS"))

# Prepare dataset
data <- as.matrix(GetAssayData(dat1, slot = "data"))
labels <- dat1$Outcome
features <- t(data) # Transpose data to have genes as features
labels_factor <- factor(labels, levels = c("NS", "S"))
labels_encoded <- as.numeric(labels_factor) - 1  # Zero-index

# Split data into training and test sets
set.seed(123)
trainIndex <- createDataPartition(labels_encoded, p = 0.8, list = FALSE)
data_train <- features[trainIndex, ]
data_test <- features[-trainIndex, ]
labels_train <- labels_encoded[trainIndex]
labels_test <- labels_encoded[-trainIndex]

# Convert data to xgb.DMatrix format
dtrain <- xgb.DMatrix(data = as.matrix(data_train), label = labels_train)
dtest <- xgb.DMatrix(data = as.matrix(data_test), label = labels_test)

# Define XGBoost parameters
params <- list(
  booster = "gbtree",
  objective = "multi:softprob",
  num_class = length(unique(labels_encoded)),
  eta = 0.3,
  gamma = 0,
  max_depth = 6,
  subsample = 1,
  colsample_bytree = 1
)

# Train XGBoost model
xgb_model <- xgb.train(params = params, data = dtrain, nrounds = 100)

# Predict probabilities
xgb_pred_probs <- predict(xgb_model, dtest)
# Reshape to a matrix with two columns
xgb_pred_probs_matrix <- matrix(xgb_pred_probs, ncol = 2, byrow = TRUE)

# Extract probabilities for the positive class
positive_class_probs <- xgb_pred_probs_matrix[, 2]

# Convert probabilities to class predictions
xgb_pred_classes <- ifelse(positive_class_probs > 0.5, 1, 0)

# Convert to factors for confusion matrix
xgb_pred_classes_factor <- factor(xgb_pred_classes, levels = c(0, 1))
labels_test_factor <- factor(labels_test, levels = c(0, 1))

# Generate and print the confusion matrix
conf_matrix <- confusionMatrix(xgb_pred_classes_factor, labels_test_factor)
print(conf_matrix)
