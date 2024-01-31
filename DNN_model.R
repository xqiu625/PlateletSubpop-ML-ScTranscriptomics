library(reticulate)
py_module_available("tensorflow")

library(Seurat)
library(keras)
library(tensorflow)
library(caret)
library(xgboost)

# Load and preprocess data
dat0 <- readRDS("PlateletsOnly_cluster_Anno.rds")
DefaultAssay(dat0) <- "RNA"
dat1 <- subset(dat0, Outcome %in% c("S", "NS"))

# Prepare dataset
data <- as.matrix(GetAssayData(dat1, slot = "data"))
labels <- dat1$Outcome
features <- t(data)  # Transpose data to have genes as features
labels_factor <- factor(labels, levels = c("NS", "S"))
labels_encoded <- as.numeric(labels_factor) - 1  # Zero-index

# Split data into training and test sets
set.seed(123)
trainIndex <- createDataPartition(labels_encoded, p = 0.8, list = FALSE)
data_train <- features[trainIndex, ]
data_test <- features[-trainIndex, ]
labels_train <- labels_encoded[trainIndex]
labels_test <- labels_encoded[-trainIndex]

# Define deep neural network model
model <- keras_model_sequential() 
model %>% 
  layer_dense(units = 128, activation = 'relu', input_shape = ncol(data_train)) %>%
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 128, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = length(unique(labels_encoded)), activation = 'softmax')

# Compile model
model %>% compile(
  loss = 'sparse_categorical_crossentropy',
  optimizer = optimizer_rmsprop(lr = 0.001),
  metrics = c('accuracy')
)

# Train model
history <- model %>% fit(
  data_train, labels_train,
  epochs = 20, 
  batch_size = 50, 
  validation_split = 0.2
)

# Predict and evaluate model
predictions <- model %>% predict(data_test)
predicted_classes <- apply(predictions, 1, which.max) - 1
labels_test_factor <- factor(labels_test, levels = unique(labels_test))
predicted_classes_factor <- factor(predicted_classes, levels = levels(labels_test_factor))
conf_matrix <- confusionMatrix(predicted_classes_factor, labels_test_factor)
print(conf_matrix)
