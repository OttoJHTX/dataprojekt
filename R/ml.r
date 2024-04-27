library(tensorflow)
library(keras)

num_rows <- 100000
num_cols <- 10

# Generate random data with a normal distribution

define_model <- function() {
  # Define the model architecture
  model <- keras_model_sequential() %>%
    layer_dense(units = 10, activation = "relu", input_shape = c(10), name = "layer1") %>%
    layer_dense(units = 40, activation = "relu", name = "layer2") %>%
    layer_dense(units = 40, activation = "relu", name = "layer3") %>%
    layer_dense(units = 3, activation = "softmax", name = "layer4")
  
  model %>% compile(
    optimizer = 'adam',
    loss = 'categorical_crossentropy',
    metrics = list('accuracy')
  )
  
  return (model)
}

train_model <- function(model, data, labels, epochs, batch_size) {
  model %>% fit(
    data,
    labels,
    epochs = epochs,
    batch_size = batch_size
  )
  
  return (model)
}

data <- matrix(rnorm(num_rows * num_cols), nrow = num_rows, ncol = num_cols)
labels <- matrix(c(0, 0, 0), nrow = num_rows, ncol = 3)

for (i in 1:nrow(labels)) {
  index <- sample(1:3, 1)
  labels[i,index] <- 1
}

network <- define_model()
network <- train_model(network, data, labels, 50, 32)

network$predict(matrix(rnorm(1 * num_cols), nrow = 1, ncol = num_cols))





