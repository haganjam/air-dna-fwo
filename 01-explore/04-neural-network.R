# Load necessary libraries
library(torch)
library(dplyr)
library(ggplot2)

# Set seed for reproducibility
set.seed(42)

# Simulation parameters
num_sites <- 100
num_species <- 10
embedding_dim <- 3  # Embedding dimensionality
lambda_reg <- 0.01  # Regularization strength

# Simulate occupancy probabilities (X)
occupancy_probs <- matrix(runif(num_sites * num_species, 0, 1), nrow = num_sites, ncol = num_species)

# Simulate true presence/absence (Y) using a threshold on occupancy_probs
true_presence_absence <- ifelse(occupancy_probs > 0.5, 1, 0)

# Simulate a trait similarity matrix (D)
trait_distances <- matrix(runif(num_species^2, 0, 1), nrow = num_species, ncol = num_species)
trait_distances <- (trait_distances + t(trait_distances)) / 2  # Make it symmetric
diag(trait_distances) <- 0  # No self-distance

# Normalize the trait distances for stability
trait_distances <- trait_distances / max(trait_distances)

# Convert data to torch tensors
X <- torch_tensor(occupancy_probs, dtype = torch_float())
Y <- torch_tensor(true_presence_absence, dtype = torch_float())
D <- torch_tensor(trait_distances, dtype = torch_float())

# Define the neural network
net <- nn_module(
  initialize = function() {
    self$fc1 <- nn_linear(num_species, 32)
    self$fc2 <- nn_linear(32, 16)
    self$output <- nn_linear(16, num_species)
    
    # Learnable embeddings for species
    self$species_embeddings <- nn_parameter(torch_randn(num_species, embedding_dim))
  },
  
  forward = function(x) {
    x <- torch_relu(self$fc1(x))
    x <- torch_relu(self$fc2(x))
    x <- torch_sigmoid(self$output(x))
    return(x)
  },
  
  # Custom loss function
  calculate_loss = function(preds, targets) {
    bce_loss <- nnf_binary_cross_entropy(preds, targets)
    
    # Trait regularization loss
    embedding_diff <- self$species_embeddings$unsqueeze(1) - self$species_embeddings
    pairwise_distances <- torch_sum(embedding_diff^2, dim = -1)
    trait_loss <- torch_sum(D * pairwise_distances)
    
    total_loss <- bce_loss + lambda_reg * trait_loss
    return(total_loss)
  }
)

# Instantiate the model
model <- net()

# Define optimizer
optimizer <- optim_adam(model$parameters, lr = 0.01)

# Training loop with loss tracking
num_epochs <- 100
batch_size <- 16
losses <- numeric(num_epochs)

for (epoch in 1:num_epochs) {
  model$train()
  
  # Forward pass
  preds <- model(X)
  
  # Calculate loss
  loss <- model$calculate_loss(preds, Y)
  losses[epoch] <- loss$item()
  
  # Backward pass and optimization
  optimizer$zero_grad()
  loss$backward()
  optimizer$step()
  
  # Print loss every 10 epochs
  if (epoch %% 10 == 0) {
    cat(sprintf("Epoch %d, Loss: %.4f\n", epoch, loss$item()))
  }
}

# Evaluate the model
model$eval()
predictions <- as.array(model(X))
as.vector(predictions)

# Convert true presence/absence and predictions to data frames for plotting
true_df <- data.frame(
  site = rep(1:num_sites, num_species),
  species = rep(1:num_species, each = num_sites),
  true_presence = as.vector(true_presence_absence)
)

pred_df <- data.frame(
  site = rep(1:num_sites, num_species),
  species = rep(1:num_species, each = num_sites),
  predicted_presence = as.vector(predictions)
)

# observed versus predicted
plot(true_df$true_presence, pred_df$predicted_presence)
plot(pred_df$predicted_presence, true_df$true_presence)

# True vs. Predicted plot (first species for simplicity)
ggplot() +
  geom_point(data = true_df, aes(x = site, y = true_presence), color = "red", size = 2) +
  geom_point(data = pred_df, aes(x = site, y = predicted_presence), color = "blue", size = 2) +
  labs(title = "True vs Predicted Presence for Species 1", x = "Site", y = "Presence Probability") +
  theme_minimal()

