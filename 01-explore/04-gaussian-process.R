install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
install.packages("StanHeaders", repos = "https://cloud.r-project.org/")

library(rstan)
example(stan_model, package = "rstan", run.dontrun = TRUE)


# Load necessary libraries
library(brms)
library(MASS)
library(ggplot2)

# Step 1: Simulate Site Locations and Asymmetrical Distance Matrix
set.seed(123)
n_sites <- 20

# Generate random site coordinates
site_coords <- data.frame(
  site = 1:n_sites,
  x = runif(n_sites, 0, 10),
  y = runif(n_sites, 0, 10)
)

# Compute base distance matrix (Euclidean distances)
euclidean_dist <- as.matrix(dist(site_coords[, c("x", "y")]))

# Add asymmetry to the distance matrix based on prevailing wind direction
# Assume a wind bias favoring movement in +x direction
wind_bias <- outer(site_coords$x, site_coords$x, "-") # Difference in x-coordinates
asym_dist <- euclidean_dist + 0.5 * (wind_bias > 0) - 0.5 * (wind_bias < 0)

# Step 2: Define the Gaussian Process Kernel with a Nugget
gp_kernel <- function(D, sigma, length_scale) {
  sigma^2 * exp(-D^2 / (2 * length_scale^2))
}

sigma <- 1         # Variance of the GP
length_scale <- 2  # Length scale of the GP

# Compute the covariance matrix
cov_matrix <- gp_kernel(asym_dist, sigma, length_scale)

# Ensure positive definiteness by adding a nugget
cov_matrix <- cov_matrix + diag(1e-6, n_sites)

# Step 3: Simulate Latent Probabilities and Binary Outcomes
latent_probs <- mvrnorm(1, mu = rep(0, n_sites), Sigma = cov_matrix, tol = 0.03) # Simulate latent function
probabilities <- pnorm(latent_probs) # Convert to probabilities using probit link
detections <- rbinom(n_sites, size = 1, prob = probabilities) # Simulate binary detections

# Combine data
data <- data.frame(
  site = site_coords$site,
  x = site_coords$x,
  y = site_coords$y,
  detection = detections
)

# Step 4: Fit the Gaussian Process Classification Model
gp_model <- brm(
  detection ~ gp(x, y),  # Gaussian Process on spatial coordinates
  data = data,
  family = bernoulli(link = "logit"), # Binary outcome with probit link
  chains = 4,
  iter = 2000,
  control = list(adapt_delta = 0.95)
)

# Step 5: Examine Results
summary(gp_model)
plot(gp_model)

# Step 6: Visualize Predictions
# Create a grid of points for prediction
x_seq <- seq(0, 10, length.out = 50)
y_seq <- seq(0, 10, length.out = 50)
prediction_grid <- expand.grid(x = x_seq, y = y_seq)

# Make predictions
predictions <- posterior_epred(gp_model, newdata = prediction_grid)
mean_probs <- rowMeans(predictions)

# Add predictions to the grid
prediction_grid$mean_prob <- mean_probs

# Plot the spatial predictions
ggplot(prediction_grid, aes(x, y, fill = mean_prob)) +
  geom_tile() +
  geom_point(data = data, aes(x = x, y = y, color = as.factor(detection)), size = 3) +
  scale_fill_viridis_c(name = "Predicted Probability") +
  scale_color_manual(name = "Detection", values = c("0" = "red", "1" = "blue")) +
  labs(title = "Predicted Detection Probabilities with Observed Detections") +
  theme_minimal()
