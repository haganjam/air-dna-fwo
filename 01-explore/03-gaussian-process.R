
# Load necessary libraries
library(tidyverse)
library(sp)
library(raster)
library(MASS)

## simulate two-layer observation-detection model

# Set simulation parameters
set.seed(123)
n_sites <- 100        # Number of sampling sites
true_sites <- 20      # Sites with experimentally placed species A
area <- 10            # Landscape area in km^2 (assumed square)
wind_direction <- pi / 4  # Wind direction in radians (e.g., NE)

# Generate random site locations
site_coords <- data.frame(
  x = runif(n_sites, 0, sqrt(area)),
  y = runif(n_sites, 0, sqrt(area))
)

# Assign true presence
true_presence <- rep(0, n_sites)
true_indices <- sample(1:n_sites, true_sites)  # Choose 20 sites for true presence
true_presence[true_indices] <- 1  # Mark as true presence

# Wind kernel function (anisotropic Gaussian)
wind_kernel <- function(dx, dy, wind_dir, sigma_x = 0.5, sigma_y = 0.2) {
  x_rot <- dx * cos(wind_dir) + dy * sin(wind_dir)
  y_rot <- -dx * sin(wind_dir) + dy * cos(wind_dir)
  exp(- (x_rot^2 / (2 * sigma_x^2) + y_rot^2 / (2 * sigma_y^2)))
}

# Compute detection probabilities influenced by wind
detection_prob <- numeric(n_sites)
for (i in 1:n_sites) {
  dx <- site_coords$x[i] - site_coords$x
  dy <- site_coords$y[i] - site_coords$y
  wind_effect <- sapply(1:n_sites, function(j) {
    wind_kernel(dx[j], dy[j], wind_direction)
  })
  detection_prob[i] <- 0.7 * true_presence[i] + 0.3 * sum(true_presence * wind_effect) / sum(wind_effect)
}

# Simulate observed detections
observed_detections <- rbinom(n_sites, 1, detection_prob)

# Combine all data into a single data frame
data <- data.frame(
  x = site_coords$x,
  y = site_coords$y,
  true_presence = true_presence,
  detection_prob = detection_prob,
  observed_detections = observed_detections
)

# Plot true presence and observed detections
ggplot(data, aes(x = x, y = y)) +
  geom_point(aes(color = as.factor(true_presence)), size = 3) +
  geom_point(aes(shape = as.factor(observed_detections)), size = 2) +
  scale_color_manual(values = c("red", "blue"), name = "True Presence") +
  scale_shape_manual(values = c(1, 19), name = "Observed Detection") +
  theme_minimal() +
  labs(title = "True Presence and Observed Detections")


## fit the stan model

# Create wind kernel matrix
wind_matrix <- matrix(0, nrow = n_sites, ncol = n_sites)
for (i in 1:n_sites) {
  for (j in 1:n_sites) {
    dx <- site_coords$x[j] - site_coords$x[i]
    dy <- site_coords$y[j] - site_coords$y[i]
    wind_matrix[i, j] <- wind_kernel(dx, dy, wind_direction)
  }
}

# Package data for Stan
stan_data <- list(
  n_sites = n_sites,
  observed_detections = observed_detections,
  wind_effect = wind_matrix
)

# load the rstan library
library(rstan)

# Fit the model
fit <- stan(
  file = "01-explore/wind_model.stan",
  data = stan_data,
  iter = 1000,
  chains = 4,
  seed = 123
)

traceplot(fit)


