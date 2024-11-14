# Load required libraries
library(sp)
library(dplyr)
library(tidyr)
library(spdep)  # For spatial modeling
library(nlme)   # For GLS models
library(ggplot2)

# Define grid size
grid_size <- 20  # 20 x 20 grid (each cell is 1 km^2)

# Generate grid cell IDs and coordinates
grid <- expand.grid(x = 1:grid_size, y = 1:grid_size) %>%
  mutate(cell_id = 1:n())

# Simulate covariates and eDNA concentration
grid <- grid %>%
  mutate(
    VegDensity = runif(n(), 0.1, 1),          # Vegetation density (0.1 - 1)
    Temperature = runif(n(), 10, 30),         # Temperature in Celsius
    Humidity = runif(n(), 0.3, 1),            # Relative humidity (0.3 - 1)
    WindSpeed = runif(n(), 0.5, 2),           # Wind speed (0.5 - 2 m/s)
    WindDir = runif(n(), 270, 360),             # Wind direction in degrees
  )

# Simulate eDNA concentration
grid$eDNA_obs <- with(grid,
                      rnorm(n = nrow(grid), 
                            mean = 1 + 0.5*VegDensity + 0.25*Humidity + 0.5*Temperature, sd = 1.5))

head(grid)

# Compute pairwise distances between grid cells
dist_matrix <- as.matrix(dist(grid[, c("x", "y")]))  # Euclidean distance between cell centers

# Define wind influence function
wind_influence <- function(wind_speed, wind_dir, angle_diff, alpha = 0.5) {
  exp(-angle_diff / alpha) * wind_speed
}

# Calculate wind influence between cells based on direction and distance
calculate_wind_factor <- function(grid, alpha = 0.5) {
  wind_matrix <- matrix(0, nrow = nrow(grid), ncol = nrow(grid))
  for (i in 1:nrow(grid)) {
    for (j in 1:nrow(grid)) {
      if (i != j) {
        # Calculate angle between cells
        dx <- grid$x[j] - grid$x[i]
        dy <- grid$y[j] - grid$y[i]
        cell_angle <- atan2(dy, dx) * (180 / pi)
        
        # Calculate angle difference with wind direction
        angle_diff <- abs(grid$WindDir[i] - cell_angle)
        angle_diff <- min(angle_diff, 360 - angle_diff)
        
        # Calculate wind influence
        wind_matrix[i, j] <- wind_influence(grid$WindSpeed[i], grid$WindDir[i], angle_diff, alpha)
      }
    }
  }
  return(wind_matrix)
}

# Generate wind influence matrix
wind_matrix <- calculate_wind_factor(grid)

# Define GLS model with covariates and distance-based eDNA decay
gls_model <- gls(
  eDNA_obs ~ VegDensity + Temperature + Humidity,
  data = grid,
  correlation = corSpatial(form = ~x + y, type = "exponential", nugget = TRUE),
  method = "REML"
)

summary(gls_model)

# Define decay function
decay_factor <- function(distance, lambda = 5) {
  exp(-distance / lambda)
}

# Calculate transported eDNA from other cells with decay and wind influence
calculate_transported_eDNA <- function(grid, dist_matrix, wind_matrix, lambda = 5) {
  transported_eDNA <- rep(0, nrow(grid))
  for (j in 1:nrow(grid)) {
    for (i in 1:nrow(grid)) {
      if (i != j) {
        transported_eDNA[j] <- transported_eDNA[j] +
          decay_factor(dist_matrix[i, j], lambda) * wind_matrix[i, j] * grid$eDNA_obs[i]
      }
    }
  }
  return(transported_eDNA)
}

# Calculate transported eDNA contribution
grid$transported_eDNA <- calculate_transported_eDNA(grid, dist_matrix, wind_matrix)

# Update GLS model to include transported eDNA
gls_model_with_transport <- gls(
  eDNA_obs ~ VegDensity + Temperature + Humidity + transported_eDNA,
  data = grid,
  correlation = corSpatial(form = ~x + y, type = "exponential", nugget = TRUE),
  method = "REML"
)

# Check the model output
summary(gls_model_with_transport)

# Assess residuals
residuals <- resid(gls_model_with_transport)
plot(residuals, main = "Residuals of the GLS Model with Transported eDNA")
hist(residuals, main = "Residual Distribution")

# plot observed versus predicted
plot(grid$eDNA_obs, predict(gls_model_with_transport))

