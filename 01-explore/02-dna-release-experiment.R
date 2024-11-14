
# Load necessary libraries
library(tidyverse)

# Parameters
set.seed(123)  # For reproducibility
initial_concentration <- 1000      # Initial DNA concentration at the release point
distances <- c(0.5, 1, 5, 10, 20, 30, 40)  # Sampling distances from the release point
times <- seq(0, 5, by = 0.5)       # Sampling times in hours
lambda_base <- 0.07                # Base spatial decay rate (distance decay)
mu_base <- 0.4                     # Base temporal decay rate (time decay)

# Function to calculate DNA concentration with decay
calculate_concentration <- function(d, t) {
  # Exponential decay model for DNA concentration
  concentration <- initial_concentration * exp(-lambda_base * d) * exp(-mu_base * t)
  return(concentration)
}

# Generate dataset
data <- expand.grid(Distance = distances, Time = times) %>%
  mutate(
    # Calculate DNA concentration using decay model
    DNA_Concentration = mapply(calculate_concentration, Distance, Time)
  )

# Add some random noise to simulate environmental variation
data$DNA_Concentration <- data$DNA_Concentration + rnorm(nrow(data), mean = 0, sd = 5)

# View the simulated data
head(data)

# Plotting DNA decay over time for each distance
ggplot(data, aes(x = Time, y = DNA_Concentration, color = as.factor(Distance))) +
  geom_line() +
  labs(title = "Simulated DNA Decay Over Time at Different Distances from Release Point",
       x = "Time (hours)", y = "DNA Concentration",
       color = "Distance (m)") +
  theme_minimal()

# Plotting DNA decay over distance for each time point
ggplot(data, aes(x = Distance, y = DNA_Concentration, color = as.factor(Time))) +
  geom_line() +
  labs(title = "Simulated DNA Decay Over Distance at Different Times from Release",
       x = "Distance (m)", y = "DNA Concentration",
       color = "Time (hours)") +
  theme_minimal()






