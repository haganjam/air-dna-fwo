# Load necessary libraries
library(tidyverse)

# Parameters
set.seed(123)  # For reproducibility
distances <- c(seq(0.5, 9.9, 0.5), 10, 20, 30, 40)  # Sampling distances from the release point
species <- LETTERS[1:5]  # Species list
treatments <- c("Control", "Treatment")  # Two treatments

# Initial concentrations and decay rates for different species and treatments
species_treatment_parameters <- tibble(
  Species = rep(species, each = 2),  # Each species has two treatments
  Treatment = rep(treatments, times = length(species)),
  initial_concentration = c(1200, 500, 1000, 350, 800, 600, 1500, 400, 900, 250),  # Lower initial concentration for "Treatment"
  lambda_base = c(0.05, 0.1, 0.07, 0.12, 0.09, 0.15, 0.06, 0.11, 0.08, 0.13)  # Higher decay rates for "Treatment"
)

# Function to calculate DNA concentration with decay for each species and treatment
calculate_concentration <- function(d, initial_concentration, lambda) {
  # Exponential decay model for DNA concentration
  concentration <- initial_concentration * exp(-lambda * d)
  return(concentration)
}

# Generate dataset
data <- expand.grid(Distance = distances, Species = species, Treatment = treatments) %>%
  left_join(species_treatment_parameters, by = c("Species", "Treatment")) %>%
  mutate(
    # Calculate DNA concentration using the decay model for each species and treatment
    DNA_Concentration = mapply(calculate_concentration, Distance, initial_concentration, lambda_base)
  )

# Add some random noise to simulate environmental variation
data$DNA_Concentration <- data$DNA_Concentration + rnorm(nrow(data), mean = 0, sd = 5)

# View the simulated data
head(data)

# set the negative values to zeros
data$DNA_Concentration <- ifelse(data$DNA_Concentration < 0, 0, data$DNA_Concentration)

# Plotting DNA decay over distance for each species and treatment

# load the custom plotting theme
source("theme-minimal.R")

# check the data
head(data)

# get the range of concentrations
range(data$DNA_Concentration)
y_lims <- c(0, 1500)

# get relevant colours
jungle_palette <- c("#145A32", "#27AE60", "#8E5A30", "#6D8E6F", "#4E342E")

# plot the control treatment (downwind)
p1 <-
  ggplot(data |> dplyr::filter(Treatment == "Control"), 
         mapping = aes(x = Distance, y = DNA_Concentration, 
                       color = Species)) +
  geom_line(size = 0.75) +
  scale_color_manual(values = jungle_palette) +
  scale_y_continuous(limits = y_lims) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "lightgrey") +
  ylab("Read abundance") +
  theme_minimal_bw() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
plot(p1)

# export the plot
ggsave(filename = "02-wp-1/figures-tables/fig-1a.pdf", p1,
       width = 7, height = 7, unit = "cm")

# plot the treatment (upwind)
p2 <-
  ggplot(data |> dplyr::filter(Treatment == "Treatment"), 
       mapping = aes(x = Distance, y = DNA_Concentration, 
                     color = Species)) +
  geom_line(size = 0.75) +
  scale_color_manual(values = jungle_palette) +
  scale_y_continuous(limits = y_lims, position = "right") +  # Move y-axis to the right
  scale_x_reverse() +  # Reverse the x-axis
  geom_hline(yintercept = 0, linetype = "dashed", colour = "lightgrey") +
  ylab("Read abundance") +
  theme_minimal_bw() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
plot(p2)

# export these plots
ggsave(filename = "02-wp-1/figures-tables/fig-1b.pdf", p2,
       width = 7, height = 7, unit = "cm")






