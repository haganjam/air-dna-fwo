
# plot a probability of occupancy distribution

# sample from the beta distribution
x <- rbeta(n = 10000, shape1 = 5, shape2 = 10)

# load ggplot2
library(ggplot2)

# load relevant theme
source("theme-minimal.R")

p1 <-
  ggplot(data = dplyr::tibble(v1 = x),
       mapping = aes(x = v1)) +
  geom_density(fill = "#27ae60", bw = 0.04, alpha = 0.6) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5) ) +
  ylab(NULL) +
  xlab(NULL) +
  theme_minimal_bw() +
  theme(panel.border = element_rect(colour = "black", size = 1, fill = NA),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# export these plots
ggsave(filename = "03-wp-2/figures-tables/fig-density.pdf", p1,
       width = 6, height = 6, unit = "cm")



sp_site <- 
  dplyr::tibble(sp1 = c(1, 0, 0, 0),
                sp2 = c(1, 1, 0, 0),
                sp3 = c(0, 1, 1, 1),
                sp4 = c(1, 0, 0, 1))

vegan::vegdist(x = sp_site, method = "bray", diag = TRUE)
vegan::veg




