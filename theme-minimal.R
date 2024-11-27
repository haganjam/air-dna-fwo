
# custom plotting theme

# load the ggplot2 library
library(ggplot2)

theme_minimal_bw <- function() {
  theme_minimal() +  # Start with minimal theme
    theme(
      plot.background = element_rect(colour = "white"),
      panel.background = element_rect(colour = "white"),
      panel.grid = element_blank(),          # Remove grid lines
      panel.border = element_blank(),        # Remove the full box border
      axis.line = element_line(color = "black"), # Add axis lines only
      text = element_text(color = "black"),  # All text in black
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(color = "black", hjust = 0.5), # Center plot title
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black")
    )
}
