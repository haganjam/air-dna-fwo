
# make a map

# load relevant packages
library(raster)
library(sf)
library(ggplot2)
library(prismatic)
library(ggspatial)

# Load the GeoTIFF exported from GEE
satellite_image <- stack("02-wp-1-2/data/site-wp2.tif")

# Convert raster to matrix for ggplot2 compatibility
satellite_matrix <- as.data.frame(satellite_image, xy = TRUE)

# rename the variables
satellite_matrix <- 
  satellite_matrix |>
  dplyr::rename(red = vis.red,   #Rename bands
                green = vis.green,
                blue = vis.blue) |>
  dplyr::filter(red != 0)

# Define the coordinates of the sampling points (same as in your GEE code)
coordinates <- data.frame(
  lon = c(5.59874, 5.58452, 5.65584, 5.64162, 5.57054, 5.64187, 5.62765, 5.65634, 
          5.62789, 5.64236, 5.61391, 5.64261, 5.61415, 5.62862, 5.62886, 5.65758, 
          5.60064, 5.58640, 5.57217, 5.57240),
  lat = c(50.90294, 50.90309, 50.91132, 50.91147, 50.91222, 50.92046, 50.92061, 50.92929, 
          50.92960, 50.93844, 50.93874, 50.94742, 50.94773, 50.95656, 50.96555, 50.97423, 
          50.97484, 50.97499, 50.97514, 50.98412)
)

# Convert to an sf object for easier handling with ggplot2
points_sf <- st_as_sf(coordinates, coords = c("lon", "lat"), crs = 4326)

# check the pairwise distances
mean(st_distance(points_sf))

# extract the bounding box
bound <- satellite_image@extent
y_lims <- c(bound[3], bound[4])
x_lims <- c(bound[1], bound[2])

# modify the x min
x_lims[1] <- 5.56

# make the plot
p1 <-
  ggplot() +
  geom_raster(data = satellite_matrix, aes(x = x, y = y),
              fill = rgb(r = satellite_matrix$red,
                         g = satellite_matrix$green,
                         b = satellite_matrix$blue,
                         maxColorValue = 255), 
              show.legend = FALSE) +
  scale_fill_identity() + 
  scale_x_continuous(limits = x_lims, expand = c(0, 0)) +
  scale_y_continuous(limits = y_lims, expand = c(0, 0)) +
  geom_spatial_point(data = coordinates, 
                     mapping = aes(x = lon, y = lat), crs = 4326, size = 3, colour = "gold") +
  coord_equal(clip = "on") +
  annotation_north_arrow(location = "tr", which_north = "true", pad_x = unit(0.05, "in"), 
                         pad_y = unit(0.05, "in"), style = north_arrow_fancy_orienteering,
                         text_col = "white") +
  xlab("Longitude (dd)") +
  ylab("Latitude (dd)") +
  ggtitle("Sampling sites: Nationaal Park Hoge Kempen") +
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        panel.background = element_rect(fill = "white"),
        axis.text = element_text(colour = "black", size = 8),
        axis.title.x = element_text(vjust = -0.2, size = 8),
        axis.title.y = element_text(vjust = 2, size = 8),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.title = element_text(hjust = 0.5, size = 10))

# export as a png file
ggsave(filename = "03-wp-2/figures-tables/fig-3a.tiff", p1,
       width = 11, height = 14, units = "cm", dpi = 1000)




