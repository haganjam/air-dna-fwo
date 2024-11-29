
# load the Oesner et al. (2019) SDM
# https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13784#support-information-section

# load necessary libraries
library(sf)
library(dplyr)
library(rnaturalearth)
library(biomod2)
library(raster)
library(spData)
library(terra)

## set-up the spatial extent data for BE, NL, FR, L and DE (to clip the bioclim data)

# load world map with country polygons
world <- ne_countries(scale = "large", returnclass = "sf")

# filter the countries to include only Belgium, Germany, and France
selected_countries <-
  world %>%
  filter(admin %in% c("Belgium", "Germany", "France", "Netherlands",
                      "Luxembourg") & region_un == "Europe")

# define the bounding box for Europe (adjust as necessary)
# here we use approximate bounds for Western Europe
europe_bbox <- st_bbox(c(xmin = -10, ymin = 40, xmax = 20, ymax = 60), crs = st_crs(selected_countries))

# convert the bounding box to an sf polygon for clipping
europe_bbox_polygon <- st_as_sfc(europe_bbox)

# clip the selected countries to the European bounding box
european_territories <- st_intersection(selected_countries, europe_bbox_polygon)

# calculate the area of each geometry and filter to keep only the largest polygon for each country
# removes territories that are not part of the mainland
european_territories <-
  european_territories %>%
  st_cast("POLYGON") %>%           # ensure each polygon is treated separately
  mutate(area = st_area(.)) %>%    # calculate the area of each polygon
  group_by(admin) %>%              # group by country
  filter(area == max(area)) %>%    # keep only the largest polygon for each country
  ungroup() %>%                    # ungroup after filtering
  dplyr::select(-area)                    # remove the temporary area column

# combine into a single multipolygon if needed
combined_polygon <- st_union(european_territories)

# display as single polygons
plot(st_geometry(european_territories))

# display the clipped polygon
plot(st_geometry(combined_polygon), col = 'lightblue')


## Oesner 2019 SDM

# load the raster
s_pres <- rast("04-wp-3/data/oesner-2019-sdm.tif")

# crop the bioclimatic data to the bounding box of the selected countries
s_pres <- crop(s_pres, extent(st_bbox(european_territories)))

# mask the bioclim data to the exact boundaries of the selected countries
s_pres <- mask(s_pres, european_territories)

# plot one of the bioclimatic layers as an example
plot(s_pres)

# plot the figure

# might need to downscale this raster
s_pres <- terra::aggregate(s_pres, fact = 4)

# create a raster data.frame for plotting
# convert to a tibble
plot_pred <- as.data.frame(s_pres, xy = TRUE) |> dplyr::as_tibble()

# normalise the probability of occurrence
plot_pred <-
  plot_pred |>
  dplyr::mutate(prod = prod / 1000)

# create a polygon over Belgium
polygon <- dplyr::tibble(xmin = 0.5,
                         xmax = 7.5,
                         ymin = 48.5,
                         ymax = 52.5)

# load ggplot2 for plotting
library(ggplot2)

# plot the data
main <-
  ggplot() +
  geom_raster(data = plot_pred, 
              mapping = aes(x = x, y = y, fill = prod)) +
  geom_sf(data = st_geometry(european_territories), fill = NA, colour = "lightgrey") +
  geom_rect(data = polygon,
            mapping = aes(xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax), colour = "red", fill = NA) +
  scale_fill_viridis_c(option = "E", end = 0.95, begin = 0.05) +
  scale_x_continuous(limits = c(-4.9, 15.05), expand = c(0, 0)) +
  scale_y_continuous(limits = c(42, 54.5), expand = c(0, 0)) +
  guides(fill = guide_colorbar(barheight = 0.5,
                               barwidth = 15,
                               ticks.colour = NA,
                               title = "Occurrence probability")) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(vjust = 0.9, colour = "black", size = 8),
        legend.text = element_text(colour = "black", size = 8),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
main

# extract the legend
main_leg <- ggpubr::get_legend(main)

# remove legend from the main plot
main <- main + theme(legend.position = "none")

inset <-
  ggplot() +
  geom_raster(data = plot_pred, 
              mapping = aes(x = x, y = y, fill = prod)) +
  geom_sf(data = st_geometry(european_territories), fill = NA, colour = "lightgrey") +
  scale_fill_viridis_c(option = "E", end = 0.95, begin = 0.05) +
  scale_x_continuous(limits = c(0.5, 7.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(48.5, 52.5), expand = c(0, 0)) +
  xlab("Longitude (dd)") +
  ylab("Latitude (dd)") +
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "red", size = 1),
        panel.background = element_rect(fill = "white"),
        axis.text = element_text(colour = "black", size = 7.5),
        axis.title.x = element_text(vjust = -0.2, size = 7.5),
        axis.title.y = element_text(vjust = 2, size = 7.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
inset

# load patchwork
library(patchwork)

# combine the main and the inset
p1 <- (main | inset)

# combine with the legend
p2 <- (p1 / main_leg) + plot_layout(heights = c(1.5, 0.1))

# export as a png file
ggsave(filename = "04-wp-3/figures-tables/fig-1.tiff", p2,
       width = 18, height = 10, units = "cm", dpi = 1000)









