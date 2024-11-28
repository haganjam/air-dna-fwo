
# reproject the Oesner et al. (2019) lynx SDM

# load relevant libraries
library(raster)
library(terra)

# load the raster
s_pres <- rast("04-wp-3/raw-data/global_ensemble_scale_integrated.tif")

# set new crs
new_crs <- "EPSG:4326"

# reproject to the relevant crs
s_pres <- project(s_pres, new_crs)

# check the data
plot(s_pres)

# save the data for later use
terra::writeRaster(s_pres , "04-wp-3/data/oesner-2019-sdm.tif", overwrite = TRUE)
