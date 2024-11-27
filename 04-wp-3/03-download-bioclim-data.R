
# download the bioclim data

# load the geodata package
library(geodata)

# download WorldClim bioclimatic data at a resolution of 2.5 minutes
worldclim_global(var = "bio", res = 2.5, path = "03-wp-3/raw-data/bio-clim-data/")