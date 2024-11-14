
# species distribution model
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14101

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
world <- ne_countries(scale = "medium", returnclass = "sf")

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


## bioclimatic data from WorldClim

# load the predictors from .tif files
path <- list.files("03-wp-3/raw-data/bio-clim-data/wc2.1_2.5m/",
                   pattern = ".tif$", 
                   full.names = TRUE)

# we use stack to make sure they are grouped and then we have the same resolution, extent etc
# this allows us to deal with all the rasters at the same time
s_pres <- raster::stack(path)
crs(s_pres) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'

# crop the bioclimatic data to the bounding box of the selected countries
s_pres <- crop(s_pres, extent(st_bbox(european_territories)))

# mask the bioclim data to the exact boundaries of the selected countries
s_pres <- mask(s_pres, as(european_territories, "Spatial"))

# plot one of the bioclimatic layers as an example
plot(s_pres[[1]], main = "Bioclimatic Variable 1 (Masked to Selected Countries)")

# convert to a raster stack
s_pres <- raster::stack(s_pres)

## load the occurrence data

# load the occurrence data
lynx_dat <- readr::read_csv("03-wp-3/data/lynx_lynx_occ.csv")
head(lynx_dat)

# extract the relevant variables
lynx_dat <- lynx_dat[, c("species", "decimalLongitude", "decimalLatitude")]

# rename the columns
names(lynx_dat) <- c("species", "x", "y")

# add a presence-absence column
lynx_dat$pa <- 1

# get an objecgt with the species name
sp <- lynx_dat$species[1]

# visualise the records
plot(lynx_dat$x, lynx_dat$y)

# for an interactive inspection of the occurrence points you can use mapview package
# pp <- st_as_sf(lynx_dat, coords=c('x', 'y'), crs = 4326)
# mapview::mapview(pp)

## sample background points

# set the seed
set.seed(10)

# sample based opn the environmental raster
back_dat <- sampleRandom(s_pres, size = 1000, xy = TRUE)[, 1:2] #sample random point using one of the environmental layers

# convert to a data.frame
back_dat <- as.data.frame(back_dat)
back_dat <- data.frame(species = sp, x = back_dat[, 1], y = back_dat[, 2], pa = 0)

# join presence data to background points
sp_dat <- dplyr::bind_rows(lynx_dat, back_dat)

# extract environmental variables
env_vars <- extract(s_pres, sp_dat[, c('x','y')])
head(env_vars)

# bind variables with species data
sp_dat <- dplyr::bind_cols(sp_dat, env_vars)

# get complete cases
sp_dat <- sp_dat[complete.cases(sp_dat), ] 
head(sp_dat)

## fit the model

# convert to a spatial object for thinning
sp_dat_thin <- vect(sp_dat, geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")

# randomly sample one point per raster cell
sp_dat_thin <- terra::spatSample(sp_dat_thin, 1, strata = rast(s_pres[[1]]))

# convert this back to a data.frame
sp_dat <- dplyr::bind_cols(dplyr::tibble(species = sp_dat_thin$species),
                           dplyr::as_tibble(terra::crds(sp_dat_thin)),
                           dplyr::tibble(pa = sp_dat_thin$pa))

# prepare the dataset into the biomod2 format
biomod_dat <-
  BIOMOD_FormatingData(
    resp.var = sp_dat$pa, #response variable (0, 1)
    expl.var = s_pres, #explanatory variables provided as a raster stack
    resp.name = sp_dat$species[1], #name of the species
    resp.xy = sp_dat[, c('x','y')], #coordinates of the response variable
  )

# how many runs?
n_runs <- 5

# run the models
models <-
  BIOMOD_Modeling(
  bm.format = biomod_dat,
  modeling.id = 'AllModels',
  models = c('RF', 'GLM'),
  CV.strategy = 'random',
  CV.nb.rep = n_runs,
  CV.perc = 0.8,
  OPT.strategy = 'bigboss',
  metric.eval = c('TSS','ROC'),
  var.import = n_runs,
  seed.val = 42
  )

# evaluate the models
eval <- get_evaluations(models)

# TSS > 0.5 is decent
# ROC > 0.75 is decent
print(eval)
View(eval)

# summarise these predictions
eval |>
  dplyr::filter(run != "allRun") |>
  dplyr::group_by(algo, metric.eval) |>
  dplyr::summarise(mean = mean(validation),
                   median = median(validation))

# show the predictions
preds_pres <- BIOMOD_Projection(
  bm.mod = models,
  new.env = s_pres,
  proj.name = 'current',
  selected.models = 'all',
  binary.meth = 'TSS',
  output.format = '.grd'
)

# plot the different model results
plot(preds_pres)

# extract the raster layers
preds_pres_stack <- get_predictions(preds_pres)

# plot some of these layers with the country borders
plot(preds_pres_stack[[5]])
plot(st_geometry(european_territories), add = TRUE, col = NA, border = "white")
points(sp_dat$x[sp_dat$pa == 1], sp_dat$y[sp_dat$pa == 1], col = "white", 
       pch = 19, cex = 0.25)

