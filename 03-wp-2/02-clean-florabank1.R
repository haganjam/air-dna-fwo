
# clean the florabank1 data

# load relevant libraries
library(ggplot2)
library(dplyr)
library(rgbif)
library(sf)

# load the gbif data
f1_dat <- occ_download_import(as.download("02-wp-1-2/raw-data/florabank1/0014629-241107131044228.zip"))

# check the data
head(f1_dat)
dim(f1_dat)

# what years are present in these data?
sort(unique(f1_dat$year))

# extract the relevant columns
f1_dat <-
  f1_dat |>
  dplyr::select(gbifID, family, genus, species, scientificName,
                decimalLongitude, decimalLatitude,
                coordinatePrecision, coordinateUncertaintyInMeters,
                year, month, day)

# remove records without coordinates
f1_dat <-
  f1_dat |>
  dplyr::filter(!is.na(decimalLongitude)) |>
  dplyr::filter(!is.na(decimalLatitude))

# remove records without year information
f1_dat <-
  f1_dat |>
  dplyr::filter(!is.na(year))

# extract records from 2020 onwards
f1_dat_rec <-
  f1_dat |>
  dplyr::filter(year %in% c(2020, 2021, 2022, 2023, 2024))
head(f1_dat_rec)
dim(f1_dat_rec)

# how many unique longitude values are there? there are clearly centroids
length(unique(f1_dat_rec$decimalLatitude))
length(unique(f1_dat_rec$decimalLongitude))

# combine the latitude and longitude into a 1 x 1 km identifier
x <- paste(f1_dat_rec$decimalLatitude, f1_dat_rec$decimalLongitude, sep = "_")
y <- as.integer(factor(x))

# hok_id
f1_dat_rec$hok_id <- as.character(y)

# sample one coordinate randomly per hok
hok_samp <-
  f1_dat_rec |>
  dplyr::select(hok_id, decimalLongitude, decimalLatitude) |>
  dplyr::group_by(hok_id) |>
  dplyr::sample_n(size = 1) |>
  dplyr::ungroup()
head(hok_samp)

# plot the samples
plot(hok_samp$decimalLongitude, hok_samp$decimalLatitude)

# for an interactive inspection of the occurrence points you can use mapview package
pp <- st_as_sf(hok_samp, coords=c('decimalLongitude', 'decimalLatitude'), crs = 4326)
mapview::mapview(pp, legend = FALSE)

# select the sites to use in Google Earth Engine
sites <- c(1198, 1123, 1122, 1194, 1193, 1267, 1268, 1348, 1350, 1436, 1434, 
           1527, 1529, 1639, 1729, 1812, 1811, 1813, 1878, 1815)

# get a table of these sites
hok_samp |>
  dplyr::filter(hok_id %in% sites) |>
  dplyr::mutate(decimalLongitude = round(decimalLongitude, 7),
                decimalLatitude = round(decimalLatitude, 7)) |>
  as.data.frame()


# export a cleaned version of the data
arrow::write_parquet(x = f1_dat_rec, "02-wp-1-2/data/clean-florabank1-2020-2024.parquet")
