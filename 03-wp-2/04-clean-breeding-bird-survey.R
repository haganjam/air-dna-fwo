
# clean breading bird survey

# load relevant libraries
library(ggplot2)
library(dplyr)
library(rgbif)
library(sf)

# load the gbif data
bbs_dat <- occ_download_import(as.download("02-wp-1-2/raw-data/breeding-bird-survey/0021646-241107131044228.zip"))

# check the data
head(bbs_dat)
dim(bbs_dat)

# what years are present in these data?
sort(unique(bbs_dat$year))

# extract the relevant columns
bbs_dat <-
  bbs_dat |>
  dplyr::select(gbifID, family, genus, species, scientificName,
                decimalLongitude, decimalLatitude,
                coordinatePrecision, coordinateUncertaintyInMeters,
                year, month, day)

# remove records without coordinates
bbs_dat <-
  bbs_dat |>
  dplyr::filter(!is.na(decimalLongitude)) |>
  dplyr::filter(!is.na(decimalLatitude))

# remove records without year information
bbs_dat_rec <-
  bbs_dat |>
  dplyr::filter(!is.na(year))
head(bbs_dat_rec)
dim(bbs_dat_rec)

# how many unique longitude values are there? there are clearly centroids
length(unique(bbs_dat_rec$decimalLatitude))
length(unique(bbs_dat_rec$decimalLongitude))

# combine the latitude and longitude into a 1 x 1 km identifier
x <- paste(bbs_dat_rec$decimalLatitude, bbs_dat_rec$decimalLongitude, sep = "_")
y <- as.integer(factor(x))

# cell_id
bbs_dat_rec$cell_id <- as.character(y)

# sample one coordinate randomly per cell
cell_samp <-
  bbs_dat_rec |>
  dplyr::select(cell_id, decimalLongitude, decimalLatitude) |>
  dplyr::group_by(cell_id) |>
  dplyr::sample_n(size = 1) |>
  dplyr::ungroup()
head(cell_samp)

# plot the samples
plot(cell_samp$decimalLongitude, cell_samp$decimalLatitude)

# for an interactive inspection of the occurrence points you can use mapview package
pp <- st_as_sf(cell_samp, coords=c('decimalLongitude', 'decimalLatitude'), crs = 4326)
mapview::mapview(pp, legend = FALSE)

# export a cleaned version of the data
arrow::write_parquet(x = bbs_dat_rec, "02-wp-1-2/data/clean-broedvogels-belgie.parquet")
