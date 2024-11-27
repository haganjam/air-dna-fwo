
# load the gbif download for cleaning
# https://ropensci.github.io/CoordinateCleaner/articles/Cleaning_GBIF_data_with_CoordinateCleaner.html

# load relevant libraries
library(ggplot2)
library(dplyr)
library(rgbif)

# load from the zip file
lynx_dat <- occ_download_import(as.download("03-wp-3/raw-data/occurrence-data/0014538-241107131044228.zip"))

# load CoordinateCleaner package

# install.packages("devtools")
# library(devtools)
# install_github("ropensci/CoordinateCleaner")

# load other relevant packages
library(countrycode)
library(CoordinateCleaner)
library(sf)

# extract the relevant variables

# select columns of interest
lynx_dat <- 
  lynx_dat %>%
  dplyr::select(species, decimalLongitude, 
                decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters,
                year, basisOfRecord, institutionCode)

# remove records without coordinates
lynx_dat <- 
  lynx_dat %>%
  dplyr::filter(!is.na(decimalLongitude)) %>%
  dplyr::filter(!is.na(decimalLatitude))

# visualise the coordinates

# plot data to get an overview
wm <- borders("world", colour = "gray50", fill = "gray50")
ggplot() +
  coord_fixed() +
  wm +
  geom_point(data = lynx_dat,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred",
             size = 0.5) +
  theme_bw()

# which countries are present?
sort(unique(lynx_dat$countryCode))

# extract BE, NL, DE, FR, LU
lynx_dat <-
  lynx_dat |>
  dplyr::filter(countryCode %in% c("NL", "DE", "FR", "LU", "BE"))
head(lynx_dat)

# convert country code from ISO2c to ISO3c
lynx_dat$countryCode <-  
  countrycode(lynx_dat$countryCode, origin =  'iso2c', destination = 'iso3c')

# flag problems
lynx_dat <- data.frame(lynx_dat)
flags <- clean_coordinates(x = lynx_dat, 
                           lon = "decimalLongitude", 
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("capitals", "centroids", "equal", "zeros"))

# check the flag summary
summary(flags)

# visualise the flags
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

# exclude problematic records
lynx_dat_cl <- lynx_dat[flags$.summary, ]

# check records with low precision
lynx_dat_cl %>% 
  mutate(Uncertainty = coordinateUncertaintyInMeters / 1000) %>% 
  ggplot(aes(x = Uncertainty)) + 
  geom_histogram() +
  xlab("Coordinate uncertainty in meters") +
  theme_bw()

# remove records with very low precision or no precision information
lynx_dat_cl <- 
  lynx_dat_cl %>%
  dplyr::filter(coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters))

# remove unsuitable data sources, especially fossils 
# which are responsible for the majority of problems in this case
table(lynx_dat_cl$basisOfRecord)

# only keep human observation, machine observation and preserved specimen
lynx_dat_cl <-
  lynx_dat_cl %>%
  dplyr::filter(basisOfRecord == "HUMAN_OBSERVATION" | 
                  basisOfRecord == "MACHINE_OBSERVATION" |
                  basisOfRecord == "PRESERVED_SPECIMEN")

# check individual count information
range(lynx_dat_cl$individualCount, na.rm = TRUE)

# check the record dates
range(lynx_dat_cl$year, na.rm = TRUE)
unique(lynx_dat_cl$year)

# how many are post 2000?
sum(lynx_dat_cl$year > 2000, na.rm = TRUE)

# only keep records after 2000
lynx_dat_cl <-
  lynx_dat_cl %>%
  dplyr::filter(year >= 2000, !is.na(year)) %>%
  dplyr::as_tibble()
head(lynx_dat_cl)

# check that we only have data on the Eurasian Lynx
unique(lynx_dat_cl$species)

# export this as a .csv file
readr::write_csv(x = lynx_dat_cl,
                 file = "03-wp-3/data/lynx_lynx_occ.csv")
