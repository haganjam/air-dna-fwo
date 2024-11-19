
# analyse the plant diversity data

# read the cleaned florabank1 data
dat <- arrow::read_parquet(file = "02-wp-1-2/data/clean-florabank1-2020-2024.parquet")
head(dat)
dim(dat)

# check for species name completeness
any(is.na(unique(dat$species)))

# get a single species entry for each hok-id
sp <-
  dat |>
  dplyr::group_by(hok_id, scientificName) |>
  dplyr::sample_n(size = 1) |>
  dplyr::ungroup() |>
  dplyr::mutate(pa = 1)

# fix the scientific names
sp$scientificName <- tolower(gsub(" ", "_", sp$scientificName))
sp$scientificName <- gsub("\\.", "", sp$scientificName)

# select the 23 hokken around the focal region
hok_targ <- c("2246", "2245", "2324", "2322", "2320",
              "2318", "2316", "2315", "2380", "2381", "2382",
              "2383", "2384", "2452", "2450", "2449",
              "2447", "2446", "2519", "2517")

# create a site-by-species matrix
sp <-
  sp |>
  dplyr::filter(hok_id %in% hok_targ) |>
  tidyr::pivot_wider(id_cols = "hok_id",
                     names_from = "scientificName",
                     values_from = "pa",
                     values_fill = 0)

# calculate some preliminary biodiversity metrics

# alpha-diversity
alpha <- apply(sp[, -1], 1, function(x) sum(x))

# extract mean and sd
alpha_m <- mean(alpha)
alpha_sd <- sd(alpha)

# gamma-diversity
gamma <- fossil::spp.est(x = t(sp[, -1]), rand = 10, abund = FALSE, counter = FALSE, max.est ='all')

# multiple-site beta diversity
beta <- betapart::beta.sample(x = sp[, -1], sites = 15,
                              index.family = "sorensen", samples = 100)

# extract the mean and sd
beta_m <- beta$mean.values[3]
beta_sd <- beta$sd.values[3]

# pull these data into a data.frame










