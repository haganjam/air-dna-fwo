
# analyse the plant and bird diversity data

## florabank

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

# select the 20 hokken around the focal region
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

# add this to the hok_id data
site <- dplyr::tibble(hok_id = sp$hok_id,
                      alpha)
head(site)

# add some noise to these data and draw from a poisson distribution
pred_m <- vector(length = nrow(site))
pred_sd <- vector(length = nrow(site))
for (i in seq_len(nrow(site))) {
  
  x <- site$alpha[i]
  y <- x + rnorm(n = 1, mean = 0, sd = 10)
  y <- ifelse(y < 0, 0.1, y)
  z <- rpois(n = 100, lambda = y)
  
  pred_m[i] <- mean(z)
  pred_sd[i] <- sd(z)
  
}

# add this to the site data
site$pred_m <- pred_m
site$pred_sd <- pred_sd

plot(site$pred_m, site$alpha)
library(ggplot2)
ggplot(data = site,
       mapping = aes(x = alpha, y = pred_m)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_errorbar(mapping = aes(ymin = pred_m-pred_sd, ymax = pred_m+pred_sd))


# gamma-diversity
gamma <- fossil::spp.est(x = t(sp[, -1]), rand = 10, abund = FALSE, counter = FALSE, max.est ='all')

# get mean and sd
gamma_m <- gamma[nrow(gamma), ][5]
gamma_high <- gamma[nrow(gamma), ][6]
gamma_low <- gamma[nrow(gamma), ][7]

# multiple-site beta diversity
beta <- betapart::beta.sample(x = sp[, -1], sites = 15,
                              index.family = "sorensen", samples = 100)

# extract the mean and sd
beta_m <- beta$mean.values[3]
beta_sd <- beta$sd.values[3]

# pull these data into a data.frame
plant_div <- dplyr::tibble(taxon = "plant",
                           metric = c("Alpha-richness", "Gamma-richness", "Beta-diversity"),
                           mean = c(alpha_m, gamma_m, beta_m),
                           error_low = c((alpha_m - alpha_sd), (gamma_low), (beta_m - beta_sd)),
                           error_high = c((alpha_m + alpha_sd), (gamma_high), (beta_m + beta_sd)))

# arrange alphabetically
plant_div <- dplyr::arrange(plant_div, metric)

## broedvogels data

# read the cleaned broedvogels data
dat <- arrow::read_parquet(file = "02-wp-1-2/data/clean-broedvogels-belgie.parquet")
head(dat)
dim(dat)

# check for species name completeness
any(is.na(unique(dat$species)))

# get a single species entry for each hok-id
sp <-
  dat |>
  dplyr::group_by(cell_id, scientificName) |>
  dplyr::sample_n(size = 1) |>
  dplyr::ungroup() |>
  dplyr::mutate(pa = 1)

# fix the scientific names
sp$scientificName <- tolower(gsub(" ", "_", sp$scientificName))
sp$scientificName <- gsub("\\.", "", sp$scientificName)

# select the cells around the focal region
cell_targ <- c("2398", "2348", "2350", "2406", "2472", "2476",
               "2352", "2570", "2482", "2411", "2412", "2485")

# check that these are the correct sites
pp <- st_as_sf(dplyr::filter(dat, cell_id %in% cell_targ),
               coords=c('decimalLongitude', 'decimalLatitude'), crs = 4326)
# mapview::mapview(pp, legend = FALSE)

# create a site-by-species matrix
sp <-
  sp |>
  dplyr::filter(cell_id %in% cell_targ) |>
  tidyr::pivot_wider(id_cols = "cell_id",
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

# get mean and sd
gamma_m <- gamma[nrow(gamma), ][5]
gamma_high <- gamma[nrow(gamma), ][6]
gamma_low <- gamma[nrow(gamma), ][7]

# multiple-site beta diversity
beta <- betapart::beta.sample(x = sp[, -1], sites = 9,
                              index.family = "sorensen", samples = 100)

# extract the mean and sd
beta_m <- beta$mean.values[3]
beta_sd <- beta$sd.values[3]

# pull these data into a data.frame
bird_div <- dplyr::tibble(taxon = "bird",
                          metric = c("Alpha-richness", "Gamma-richness", "Beta-diversity"),
                          mean = c(alpha_m, gamma_m, beta_m),
                          error_low = c((alpha_m - alpha_sd), (gamma_low), (beta_m - beta_sd)),
                          error_high = c((alpha_m + alpha_sd), (gamma_high), (beta_m + beta_sd)))

# arrange alphabetically
bird_div <- dplyr::arrange(bird_div, metric)

## plot for proposal

# titles
titles <- c(
  expression(alpha * "-richness"),
  expression(beta * "-richness"),
  expression(gamma * "-richness")
)

# custom theme with only x-axis
theme_only_x_axis <- theme(
  strip.background = element_blank(),
  panel.background = element_blank(),         # Remove panel background
  panel.grid = element_blank(),               # Remove grid lines
  axis.title.y = element_blank(),             # Remove y-axis title
  axis.text.y = element_blank(),              # Remove y-axis text
  axis.ticks.y = element_blank(),             # Remove y-axis ticks
  axis.line.y = element_blank(),              # Remove y-axis line
  axis.title.x = element_text(),              # Keep x-axis title
  axis.text.x = element_text(),               # Keep x-axis text
  axis.ticks.x = element_line(),              # Keep x-axis ticks
  axis.line.x = element_line(),               # Keep x-axis line        # Remove plot background
  legend.position = "none"                    # Remove legend
)


# get modelled plant data
plant_list <- vector("list", length = nrow(plant_div))
for (i in seq_len(nrow(plant_div))) {
  
  # sample some random error
  e <- runif(n = 1, min = 0.05, max = 0.10)
  e <- e * sample(x = c(-1, 1), size = 1)
  
  # add this random error to the mean
  m <- plant_div[i, ]$mean
  m <- m + (m*e)
  
  # sample 100 realisations from the poisson distribution
  est <- 
    if (m >= 1) {
      rpois(n = 100, m)
    } else {
      truncnorm::rtruncnorm(n = 100, m = m, sd = 0.05, a = 0, b = 1)
    }
  
  # bind into a data.frame
  df <- dplyr::tibble(taxon = plant_div[i, ]$taxon,
                      metric = plant_div[i, ]$metric,
                      mean = est)
  
  # bind into output list
  plant_list[[i]] <- df
  
}

# make the plant plots
plant_plots <- vector("list", length = length(plant_list))
for (i in seq_along(plant_list)) {
  
  p1 <-
    ggplot() +
    geom_point(data = plant_div[i,],
               mapping = aes(x = mean, y = "B")) +
    geom_errorbarh(data = plant_div[i, ],
                   mapping = aes(xmin = error_low, xmax = error_high, y = "B"),
                   height = 0) +
    ggbeeswarm::geom_quasirandom(data = plant_list[[i]],
                                 mapping = aes(x = mean, y = "A"), width = 0.1,
                                 colour = "red", alpha = 0.25) +
    ggtitle(titles[i]) +
    xlab(NULL) +
    theme_only_x_axis +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(colour = "black"))
  
  plant_plots[[i]] <- p1
  
}

# bind into a grid
p_plot <- cowplot::plot_grid(plotlist = plant_plots, nrow = 1, ncol = 3)

# get modelled bird data
bird_list <- vector("list", length = nrow(bird_div))
for (i in seq_len(nrow(bird_div))) {
  
  # sample some random error
  e <- runif(n = 1, min = 0.05, max = 0.10)
  e <- e * sample(x = c(-1, 1), size = 1)
  
  # add this random error to the mean
  m <- bird_div[i, ]$mean
  m <- m + (m*e)
  
  # sample 100 realisations from the poisson distribution
  est <- 
    if (m >= 1) {
      rpois(n = 100, m)
    } else {
      truncnorm::rtruncnorm(n = 100, m = m, sd = 0.05, a = 0, b = 1)
    }
  
  # bind into a data.frame
  df <- dplyr::tibble(taxon = bird_div[i, ]$taxon,
                      metric = bird_div[i, ]$metric,
                      mean = est)
  
  # bind into output list
  bird_list[[i]] <- df
  
}

# make the plant plots
bird_plots <- vector("list", length = length(bird_list))
for (i in seq_along(bird_list)) {
  
  p1 <-
    ggplot() +
    geom_point(data = bird_div[i,],
               mapping = aes(x = mean, y = "B")) +
    geom_errorbarh(data = bird_div[i, ],
                   mapping = aes(xmin = error_low, xmax = error_high, y = "B"),
                   height = 0) +
    ggbeeswarm::geom_quasirandom(data = bird_list[[i]],
                                 mapping = aes(x = mean, y = "A"), width = 0.1,
                                 colour = "red", alpha = 0.25) +
    ggtitle("") +
    xlab(NULL) +
    theme_only_x_axis +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(colour = "black"))
  
  bird_plots[[i]] <- p1
  
}

# bind into a grid
b_plot <- cowplot::plot_grid(plotlist = bird_plots, nrow = 1, ncol = 3)

# combine the two plots
figx <- cowplot::plot_grid(p_plot, b_plot, nrow = 2, ncol = 1)

# export
ggsave(filename = "02-wp-1-2/figures-tables/fig_x.png", figx,
       dpi = 400, units = "cm", width = 18, height = 9)



