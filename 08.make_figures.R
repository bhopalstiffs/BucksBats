# =============================================================================
# 08.make_figures.R
# Bucks Bats: generate manuscript figures from final model outputs
# =============================================================================
#
# This script is downstream of model fitting.
#
# Inputs:
#   data/formatted_data.RData
#   model_fits/m_agg_season_final.rds
#   model_fits/m_agg_env_final.rds
#   model_fits/m_low_final.rds
#   model_fits/m_mid_final.rds
#   model_fits/m_high_final.rds
#
# Derived posterior predictions:
#   derived/posterior_predictions/
#
# Figure output:
#   figs/
#
# Set rebuild_posterior_predictions <- TRUE only when the fitted models,
# prediction definitions, or prediction grids have changed. Cosmetic figure
# changes should use the saved posterior draws and do not require recomputing
# predictions.
# =============================================================================


# 0. Setup --------------------------------------------------------------------

library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)
library(RColorBrewer)
library(lubridate)

source("scripts/summarize_functions.R")
source("scripts/figure_functions.R")

data_file     <- file.path("data", "formatted_data.RData")
model_dir     <- "model_fits"
fig_dir       <- "figs"
supp_fig_dir  <- file.path("figs", "supplementary")
posterior_dir <- file.path("derived", "posterior_predictions")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(posterior_dir, recursive = TRUE, showWarnings = FALSE)

rebuild_posterior_predictions <- TRUE

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path, call. = FALSE)
  }
  invisible(path)
}


# 1. Load formatted analytical data -------------------------------------------

require_file(data_file)
load(data_file)

freq_data <- alldata_freq_group |>
  mutate(
    Site = factor(Site),
    freq_group = factor(freq_group, levels = c("low", "mid", "high"))
  )

low_data  <- freq_data |> filter(freq_group == "low")
mid_data  <- freq_data |> filter(freq_group == "mid")
high_data <- freq_data |> filter(freq_group == "high")


# 2. Load final fitted models --------------------------------------------------

model_files <- c(
  agg_season = file.path(model_dir, "m_agg_season_final.rds"),
  agg_env    = file.path(model_dir, "m_agg_env_final.rds"),
  low        = file.path(model_dir, "m_low_final.rds"),
  mid        = file.path(model_dir, "m_mid_final.rds"),
  high       = file.path(model_dir, "m_high_final.rds")
)

walk(model_files, require_file)

m_agg_season <- readRDS(model_files["agg_season"])
m_agg_env    <- readRDS(model_files["agg_env"])
m_low        <- readRDS(model_files["low"])
m_mid        <- readRDS(model_files["mid"])
m_high       <- readRDS(model_files["high"])


# 3. Posterior prediction files -----------------------------------------------

posterior_files <- c(
  season_total       = file.path(posterior_dir, "season_total_draws.rds"),
  season_conditional = file.path(posterior_dir, "season_conditional_draws.rds"),
  season_sites       = file.path(posterior_dir, "season_site_draws.rds"),
  temp_low           = file.path(posterior_dir, "temp_low_draws.rds"),
  temp_mid           = file.path(posterior_dir, "temp_mid_draws.rds"),
  temp_high          = file.path(posterior_dir, "temp_high_draws.rds"),
  rh_low             = file.path(posterior_dir, "rh_low_draws.rds"),
  rh_mid             = file.path(posterior_dir, "rh_mid_draws.rds"),
  rh_high            = file.path(posterior_dir, "rh_high_draws.rds"),
  lunar_low          = file.path(posterior_dir, "lunar_low_draws.rds"),
  lunar_mid          = file.path(posterior_dir, "lunar_mid_draws.rds"),
  lunar_high         = file.path(posterior_dir, "lunar_high_draws.rds")
)


# 4. Rebuild posterior predictions when requested -----------------------------

if (rebuild_posterior_predictions) {

  message("Rebuilding posterior prediction files...")

  doy_seq <- seq(
    min(agg_data$doy, na.rm = TRUE),
    max(agg_data$doy, na.rm = TRUE),
    by = 1
  )

  post_season_total <- get_predictions(
    model = m_agg_season, data = agg_data,
    focal = "doy", focal_seq = doy_seq,
    prediction = "site_average"
  )

  post_season_conditional <- get_predictions(
    model = m_agg_env, data = agg_data,
    focal = "doy", focal_seq = doy_seq,
    prediction = "site_average"
  )

  post_season_sites <- get_site_season_predictions(
    model = m_agg_season, data = agg_data, doy_seq = doy_seq
  )

  saveRDS(post_season_total,       posterior_files["season_total"])
  saveRDS(post_season_conditional, posterior_files["season_conditional"])
  saveRDS(post_season_sites,       posterior_files["season_sites"])

  temp_seq <- seq(
    min(agg_data$avetemp, na.rm = TRUE),
    max(agg_data$avetemp, na.rm = TRUE),
    length.out = 100
  )

  rh_seq <- seq(
    min(agg_data$averh, na.rm = TRUE),
    max(agg_data$averh, na.rm = TRUE),
    length.out = 100
  )

  lunar_seq <- seq(
    min(agg_data$Lunar, na.rm = TRUE),
    max(agg_data$Lunar, na.rm = TRUE),
    length.out = 100
  )

  rh_scaled_seq <- (
    rh_seq - mean(agg_data$averh, na.rm = TRUE)
  ) / sd(agg_data$averh, na.rm = TRUE)

  post_temp_low <- get_predictions_chunked(
    model = m_low, data = low_data,
    focal = "scale_avetemp", focal_seq = temp_seq,
    scaled_var = "scale_avetemp", raw_var = "avetemp",
    prediction = "site_average", season = "average", chunk_size = 5
  )

  post_temp_mid <- get_predictions_chunked(
    model = m_mid, data = mid_data,
    focal = "scale_avetemp", focal_seq = temp_seq,
    scaled_var = "scale_avetemp", raw_var = "avetemp",
    prediction = "site_average", season = "average", chunk_size = 5
  )

  post_temp_high <- get_predictions_chunked(
    model = m_high, data = high_data,
    focal = "scale_avetemp", focal_seq = temp_seq,
    scaled_var = "scale_avetemp", raw_var = "avetemp",
    prediction = "site_average", season = "average", chunk_size = 5
  )

  saveRDS(post_temp_low,  posterior_files["temp_low"])
  saveRDS(post_temp_mid,  posterior_files["temp_mid"])
  saveRDS(post_temp_high, posterior_files["temp_high"])

  post_rh_low <- get_predictions_chunked(
    model = m_low, data = low_data,
    focal = "scale_averh", focal_seq = rh_scaled_seq,
    prediction = "site_average", season = "average"
  ) |>
    mutate(
      averh = scale_averh * sd(agg_data$averh, na.rm = TRUE) +
        mean(agg_data$averh, na.rm = TRUE)
    )

  post_rh_mid <- get_predictions_chunked(
    model = m_mid, data = mid_data,
    focal = "scale_averh", focal_seq = rh_scaled_seq,
    prediction = "site_average", season = "average"
  ) |>
    mutate(
      averh = scale_averh * sd(agg_data$averh, na.rm = TRUE) +
        mean(agg_data$averh, na.rm = TRUE)
    )

  post_rh_high <- get_predictions_chunked(
    model = m_high, data = high_data,
    focal = "scale_averh", focal_seq = rh_scaled_seq,
    prediction = "site_average", season = "average"
  ) |>
    mutate(
      averh = scale_averh * sd(agg_data$averh, na.rm = TRUE) +
        mean(agg_data$averh, na.rm = TRUE)
    )

  saveRDS(post_rh_low,  posterior_files["rh_low"])
  saveRDS(post_rh_mid,  posterior_files["rh_mid"])
  saveRDS(post_rh_high, posterior_files["rh_high"])

  post_lunar_low <- get_predictions_chunked(
    model = m_low, data = low_data,
    focal = "Lunar", focal_seq = lunar_seq,
    prediction = "site_average", season = "average"
  )

  post_lunar_mid <- get_predictions_chunked(
    model = m_mid, data = mid_data,
    focal = "Lunar", focal_seq = lunar_seq,
    prediction = "site_average", season = "average"
  )

  post_lunar_high <- get_predictions_chunked(
    model = m_high, data = high_data,
    focal = "Lunar", focal_seq = lunar_seq,
    prediction = "site_average", season = "average"
  )

  saveRDS(post_lunar_low,  posterior_files["lunar_low"])
  saveRDS(post_lunar_mid,  posterior_files["lunar_mid"])
  saveRDS(post_lunar_high, posterior_files["lunar_high"])

  message("Posterior prediction files rebuilt.")

} else {

  # Cosmetic figure changes should normally use this branch.
  walk(posterior_files, require_file)

  post_season_total       <- readRDS(posterior_files["season_total"])
  post_season_conditional <- readRDS(posterior_files["season_conditional"])
  post_season_sites       <- readRDS(posterior_files["season_sites"])

  post_temp_low  <- readRDS(posterior_files["temp_low"])
  post_temp_mid  <- readRDS(posterior_files["temp_mid"])
  post_temp_high <- readRDS(posterior_files["temp_high"])

  post_rh_low  <- readRDS(posterior_files["rh_low"])
  post_rh_mid  <- readRDS(posterior_files["rh_mid"])
  post_rh_high <- readRDS(posterior_files["rh_high"])

  post_lunar_low  <- readRDS(posterior_files["lunar_low"])
  post_lunar_mid  <- readRDS(posterior_files["lunar_mid"])
  post_lunar_high <- readRDS(posterior_files["lunar_high"])

  message("Loaded saved posterior prediction files.")
}


# 5. Main-text figures ---------------------------------------------------------

fig_season <- make_figure3(
  post_total       = post_season_total,
  post_conditional = post_season_conditional,
  post_sites       = post_season_sites,
  fig_dir          = fig_dir,
  filename         = "Figure3_seasonality.PNG"
)

fig_rh <- make_rh_figure(
  post_low  = post_rh_low,
  post_mid  = post_rh_mid,
  post_high = post_rh_high,
  fig_dir   = fig_dir,
  filename  = "Figure4_humidity.PNG"
)

fig_temp <- make_temperature_figure(
  post_low  = post_temp_low,
  post_mid  = post_temp_mid,
  post_high = post_temp_high,
  fig_dir   = fig_dir,
  filename  = "Figure5_temperature.PNG"
)

fig_lunar <- make_lunar_figure(
  post_low  = post_lunar_low,
  post_mid  = post_lunar_mid,
  post_high = post_lunar_high,
  fig_dir   = fig_dir,
  filename  = "Figure6_lunar.PNG"
)

if (interactive()) {
  print(fig_season)
  print(fig_temp)
  print(fig_rh)
  print(fig_lunar)
}

message("Main-text figures written to: ", fig_dir)


# 6. Supplementary figures -----------------------------------------------------

fig_s1 <- make_figure_s1(
  post_sites   = post_season_sites,
  post_average = post_season_total,
  fig_dir      = supp_fig_dir
)

fig_s2 <- make_figure_s2(
  models = list(
    Aggregate = m_agg_env,
    Low       = m_low,
    Mid       = m_mid,
    High      = m_high
  ),
  data_list = list(
    Aggregate = agg_data,
    Low       = low_data,
    Mid       = mid_data,
    High      = high_data
  ),
  fig_dir = supp_fig_dir
)

fig_s3 <- make_figure_s3(
  model   = m_agg_env,
  fig_dir = supp_fig_dir
)

fig_s4 <- make_figure_s4(
  model   = m_agg_env,
  data    = agg_data,
  fig_dir = supp_fig_dir
)

fig_s5 <- make_figure_s5(
  models = list(
    Low  = m_low,
    Mid  = m_mid,
    High = m_high
  ),
  fig_dir = supp_fig_dir
)

if (interactive()) {
  print(fig_s1)
  print(fig_s2)
  print(fig_s3)
  print(fig_s4)
  print(fig_s5)
}

message("Supplementary figures written to: ", supp_fig_dir)
