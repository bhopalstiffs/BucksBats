# =============================================================================
# Bucks Bats: reproducible analysis pipeline
# =============================================================================
# This script:
#   (1) rebuilds or loads effort-corrected analytical data
#   (2) validates the analytical sampling frame
#   (3) optionally runs short model-development diagnostics
#   (4) fits or loads final Bayesian GAMs
#   (5) evaluates final-model fit and residual temporal structure
#   (6) records software versions and package citations
# =============================================================================

source("scripts/summarize_functions.R")
source("scripts/diagnostic_functions.R")
source("scripts/bibliography_functions.R")

# 0. Setup --------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(brms)
library(cmdstanr)
library(tidybayes)
library(patchwork)
library(RColorBrewer)
library(flextable)
library(officer)
library(posterior)

options(mc.cores = parallel::detectCores())

# Directory layout (relative to analysis/) -----------------------------------
data_dir       <- "data/"
model_fits_dir <- "model_fits/"
table_dir      <- "tables/"
script_dir     <- "scripts/"
derived_dir    <- "derived/"
manuscript_dir <- "manuscript/"
reference_dir  <- file.path(manuscript_dir, "references")

walk(
  c(
    data_dir,
    model_fits_dir,
    table_dir,
    derived_dir,
    reference_dir
  ),
  ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE)
)

# -----------------------------------------------------------------------------
# Workflow toggles
# -----------------------------------------------------------------------------

# TRUE  = rebuild formatted data from the original raw acoustic/environment files
# FALSE = load the previously formatted data from data/formatted_data.RData
format_raw_data <- FALSE

# Short one-chain model-development fits
run_aggregate_diagnostic <- FALSE
run_frequency_diagnostics <- FALSE

# Final production fits
fit_aggregate_models <- FALSE
fit_frequency_models <- FALSE

update_package_bibliography <- FALSE

# =============================================================================
# 1-5. Data preparation
# =============================================================================

formatted_file <- file.path(data_dir, "formatted_data.RData")

if (format_raw_data) {
  
  message("Rebuilding formatted data from raw files...")
  
# ---------------------------------------------------------------------------
# 1. Load raw data
# ---------------------------------------------------------------------------
  
bat_dir <- "D:/Stillwater/BucksBats/classified data/BatData2021"
env_path <- "D:/Stillwater/BucksBats/CovariateData/envirodata.xlsx"

data <- list.files(bat_dir, pattern = "*.csv", full.names = TRUE) |>
  map_df(~ read_csv(.x))

envirodatafull <- env_path |>
  excel_sheets() |>
  set_names() |>
  map_df(~ read_excel(env_path, sheet = .x))


# 2. Format bat acoustic data -------------------------------------------------

cdata_species <- data |>
  select(Site, DATE:`HOUR-12`, `AUTO ID*`) |>
  rename(Species = `AUTO ID*`) |>
  filter(!Species %in% c("NoID", "Noise")) |>
  mutate(
    Time     = str_sub(TIME, 3),
    Time1    = paste0(HOUR, Time),
    DateTime = as.POSIXct(paste(DATE, Time1),
                          format = "%m/%d/%Y %H:%M:%S"),
    NiteDay  = as_date(if_else(hour(DateTime) < 12,
                               DateTime - days(1), DateTime)),
    tc       = cut(DateTime, breaks = "5 min")
  ) |>
  group_by(Site, NiteDay, tc, Species) |> tally() |>
  group_by(Site, NiteDay, Species)      |> tally() |>
  ungroup() |>
  mutate(datekey = strftime(NiteDay, format = "%Y-%m-%d"))

# Species --> characteristic frequency group
species_freq <- tribble(
  ~Species, ~freq_group, ~char_freq_kHz, ~common_name,
  # HIGH (>40 kHz) — most strongly affected by humidity attenuation
  "MYOYUM", "high", 48, "Yuma myotis",
  "MYOLUC", "high", 42, "little brown",
  "MYOCAL", "high", 50, "California myotis",
  "MYOCIL", "high", 45, "western small-footed",
  "MYOVOL", "high", 40, "long-legged myotis",
  "PARHES", "high", 45, "canyon bat",
  "LASBLO", "high", 42, "western red",
  # MID (25-40 kHz)
  "MYOEVO", "mid",  35, "long-eared myotis",
  "MYOTHY", "mid",  32, "fringed myotis",
  "EPTFUS", "mid",  27, "big brown",
  "LASNOC", "mid",  26, "silver-haired",
  "TADBRA", "mid",  25, "Mexican free-tailed",
  "ANTPAL", "mid",  30, "pallid bat",
  "CORTOW", "mid",  28, "Townsend's big-eared",
  # LOW (<25 kHz)
  "LASCIN", "low",  20, "hoary",
  "EUMPER", "low",  12, "western mastiff",
  "EUDMAC", "low",  11, "spotted bat"
)


# 3. Format environmental data ------------------------------------------------

envirodata <- envirodatafull |>
  rename(avetemp = `Average of AIRC`, averh = `Average of RH`) |>
  select(-c(`Max of AIRC`, `Min of AIRC`, `Max of RH`, `Min of RH`)) |>
  mutate(across(c(avetemp, averh),
                .names = "scale_{.col}",
                ~ scale(.) |> as.vector())) |>
  mutate(
    habitat1 = case_when(
      Site == "BUC3"  ~ "river",
      Site == "BUL11" ~ "reservoir",
      Site == "LBL5"  ~ "reservoir",
      Site == "GRP5"  ~ "forest",
      Site == "GRP1"  ~ "river",
      Site == "BUL18" ~ "reservoir",
      Site == "THL1"  ~ "reservoir",
      Site == "BUL4"  ~ "reservoir"
    ),
    habitat2 = case_when(
      Site == "BUC3"  ~ "forest",
      Site == "BUL11" ~ "forest/reservoir",
      Site == "LBL5"  ~ "reservoir/forest",
      Site == "GRP5"  ~ "forest/shrub",
      Site == "GRP1"  ~ "forest",
      Site == "BUL18" ~ "forest/shrub",
      Site == "THL1"  ~ "montane forest",
      Site == "BUL4"  ~ "forest/shrub"
    ),
    elevation = case_when(
      Site == "BUC3"  ~ "low",
      Site == "BUL11" ~ "middle",
      Site == "LBL5"  ~ "middle",
      Site == "GRP5"  ~ "middle",
      Site == "GRP1"  ~ "middle",
      Site == "BUL18" ~ "high",
      Site == "THL1"  ~ "high",
      Site == "BUL4"  ~ "middle"
    ),
    datekey = strftime(Date, format = "%Y-%m-%d")
  )

# 4. Join and aggregate -------------------------------------------------------

species_list <- species_freq$Species

# Identify site-nights with evidence that the detector was operating
detector_effort <- data |>
  mutate(
    Time = str_sub(TIME, 3),
    Time1 = paste0(HOUR, Time),
    DateTime = as.POSIXct(
      paste(DATE, Time1),
      format = "%m/%d/%Y %H:%M:%S"
    ),
    NiteDay = as_date(
      if_else(
        hour(DateTime) < 12,
        DateTime - days(1),
        DateTime
      )
    ),
    datekey = strftime(NiteDay, format = "%Y-%m-%d")
  ) |>
  group_by(Site, datekey) |>
  summarise(
    n_raw_files = n(),
    .groups = "drop"
  )

# Deployment skeleton contains only nights with evidence
# that the acoustic detector was operating
deployment_skeleton <- envirodata |>
  select(Site, datekey) |>
  distinct() |>
  inner_join(
    detector_effort |> select(Site, datekey),
    by = c("Site", "datekey")
  ) |>
  crossing(Species = species_list)

# Fill nondetections with true zeros only on sampled nights
species_activity <- deployment_skeleton |>
  left_join(
    cdata_species |> select(Site, datekey, Species, n),
    by = c("Site", "datekey", "Species")
  ) |>
  mutate(n = replace_na(n, 0))

# Add environmental covariates and frequency-group lookup
alldata_species <- species_activity |>
  left_join(envirodata,   by = c("Site", "datekey")) |>
  left_join(species_freq, by = "Species") |>
  mutate(across(c(Site, habitat1, habitat2, elevation, freq_group),
                as.factor),
         NiteDayDate = as_date(datekey),
         doy         = yday(NiteDayDate)) |>
  na.omit()

# Sum activity across species within each frequency group, per site-night
alldata_freq_group <- alldata_species |>
  group_by(Site, datekey, freq_group, NiteDayDate, doy,
           avetemp, averh, scale_avetemp, scale_averh,
           Lunar, habitat1, habitat2, elevation) |>
  summarise(n = sum(n), .groups = "drop") |>
  mutate(freq_group = factor(freq_group, levels = c("low", "mid", "high")))

# Community-total activity, per site-night
agg_data <- alldata_freq_group |>
  group_by(Site, datekey, NiteDayDate, doy,
           avetemp, averh, scale_avetemp, scale_averh,
           Lunar, habitat1, habitat2, elevation) |>
  summarise(n = sum(n), .groups = "drop")
  
  
# ---------------------------------------------------------------------------
# 5. Save formatted data
# ---------------------------------------------------------------------------
  
  save(
    cdata_species,
    species_freq,
    envirodata,
    alldata_species,
    alldata_freq_group,
    agg_data,
    file = formatted_file
  )
  
  message("Formatted data saved to: ", formatted_file)
  
} else {
  
# ---------------------------------------------------------------------------
# Load previously formatted data
# ---------------------------------------------------------------------------
  
  if (!file.exists(formatted_file)) {
    stop(
      "Formatted data file not found: ", formatted_file,
      "\nEither restore the formatted data file or set format_raw_data <- TRUE ",
      "and provide access to the raw data."
    )
  }
  
  load(formatted_file)
  
  message("Loaded formatted data from: ", formatted_file)
}

# -----------------------------------------------------------------------------
# Validate formatted analytical data
# -----------------------------------------------------------------------------

expected_site_nights <- 1056L
expected_total_activity <- 103681L
expected_n_sites <- 8L
expected_n_species <- 17L

if (nrow(agg_data) != expected_site_nights) {
  stop(
    "agg_data contains ", nrow(agg_data),
    " site-nights; expected ", expected_site_nights,
    ". Check detector-effort filtering."
  )
}

if (n_distinct(agg_data$Site) != expected_n_sites) {
  stop("Unexpected number of sites in agg_data.")
}

if (sum(agg_data$n) != expected_total_activity) {
  stop("Unexpected total activity in agg_data.")
}

if (nrow(alldata_freq_group) != expected_site_nights * 3L) {
  stop("Unexpected number of rows in alldata_freq_group.")
}

if (n_distinct(alldata_species$Species) != expected_n_species) {
  stop("Unexpected number of species in alldata_species.")
}

message(
  "Data integrity check passed: ",
  expected_site_nights,
  " sampled site-nights; ",
  expected_total_activity,
  " activity intervals."
)


# 6. Fit (or load) models -----------------------------------------------------

# Final model-development decisions:
#   * Negative-binomial likelihood rather than Poisson to accommodate
#     extra-Poisson variation and improve HMC geometry.
#   * Site-specific factor smooths of day of year to represent differences in
#     seasonal phenology among sites.
#   * No AR(1) term: calendar-aware residual diagnostics showed that the long
#     positive autocorrelation tail largely disappeared after site-specific
#     seasonality was represented. Moderate short-lag dependence remains.
#
# IMPORTANT: the aggregate conditional model has already been fitted and should
# be saved as m_agg_env_final.rds rather than unnecessarily refitted.

final_prior <- prior(normal(0, 10), class = Intercept)

final_control <- list(
  adapt_delta   = 0.995,
  max_treedepth = 12
)

final_args <- list(
  family  = negbinomial(),
  prior   = final_prior,
  chains  = 4,
  iter    = 4000,
  warmup  = 1000,
  cores   = 4,
  backend = "cmdstanr",
  control = final_control
)

# Primary conditional environmental formula.
# The factor smooth is fully penalized, so a separate (1 | Site) term is not
# included.
env_formula <- bf(
  n ~
    s(doy) +
    s(doy, Site, bs = "fs") +
    s(scale_avetemp) +
    s(scale_averh) +
    s(Lunar)
)

# Total seasonal formula. This model is intentionally not conditioned on
# temperature or RH because those variables lie on seasonal pathways of
# interest. Fit only after the short diagnostic fit below behaves acceptably.
season_formula <- bf(
  n ~
    s(doy) +
    s(doy, Site, bs = "fs")
)

# -----------------------------------------------------------------------------
# 6a. Optional short diagnostic fits
# -----------------------------------------------------------------------------
# One-chain fits used only for model development and diagnostics.
# Biological inference is based on the final four-chain production fits.

diagnostic_args <- list(
  family  = negbinomial(),
  prior   = final_prior,
  chains  = 1,
  iter    = 1500,
  warmup  = 750,
  cores   = 1,
  backend = "cmdstanr",
  control = final_control
)

fit_diagnostic_env <- function(dat) {
  do.call(
    brm,
    c(
      list(
        formula = env_formula,
        data    = dat
      ),
      diagnostic_args
    )
  )
}


# Aggregate seasonal diagnostic -----------------------------------------------

if (run_aggregate_diagnostic) {
  
  m_agg_season_test <- do.call(
    brm,
    c(
      list(
        formula = season_formula,
        data    = agg_data
      ),
      diagnostic_args
    )
  )
}


# Frequency-group datasets ----------------------------------------------------

freq_data <- alldata_freq_group |>
  mutate(
    Site = factor(Site),
    freq_group = factor(
      freq_group,
      levels = c("low", "mid", "high")
    )
  )

low_data <- freq_data |>
  filter(freq_group == "low")

mid_data <- freq_data |>
  filter(freq_group == "mid")

high_data <- freq_data |>
  filter(freq_group == "high")


# Frequency-group diagnostics -------------------------------------------------

if (run_frequency_diagnostics) {
  
  m_low_test <- fit_diagnostic_env(low_data)
  
  m_mid_test <- update(
    m_low_test,
    newdata = mid_data,
    recompile = FALSE
  )
  
  m_high_test <- update(
    m_low_test,
    newdata = high_data,
    recompile = FALSE
  )
}

# -----------------------------------------------------------------------------
# 6b. Production aggregate models
# -----------------------------------------------------------------------------

if (fit_aggregate_models) {

  m_agg_season <- do.call(
    brm,
    c(
      list(
        formula = season_formula,
        data    = agg_data
      ),
      final_args
    )
  )

  m_agg_env <- do.call(
    brm,
    c(
      list(
        formula = env_formula,
        data    = agg_data
      ),
      final_args
    )
  )

  saveRDS(m_agg_season,
          file.path(model_fits_dir, "m_agg_season_final.rds"))
  saveRDS(m_agg_env,
          file.path(model_fits_dir, "m_agg_env_final.rds"))

} else {

  # Load whichever final aggregate models currently exist. The conditional
  # environmental model should exist after saving the validated fit from the
  # model-development session. The seasonal model can be added after its
  # diagnostic fit has been checked.

  env_file <- file.path(model_fits_dir, "m_agg_env_final.rds")
  season_file <- file.path(model_fits_dir, "m_agg_season_final.rds")

  if (!file.exists(env_file)) {
    stop(
      "Final aggregate environmental model not found: ", env_file,
      "\nSave the validated model with saveRDS(m_env_conditional_nb_fs, '",
      env_file, "')."
    )
  }

  m_agg_env <- readRDS(env_file)

  m_agg_season <- if (file.exists(season_file)) {
    readRDS(season_file)
  } else {
    NULL
  }
}

# -----------------------------------------------------------------------------
# 6c. Production frequency-group models
# -----------------------------------------------------------------------------
# Frequency-group models are sensitivity analyses used to evaluate whether
# environmental responses in the aggregate activity model are broadly
# consistent across acoustic frequency groups.
#
# Model structure was retained from the final aggregate environmental model:
# negative-binomial likelihood, population-level seasonal smooth, and
# site-specific seasonal factor smooths.
#
# Short one-chain diagnostic fits indicated acceptable model behavior for all
# three frequency groups. Remaining posterior predictive discrepancies are
# treated as group-specific distributional heterogeneity rather than grounds
# for further model expansion.

freq_model_files <- c(
  low  = file.path(model_fits_dir, "m_low_final.rds"),
  mid  = file.path(model_fits_dir, "m_mid_final.rds"),
  high = file.path(model_fits_dir, "m_high_final.rds")
)

if (fit_frequency_models) {
  
  message("Fitting final low-frequency model...")
  
  m_low <- do.call(
    brm,
    c(
      list(
        formula = env_formula,
        data    = low_data
      ),
      final_args
    )
  )
  
  saveRDS(
    m_low,
    freq_model_files["low"]
  )
  
  
  message("Fitting final mid-frequency model...")
  
  m_mid <- do.call(
    brm,
    c(
      list(
        formula = env_formula,
        data    = mid_data
      ),
      final_args
    )
  )
  
  saveRDS(
    m_mid,
    freq_model_files["mid"]
  )
  
  
  message("Fitting final high-frequency model...")
  
  m_high <- do.call(
    brm,
    c(
      list(
        formula = env_formula,
        data    = high_data
      ),
      final_args
    )
  )
  
  saveRDS(
    m_high,
    freq_model_files["high"]
  )
  
  message("Final frequency-group models fitted and saved.")
  
} else {
  
  missing_freq_models <- freq_model_files[
    !file.exists(freq_model_files)
  ]
  
  if (length(missing_freq_models) > 0) {
    
    message(
      "Final frequency-group models not yet available: ",
      paste(names(missing_freq_models), collapse = ", "),
      ". Set fit_frequency_models <- TRUE to fit them."
    )
    
    m_low  <- NULL
    m_mid  <- NULL
    m_high <- NULL
    
  } else {
    
    m_low  <- readRDS(freq_model_files["low"])
    m_mid  <- readRDS(freq_model_files["mid"])
    m_high <- readRDS(freq_model_files["high"])
    
    message("Loaded final frequency-group models.")
  }
}

# =============================================================================
# 8. Final-model diagnostics
# =============================================================================

# -----------------------------------------------------------------------------
# Aggregate seasonal model
# -----------------------------------------------------------------------------

acf_agg_season <- calendar_acf(
  m_agg_season,
  agg_data
)

write_csv(
  acf_agg_season,
  file.path(table_dir, "aggregate_seasonal_residual_acf.csv")
)

pp_season_density <- pp_check(
  m_agg_season,
  type = "dens_overlay",
  ndraws = 100
)

pp_season_mean <- pp_check(
  m_agg_season,
  type = "stat",
  stat = "mean",
  ndraws = 100
)

pp_season_sd <- pp_check(
  m_agg_season,
  type = "stat",
  stat = "sd",
  ndraws = 100
)

pp_season_zero <- pp_check(
  m_agg_season,
  type = "stat",
  stat = function(x) mean(x == 0),
  ndraws = 100
)


# -----------------------------------------------------------------------------
# Aggregate conditional environmental model
# -----------------------------------------------------------------------------

acf_agg_env <- calendar_acf(
  m_agg_env,
  agg_data
)

write_csv(
  acf_agg_env,
  file.path(table_dir, "aggregate_environmental_residual_acf.csv")
)

pp_env_density <- pp_check(
  m_agg_env,
  type = "dens_overlay",
  ndraws = 100
)

pp_env_mean <- pp_check(
  m_agg_env,
  type = "stat",
  stat = "mean",
  ndraws = 100
)

pp_env_sd <- pp_check(
  m_agg_env,
  type = "stat",
  stat = "sd",
  ndraws = 100
)

pp_env_zero <- pp_check(
  m_agg_env,
  type = "stat",
  stat = function(x) mean(x == 0),
  ndraws = 100
)


# 10. Output generation -------------------------------------------------------

write_csv(
  acf_low_test,
  file.path(table_dir, "acf_low_test.csv")
)

write_csv(
  acf_mid_test,
  file.path(table_dir, "acf_mid_test.csv")
)

write_csv(
  acf_high_test,
  file.path(table_dir, "acf_high_test.csv")
)
# =============================================================================
# 9. Software citations and versions
# =============================================================================

manuscript_packages <- c(
  "tidyverse",
  "brms",
  "cmdstanr",
  "posterior",
  "tidybayes",
  "mgcv",
  "bayesplot",
  "ggplot2"
)

if (update_package_bibliography) {
  
  software_versions <- write_package_bib(
    packages = manuscript_packages,
    bib_file = file.path(
      reference_dir,
      "packages.bib"
    ),
    version_file = file.path(
      table_dir,
      "software_versions.csv"
    )
  )
}
