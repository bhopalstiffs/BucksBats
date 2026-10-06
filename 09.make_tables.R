# =============================================================================
# 09.make_tables.R
#
# Generate main-text and supplementary tables for the BucksBats manuscript.
#
# Main tables:
#   Table 1 - Study design and environmental conditions
#   Table 2 - Environmental contrasts by frequency group
#
# Supplementary tables:
#   Table S1 - Final model specifications
#   Table S2 - Aggregate-model sampling diagnostics
#   Table S3 - Frequency-group sampling diagnostics
#   Table S4 - Species composition and frequency classification
#   Table S5 - Retained sampling effort by site
#   Table S6 - Software and computational environment
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------

library(tidyverse)
library(flextable)
library(officer)
library(posterior)
library(brms)


# -----------------------------------------------------------------------------
# 2. Project paths
# -----------------------------------------------------------------------------

data_dir      <- "data"
model_dir     <- "model_fits"
posterior_dir <- file.path("derived", "posterior_predictions")
table_dir <- file.path("derived", "tables")

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------------------------------------------------------
# 3. Functions
# -----------------------------------------------------------------------------

source(file.path("scripts", "summarize_functions.R"))
source(file.path("scripts", "table_functions.R"))


# -----------------------------------------------------------------------------
# 4. Load formatted analysis data
# -----------------------------------------------------------------------------

load(file.path(data_dir, "formatted_data.RData"))


# -----------------------------------------------------------------------------
# 5. Load final models
# -----------------------------------------------------------------------------

# Aggregate models
m_agg_season <- readRDS(
  file.path(model_dir, "m_agg_season_final.rds")
)

m_agg_env <- readRDS(
  file.path(model_dir, "m_agg_env_final.rds")
)

models_aggregate <- list(
  season_total    = m_agg_season,
  env_conditional = m_agg_env
)


# Frequency-group conditional environmental models
m_low <- readRDS(
  file.path(model_dir, "m_low_final.rds")
)

m_mid <- readRDS(
  file.path(model_dir, "m_mid_final.rds")
)

m_high <- readRDS(
  file.path(model_dir, "m_high_final.rds")
)

models_frequency <- list(
  low  = m_low,
  mid  = m_mid,
  high = m_high
)


# -----------------------------------------------------------------------------
# 6. Load posterior predictions
# -----------------------------------------------------------------------------

# Temperature
post_temp_low <- readRDS(
  file.path(posterior_dir, "temp_low_draws.rds")
)

post_temp_mid <- readRDS(
  file.path(posterior_dir, "temp_mid_draws.rds")
)

post_temp_high <- readRDS(
  file.path(posterior_dir, "temp_high_draws.rds")
)


# Relative humidity
post_rh_low <- readRDS(
  file.path(posterior_dir, "rh_low_draws.rds")
)

post_rh_mid <- readRDS(
  file.path(posterior_dir, "rh_mid_draws.rds")
)

post_rh_high <- readRDS(
  file.path(posterior_dir, "rh_high_draws.rds")
)


# Lunar illumination
post_lunar_low <- readRDS(
  file.path(posterior_dir, "lunar_low_draws.rds")
)

post_lunar_mid <- readRDS(
  file.path(posterior_dir, "lunar_mid_draws.rds")
)

post_lunar_high <- readRDS(
  file.path(posterior_dir, "lunar_high_draws.rds")
)


# -----------------------------------------------------------------------------
# 7. Environmental contrast summaries
# -----------------------------------------------------------------------------
#
# Reconstruct the manuscript contrasts directly from the saved posterior
# predictions so that Table 2 is completely reproducible from saved outputs.
# -----------------------------------------------------------------------------

# --- Temperature: observed 10th to 90th percentile --------------------------

temp_quantiles <- quantile(
  agg_data$avetemp,
  probs = c(0.10, 0.90),
  na.rm = TRUE
)

temp_grid <- post_temp_low |>
  distinct(avetemp) |>
  pull(avetemp)

temp_check_values <- map_dbl(
  temp_quantiles,
  \(x) temp_grid[which.min(abs(temp_grid - x))]
)

temp_contrasts <- bind_rows(
  low  = post_temp_low,
  mid  = post_temp_mid,
  high = post_temp_high,
  .id = "frequency_group"
) |>
  ungroup() |>
  filter(avetemp %in% temp_check_values) |>
  select(
    frequency_group,
    .draw,
    avetemp,
    .epred
  ) |>
  mutate(
    endpoint = case_when(
      avetemp == min(temp_check_values) ~ "q10",
      avetemp == max(temp_check_values) ~ "q90"
    )
  ) |>
  select(-avetemp) |>
  pivot_wider(
    names_from = endpoint,
    values_from = .epred
  ) |>
  mutate(
    difference = q90 - q10,
    ratio = q90 / q10
  )

temp_contrast_summary <- temp_contrasts |>
  group_by(frequency_group) |>
  summarise(
    q10_median = median(q10),
    q90_median = median(q90),
    
    diff_median = median(difference),
    diff_lower = quantile(difference, 0.055),
    diff_upper = quantile(difference, 0.945),
    
    ratio_median = median(ratio),
    ratio_lower = quantile(ratio, 0.055),
    ratio_upper = quantile(ratio, 0.945),
    
    prob_increase = mean(difference > 0),
    
    .groups = "drop"
  )


# --- Relative humidity: observed 10th to 90th percentile --------------------

rh_quantiles <- quantile(
  agg_data$averh,
  probs = c(0.10, 0.90),
  na.rm = TRUE
)

rh_grid <- post_rh_low |>
  distinct(averh) |>
  pull(averh)

rh_check_values <- map_dbl(
  rh_quantiles,
  \(x) rh_grid[which.min(abs(rh_grid - x))]
)

rh_contrasts <- bind_rows(
  low  = post_rh_low,
  mid  = post_rh_mid,
  high = post_rh_high,
  .id = "frequency_group"
) |>
  ungroup() |>
  filter(averh %in% rh_check_values) |>
  select(
    frequency_group,
    .draw,
    averh,
    .epred
  ) |>
  mutate(
    endpoint = case_when(
      averh == min(rh_check_values) ~ "q10",
      averh == max(rh_check_values) ~ "q90"
    )
  ) |>
  select(-averh) |>
  pivot_wider(
    names_from = endpoint,
    values_from = .epred
  ) |>
  mutate(
    difference = q90 - q10,
    ratio = q90 / q10
  )

rh_contrast_summary <- rh_contrasts |>
  group_by(frequency_group) |>
  summarise(
    q10_median = median(q10),
    q90_median = median(q90),
    
    diff_median = median(difference),
    diff_lower = quantile(difference, 0.055),
    diff_upper = quantile(difference, 0.945),
    
    ratio_median = median(ratio),
    ratio_lower = quantile(ratio, 0.055),
    ratio_upper = quantile(ratio, 0.945),
    
    prob_decrease = mean(difference < 0),
    
    .groups = "drop"
  )


# --- Lunar illumination: new moon to full moon ------------------------------

lunar_targets <- c(0, 1)

lunar_grid <- post_lunar_low |>
  distinct(Lunar) |>
  pull(Lunar)

lunar_check_values <- map_dbl(
  lunar_targets,
  \(x) lunar_grid[which.min(abs(lunar_grid - x))]
)

lunar_contrasts <- bind_rows(
  low  = post_lunar_low,
  mid  = post_lunar_mid,
  high = post_lunar_high,
  .id = "frequency_group"
) |>
  ungroup() |>
  filter(Lunar %in% lunar_check_values) |>
  select(
    frequency_group,
    .draw,
    Lunar,
    .epred
  ) |>
  mutate(
    endpoint = case_when(
      Lunar == min(lunar_check_values) ~ "lunar_000",
      Lunar == max(lunar_check_values) ~ "lunar_100"
    )
  ) |>
  select(-Lunar) |>
  pivot_wider(
    names_from = endpoint,
    values_from = .epred
  ) |>
  mutate(
    diff_00_100 = lunar_100 - lunar_000,
    ratio_00_100 = lunar_100 / lunar_000
  )

lunar_contrast_summary <- lunar_contrasts |>
  group_by(frequency_group) |>
  summarise(
    diff_00_100_median = median(diff_00_100),
    diff_00_100_lower = quantile(diff_00_100, 0.055),
    diff_00_100_upper = quantile(diff_00_100, 0.945),
    
    prob_00_100_increase = mean(diff_00_100 > 0),
    
    ratio_00_100_median = median(ratio_00_100),
    ratio_00_100_lower = quantile(ratio_00_100, 0.055),
    ratio_00_100_upper = quantile(ratio_00_100, 0.945),
    
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 8. Main-text tables
# -----------------------------------------------------------------------------

table_1 <- make_table_1(
  agg_data = agg_data,
  cdata_species = cdata_species
)

saveRDS(
  table_1,
  file.path(table_dir, "table_1.rds")
)

table_2 <- make_table_2(
  ref_data = agg_data,
  posterior_dir = posterior_dir,
  table_dir = table_dir
)

saveRDS(
  table_2,
  file.path(table_dir, "table_2.rds")
)


# -----------------------------------------------------------------------------
# 9. Supplementary tables
# -----------------------------------------------------------------------------

table_s1 <- make_table_s1(
  models_aggregate = models_aggregate,
  models_frequency = models_frequency
)

table_s2 <- make_table_s2(
  models_aggregate = models_aggregate
)

table_s3 <- make_table_s3(
  models_frequency = models_frequency
)

table_s4 <- make_table_s4(
  cdata_species = cdata_species
)

table_s5 <- make_table_s5(
  agg_data = agg_data
)

table_s6 <- make_table_s6()

supp_tables <- list(
  table_s1 = table_s1,
  table_s2 = table_s2,
  table_s3 = table_s3,
  table_s4 = table_s4,
  table_s5 = table_s5,
  table_s6 = table_s6
)

# Fail clearly if the diagnostic extraction did not return finite R-hat/ESS values.
diag_check <- bind_rows(table_s2, table_s3)

if (any(!is.finite(diag_check$`Maximum R-hat`)) ||
    any(!is.finite(diag_check$`Minimum bulk ESS`)) ||
    any(!is.finite(diag_check$`Minimum tail ESS`))) {
  stop(
    "Non-finite R-hat or ESS values remain in Tables S2/S3. ",
    "Check model_diagnostics() before saving supplementary tables.",
    call. = FALSE
  )
}

print(table_s1)
print(table_s2)
print(table_s3)

iwalk(
  supp_tables,
  ~ saveRDS(.x, file.path(table_dir, paste0(.y, ".rds")))
)


# -----------------------------------------------------------------------------
# 10. Finished
# -----------------------------------------------------------------------------

message(
  "Main-text and supplementary tables written to: ",
  normalizePath(table_dir)
)
