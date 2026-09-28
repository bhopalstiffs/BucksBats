# =============================================================================
# Bucks Bats: posterior summarization helpers
# =============================================================================
# Functions used by the table- and figure-building scripts.
#
# This version is written for the final negative-binomial GAM structure:
#
#   n ~ s(doy) + s(doy, Site, bs = "fs") +
#       s(scale_avetemp) + s(scale_averh) + s(Lunar)
#
# The key change from the original helpers is that predictions from models with
# a site-specific factor smooth are generated for every observed site and then
# averaged WITHIN posterior draw. This retains uncertainty correctly and makes
# the estimand explicit.
#
# Required packages (loaded in test.R):
#   tidyverse, tidybayes, posterior, brms
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Prediction-grid helpers
# -----------------------------------------------------------------------------

# Build a reference grid for a focal predictor.
#
# Non-focal environmental covariates are held at interpretable reference values:
#   * doy: median observed day of year (REFERENCE mode only)
#   * standardized temperature: 0 (observed mean)
#   * standardized RH: 0 (observed mean)
#   * Lunar: median observed value
#
# If a raw/scaled pair is supplied (e.g., avetemp/scale_avetemp), focal_seq is
# interpreted on the RAW scale and both columns are populated.
build_pred_grid <- function(focal, focal_seq, data,
                            scaled_var = NULL, raw_var = NULL) {

  pred_grid <- tibble(
    doy           = median(data$doy, na.rm = TRUE),
    scale_avetemp = 0,
    scale_averh   = 0,
    Lunar         = median(data$Lunar, na.rm = TRUE)
  ) |>
    slice(rep(1, length(focal_seq))) |>
    mutate(grid_id = row_number())

  pred_grid[[focal]] <- focal_seq

  if (!is.null(scaled_var) && !is.null(raw_var)) {
    pred_grid[[raw_var]] <- focal_seq
    pred_grid[[scaled_var]] <-
      (focal_seq - mean(data[[raw_var]], na.rm = TRUE)) /
      sd(data[[raw_var]], na.rm = TRUE)
  }

  pred_grid
}


# Common manuscript endpoints for environmental contrasts.
#
# IMPORTANT: supply the SAME reference_data object for low-, mid-, and
# high-frequency models. This guarantees that all frequency groups are compared
# over identical environmental ranges.
#
# Temperature and RH use the pooled 10th and 90th percentiles; lunar phase uses
# its natural endpoints, 0 (new moon) and 1 (full moon).
get_manuscript_endpoints <- function(reference_data) {

  tibble(
    predictor = c("temperature", "relative_humidity", "lunar"),
    value_1 = c(
      unname(quantile(reference_data$avetemp, 0.10, na.rm = TRUE)),
      unname(quantile(reference_data$averh,   0.10, na.rm = TRUE)),
      0
    ),
    value_2 = c(
      unname(quantile(reference_data$avetemp, 0.90, na.rm = TRUE)),
      unname(quantile(reference_data$averh,   0.90, na.rm = TRUE)),
      1
    )
  )
}


# -----------------------------------------------------------------------------
# 2. Posterior expected predictions
# -----------------------------------------------------------------------------

# Site-averaged posterior expected predictions at a REFERENCE day of year.
#
# This retains the original behavior and is useful when a prediction is meant to
# be conditional on a representative point in the season. Environmental
# manuscript figures/contrasts should instead use season = "average" below.
get_site_averaged_predictions_reference <- function(model, data, focal, focal_seq,
                                                     scaled_var = NULL,
                                                     raw_var = NULL) {

  pred_grid <- build_pred_grid(
    focal      = focal,
    focal_seq  = focal_seq,
    data       = data,
    scaled_var = scaled_var,
    raw_var    = raw_var
  )

  site_levels <- levels(factor(data$Site))

  site_grid <- tidyr::crossing(
    pred_grid,
    Site = site_levels
  ) |>
    mutate(Site = factor(Site, levels = site_levels))

  site_grid |>
    add_epred_draws(model, re_formula = NULL) |>
    group_by(.draw, grid_id) |>
    summarise(.epred = mean(.epred), .groups = "drop") |>
    left_join(pred_grid, by = "grid_id")
}


# Site- and season-averaged posterior expected predictions.
#
# PRIMARY estimand for environmental response curves and contrasts.
#
# For each focal environmental value, predictions are generated at every
# OBSERVED Site x DOY combination in the analysis data. The averaging is done in
# two stages WITHIN each posterior draw:
#   1. average over the observed DOY distribution within each site;
#   2. average those site-specific means equally across sites.
#
# Thus, dates actually sampled at a site define that site's seasonal average,
# while sites with more sampled nights do not receive greater weight in the
# final ecological summary.
get_site_season_averaged_predictions <- function(model, data, focal, focal_seq,
                                                 scaled_var = NULL,
                                                 raw_var = NULL) {

  if (focal == "doy") {
    stop("Season-averaged predictions are for environmental predictors; use season = 'reference' when focal = 'doy'.")
  }

  pred_grid <- build_pred_grid(
    focal      = focal,
    focal_seq  = focal_seq,
    data       = data,
    scaled_var = scaled_var,
    raw_var    = raw_var
  )

  site_levels <- levels(factor(data$Site))

  observed_site_doy <- data |>
    transmute(
      Site = factor(Site, levels = site_levels),
      doy  = doy
    ) |>
    filter(!is.na(Site), !is.na(doy)) |>
    distinct(Site, doy)

  # Cross each focal value with the Site x DOY combinations that were actually
  # sampled. Non-focal environmental predictors remain at reference values.
  season_grid <- tidyr::crossing(
    pred_grid |>
      select(-doy),
    observed_site_doy
  ) |>
    mutate(Site = factor(Site, levels = site_levels))

  season_grid |>
    add_epred_draws(model, re_formula = NULL) |>
    group_by(.draw, grid_id, Site) |>
    summarise(.site_epred = mean(.epred), .groups = "drop") |>
    group_by(.draw, grid_id) |>
    summarise(.epred = mean(.site_epred), .groups = "drop") |>
    left_join(pred_grid |> select(-doy), by = "grid_id")
}


# Back-compatible alias for reference-day site averaging.
get_site_averaged_predictions <- get_site_averaged_predictions_reference


# Population/reference predictions for models that do NOT require Site inside a
# smooth term. Retained for older/reference models only.
get_population_predictions <- function(model, data, focal, focal_seq,
                                       scaled_var = NULL, raw_var = NULL) {

  build_pred_grid(
    focal      = focal,
    focal_seq  = focal_seq,
    data       = data,
    scaled_var = scaled_var,
    raw_var    = raw_var
  ) |>
    add_epred_draws(
      model,
      re_formula = NA,
      allow_new_levels = TRUE
    )
}


# Unified prediction interface.
#
# For final factor-smooth models:
#   prediction = "site_average", season = "average"
# is the manuscript default for ENVIRONMENTAL effects.
#
# season = "reference" retains the older median-DOY behavior and is required
# when focal = "doy".
get_predictions <- function(model, data, focal, focal_seq,
                            scaled_var = NULL, raw_var = NULL,
                            prediction = c("site_average", "population"),
                            season = c("average", "reference")) {

  prediction <- match.arg(prediction)
  season <- match.arg(season)

  if (prediction == "population") {
    return(get_population_predictions(
      model, data, focal, focal_seq,
      scaled_var = scaled_var,
      raw_var    = raw_var
    ))
  }

  if (focal == "doy") {
    season <- "reference"
  }

  if (season == "average") {
    get_site_season_averaged_predictions(
      model, data, focal, focal_seq,
      scaled_var = scaled_var,
      raw_var    = raw_var
    )
  } else {
    get_site_averaged_predictions_reference(
      model, data, focal, focal_seq,
      scaled_var = scaled_var,
      raw_var    = raw_var
    )
  }
}


# Summarize posterior expected predictions at each focal value.
summarize_predictions <- function(model, data, focal, focal_seq,
                                  scaled_var = NULL, raw_var = NULL,
                                  prediction = c("site_average", "population"),
                                  season = c("average", "reference"),
                                  .width = 0.89) {

  prediction <- match.arg(prediction)
  season <- match.arg(season)

  get_predictions(
    model      = model,
    data       = data,
    focal      = focal,
    focal_seq  = focal_seq,
    scaled_var = scaled_var,
    raw_var    = raw_var,
    prediction = prediction,
    season     = season
  ) |>
    group_by(grid_id) |>
    median_qi(.epred, .width = .width) |>
    ungroup()
}


# Range of posterior interval bounds across a focal sequence. Used to construct
# shared y-axes across panels.
get_epred_range <- function(model, data, focal, focal_seq,
                            scaled_var = NULL, raw_var = NULL,
                            prediction = c("site_average", "population"),
                            season = c("average", "reference"),
                            .width = 0.89) {

  prediction <- match.arg(prediction)
  season <- match.arg(season)

  preds_sum <- summarize_predictions(
    model      = model,
    data       = data,
    focal      = focal,
    focal_seq  = focal_seq,
    scaled_var = scaled_var,
    raw_var    = raw_var,
    prediction = prediction,
    season     = season,
    .width     = .width
  )

  range(c(preds_sum$.lower, preds_sum$.upper), na.rm = TRUE)
}


# -----------------------------------------------------------------------------
# 3. Site-specific seasonal trajectories
# -----------------------------------------------------------------------------

# Draw expected seasonal trajectories separately for each observed site.
# Environmental covariates are held at their reference values while doy varies.
# This is intended for the figure showing heterogeneity in seasonal phenology.
get_site_season_predictions <- function(model, data, doy_seq) {

  pred_grid <- build_pred_grid(
    focal     = "doy",
    focal_seq = doy_seq,
    data      = data
  )

  site_levels <- levels(factor(data$Site))

  tidyr::crossing(
    pred_grid,
    Site = site_levels
  ) |>
    mutate(Site = factor(Site, levels = site_levels)) |>
    add_epred_draws(model, re_formula = NULL)
}


# Summarize each site's seasonal trajectory.
summarize_site_season_predictions <- function(model, data, doy_seq,
                                              .width = 0.89) {

  get_site_season_predictions(model, data, doy_seq) |>
    group_by(Site, grid_id, doy) |>
    median_qi(.epred, .width = .width) |>
    ungroup()
}


# -----------------------------------------------------------------------------
# 4. Posterior contrasts from effect curves
# -----------------------------------------------------------------------------

# Contrast two focal values on the expected-response scale.
# Returns posterior draws of the absolute difference and response ratio.
# Useful for manuscript-ready summaries that are more interpretable than GAM
# basis coefficients.
get_epred_contrast <- function(model, data, focal, values,
                               scaled_var = NULL, raw_var = NULL,
                               prediction = c("site_average", "population"),
                               season = c("average", "reference")) {

  prediction <- match.arg(prediction)
  season <- match.arg(season)

  if (length(values) != 2) {
    stop("`values` must contain exactly two focal values.")
  }

  preds <- get_predictions(
    model      = model,
    data       = data,
    focal      = focal,
    focal_seq  = values,
    scaled_var = scaled_var,
    raw_var    = raw_var,
    prediction = prediction,
    season     = season
  ) |>
    select(.draw, grid_id, .epred) |>
    tidyr::pivot_wider(
      names_from  = grid_id,
      values_from = .epred,
      names_prefix = "pred_"
    )

  preds |>
    transmute(
      .draw,
      value_1    = values[1],
      value_2    = values[2],
      epred_1    = pred_1,
      epred_2    = pred_2,
      difference = pred_2 - pred_1,
      ratio      = pred_2 / pred_1
    )
}


# Summarize an expected-response contrast with an 89% credible interval.
summarize_epred_contrast <- function(model, data, focal, values,
                                     scaled_var = NULL, raw_var = NULL,
                                     prediction = c("site_average", "population"),
                                     season = c("average", "reference"),
                                     .width = 0.89) {

  prediction <- match.arg(prediction)
  season <- match.arg(season)

  draws <- get_epred_contrast(
    model      = model,
    data       = data,
    focal      = focal,
    values     = values,
    scaled_var = scaled_var,
    raw_var    = raw_var,
    prediction = prediction,
    season     = season
  )

  bind_rows(
    draws |>
      median_qi(difference, .width = .width) |>
      mutate(contrast = "difference"),
    draws |>
      median_qi(ratio, .width = .width) |>
      mutate(contrast = "ratio")
  )
}


# -----------------------------------------------------------------------------
# 5. Model-parameter summaries
# -----------------------------------------------------------------------------

# Extract a diagnostic/model-specification summary from a brms model.
#
# GAM basis coefficients are retained for completeness, but they should not be
# interpreted individually as biological effects. Biological inference should
# come from posterior expected-response curves and contrasts.
#
# Includes the negative-binomial shape parameter when present.
extract_model_summary <- function(model, model_name) {

  draws <- as_draws_df(model)
  vars  <- variables(draws)

  vars_keep <- vars[
    (
      grepl("^(b_|sds_|sd_)", vars) |
        vars %in% c("shape", "b_shape")
    ) &
      !grepl("(^zs_|^r_|^z_|^prior_|^lprior|^lp__)", vars)
  ]

  summarise_draws(
    subset_draws(draws, variable = vars_keep),
    median,
    ~ quantile(.x, probs = 0.055),
    ~ quantile(.x, probs = 0.945),
    rhat,
    ess_bulk,
    ess_tail
  ) |>
    rename(
      Median      = median,
      `Lower 89%` = `5.5%`,
      `Upper 89%` = `94.5%`,
      Rhat        = rhat,
      `Bulk ESS`  = ess_bulk,
      `Tail ESS`  = ess_tail
    ) |>
    mutate(
      Median      = round(Median, 2),
      `Lower 89%` = round(`Lower 89%`, 2),
      `Upper 89%` = round(`Upper 89%`, 2),
      Rhat        = round(Rhat, 3),
      `Bulk ESS`  = round(`Bulk ESS`),
      `Tail ESS`  = round(`Tail ESS`),
      Model       = model_name
    ) |>
    select(
      Model, variable, Median,
      `Lower 89%`, `Upper 89%`,
      Rhat, `Bulk ESS`, `Tail ESS`
    )
}


# Map brms parameter names to readable display labels.
# Factor-smooth names can vary slightly with brms/mgcv versions, so pattern
# matching is used rather than relying entirely on exact names.
prettify_var <- function(x) {

  case_when(
    x == "b_Intercept" ~ "Intercept",

    x == "b_sdoy_1" ~ "s(doy), basis coefficient",
    x == "b_sscale_avetemp_1" ~ "s(temperature), basis coefficient",
    x == "b_sscale_averh_1" ~ "s(humidity), basis coefficient",
    x == "b_sLunar_1" ~ "s(lunar), basis coefficient",

    x == "sds_sdoy_1" ~ "s(doy), smooth SD",
    grepl("^sds_sdoySite", x) ~ "Site-specific seasonal smooth SD",
    x == "sds_sscale_avetemp_1" ~ "s(temperature), smooth SD",
    x == "sds_sscale_averh_1" ~ "s(humidity), smooth SD",
    x == "sds_sLunar_1" ~ "s(lunar), smooth SD",

    x == "sd_Site__Intercept" ~ "Site random-intercept SD",
    x %in% c("shape", "b_shape") ~ "Negative-binomial shape",

    TRUE ~ x
  )
}


# -----------------------------------------------------------------------------
# 6. Model-specification table helper
# -----------------------------------------------------------------------------

# Build one row of the model-specification table.
# `adapt_delta` and `max_treedepth` are supplied explicitly because they are
# fitting controls rather than posterior parameters.
build_spec_row <- function(model, model_label, dataset_label,
                           adapt_delta = NA_real_,
                           max_treedepth = NA_integer_) {

  fmla <- formula(model)$formula |>
    deparse() |>
    paste(collapse = " ") |>
    gsub("\\s+", " ", x = _)

  fam <- family(model)

  tibble(
    Model               = model_label,
    Formula             = fmla,
    Family              = fam$family,
    Link                = fam$link,
    Dataset             = dataset_label,
    N                   = nobs(model),
    adapt_delta         = adapt_delta,
    max_treedepth       = max_treedepth,
    `Post-warmup draws` = ndraws(model)
  )
}

# -----------------------------------------------------------------------------
# Memory-efficient posterior predictions over a focal grid
# -----------------------------------------------------------------------------

get_predictions_chunked <- function(model,
                                    data,
                                    focal,
                                    focal_seq,
                                    scaled_var = NULL,
                                    raw_var = NULL,
                                    prediction = "site_average",
                                    season = "average",
                                    chunk_size = 5) {
  
  chunks <- split(
    focal_seq,
    ceiling(seq_along(focal_seq) / chunk_size)
  )
  
  out <- vector("list", length(chunks))
  
  for (i in seq_along(chunks)) {
    
    message(
      "Prediction chunk ",
      i, " of ", length(chunks),
      " (", length(chunks[[i]]), " focal values)"
    )
    
    out[[i]] <- get_predictions(
      model      = model,
      data       = data,
      focal      = focal,
      focal_seq  = chunks[[i]],
      scaled_var = scaled_var,
      raw_var    = raw_var,
      prediction = prediction,
      season     = season
    )
    
    # Encourage R to release temporary objects
    gc()
  }
  
  dplyr::bind_rows(out)
}
