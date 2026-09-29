# =============================================================================
# Supplementary-table generators
# =============================================================================
# Each `make_table_sX()` function builds the corresponding flextable, saves it
# as a Word document, and returns the flextable invisibly for previewing.
#
# Depends on: summarize_functions.R (extract_model_summary, prettify_var,
# build_spec_row).
# Required packages (load in test.R): tidyverse, flextable, officer, posterior, brms
# =============================================================================

# =============================================================================
# Main-text table generators
# =============================================================================


# -----------------------------------------------------------------------------
# Table 1: Study design and environmental conditions
# -----------------------------------------------------------------------------

make_table_1 <- function(agg_data, cdata_species) {
  
  tibble(
    Characteristic = c(
      "Sites",
      "Sampling nights",
      "Site-nights",
      "Bat species",
      "Sampling period",
      "Day of year",
      "Nightly mean temperature (°C)",
      "Nightly mean relative humidity (%)",
      "Lunar illumination"
    ),
    
    Value = c(
      format(n_distinct(agg_data$Site), big.mark = ","),
      format(n_distinct(agg_data$NiteDayDate), big.mark = ","),
      format(nrow(agg_data), big.mark = ","),
      format(n_distinct(cdata_species$Species), big.mark = ","),
      
      paste0(
        format(min(agg_data$NiteDayDate, na.rm = TRUE), "%B %d"),
        "–",
        format(max(agg_data$NiteDayDate, na.rm = TRUE), "%B %d")
      ),
      
      paste0(
        min(agg_data$doy, na.rm = TRUE),
        "–",
        max(agg_data$doy, na.rm = TRUE)
      ),
      
      sprintf(
        "%.1f–%.1f",
        min(agg_data$avetemp, na.rm = TRUE),
        max(agg_data$avetemp, na.rm = TRUE)
      ),
      
      sprintf(
        "%.1f–%.1f",
        min(agg_data$averh, na.rm = TRUE),
        max(agg_data$averh, na.rm = TRUE)
      ),
      
      sprintf(
        "%.2f–%.2f",
        min(agg_data$Lunar, na.rm = TRUE),
        max(agg_data$Lunar, na.rm = TRUE)
      )
    )
  )
}


# -----------------------------------------------------------------------------
# Table 2: Environmental contrasts by frequency group
# -----------------------------------------------------------------------------

make_table_2 <- function(
    ref_data,
    posterior_dir = "derived/posterior_predictions",
    table_dir = "derived/tables"
) {
  
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------------------
  # Helper: read posterior draws
  # ---------------------------------------------------------------------------
  
  read_group_draws <- function(predictor) {
    
    bind_rows(
      low  = readRDS(file.path(
        posterior_dir,
        paste0(predictor, "_low_draws.rds")
      )),
      mid  = readRDS(file.path(
        posterior_dir,
        paste0(predictor, "_mid_draws.rds")
      )),
      high = readRDS(file.path(
        posterior_dir,
        paste0(predictor, "_high_draws.rds")
      )),
      .id = "frequency_group"
    ) |>
      ungroup()
  }
  
  
  # ---------------------------------------------------------------------------
  # Temperature: observed q10 -> q90
  # ---------------------------------------------------------------------------
  
  temp_draws <- read_group_draws("temp")
  
  temp_quantiles <- quantile(
    ref_data$avetemp,
    probs = c(0.10, 0.90),
    na.rm = TRUE
  )
  
  
  # Use nearest values represented on prediction grid
  temp_grid <- sort(unique(temp_draws$avetemp))
  
  temp_endpoints <- map_dbl(
    temp_quantiles,
    ~ temp_grid[which.min(abs(temp_grid - .x))]
  )
  
  temp_contrasts <- temp_draws |>
    filter(avetemp %in% temp_endpoints) |>
    select(
      frequency_group,
      .draw,
      avetemp,
      .epred
    ) |>
    mutate(
      endpoint = if_else(
        avetemp == min(temp_endpoints),
        "low",
        "high"
      )
    ) |>
    select(-avetemp) |>
    pivot_wider(
      names_from = endpoint,
      values_from = .epred
    ) |>
    mutate(
      difference = high - low,
      ratio = high / low
    )
  
  temp_tab <- temp_contrasts |>
    group_by(frequency_group) |>
    summarise(
      predicted_low = median(low),
      predicted_high = median(high),
      
      diff_median = median(difference),
      diff_lower = quantile(difference, 0.055),
      diff_upper = quantile(difference, 0.945),
      
      ratio_median = median(ratio),
      ratio_lower = quantile(ratio, 0.055),
      ratio_upper = quantile(ratio, 0.945),
      
      prob_increase = mean(difference > 0),
      
      .groups = "drop"
    ) |>
    mutate(
      predictor = "Temperature",
      contrast = "10th–90th percentile"
    )
  
  
  # ---------------------------------------------------------------------------
  # Relative humidity: observed q10 -> q90
  # ---------------------------------------------------------------------------
  
  rh_draws <- read_group_draws("rh")
  
  rh_quantiles <- quantile(
    ref_data$averh,
    probs = c(0.10, 0.90),
    na.rm = TRUE
  )
  
  rh_grid <- sort(unique(rh_draws$averh))
  
  rh_endpoints <- map_dbl(
    rh_quantiles,
    ~ rh_grid[which.min(abs(rh_grid - .x))]
  )
  
  rh_contrasts <- rh_draws |>
    filter(averh %in% rh_endpoints) |>
    select(
      frequency_group,
      .draw,
      averh,
      .epred
    ) |>
    mutate(
      endpoint = if_else(
        averh == min(rh_endpoints),
        "low",
        "high"
      )
    ) |>
    select(-averh) |>
    pivot_wider(
      names_from = endpoint,
      values_from = .epred
    ) |>
    mutate(
      difference = high - low,
      ratio = high / low
    )
  
  rh_tab <- rh_contrasts |>
    group_by(frequency_group) |>
    summarise(
      predicted_low = median(low),
      predicted_high = median(high),
      
      diff_median = median(difference),
      diff_lower = quantile(difference, 0.055),
      diff_upper = quantile(difference, 0.945),
      
      ratio_median = median(ratio),
      ratio_lower = quantile(ratio, 0.055),
      ratio_upper = quantile(ratio, 0.945),
      
      prob_increase = mean(difference > 0),
      
      .groups = "drop"
    ) |>
    mutate(
      predictor = "Relative humidity",
      contrast = "10th–90th percentile"
    )
  
  
  # ---------------------------------------------------------------------------
  # Lunar illumination: new moon -> full moon
  # ---------------------------------------------------------------------------
  
  lunar_draws <- read_group_draws("lunar")
  
  lunar_contrasts <- lunar_draws |>
    filter(Lunar %in% c(0, 1)) |>
    select(
      frequency_group,
      .draw,
      Lunar,
      .epred
    ) |>
    mutate(
      endpoint = if_else(
        Lunar == 0,
        "low",
        "high"
      )
    ) |>
    select(-Lunar) |>
    pivot_wider(
      names_from = endpoint,
      values_from = .epred
    ) |>
    mutate(
      difference = high - low,
      ratio = high / low
    )
  
  lunar_tab <- lunar_contrasts |>
    group_by(frequency_group) |>
    summarise(
      predicted_low = median(low),
      predicted_high = median(high),
      
      diff_median = median(difference),
      diff_lower = quantile(difference, 0.055),
      diff_upper = quantile(difference, 0.945),
      
      ratio_median = median(ratio),
      ratio_lower = quantile(ratio, 0.055),
      ratio_upper = quantile(ratio, 0.945),
      
      prob_increase = mean(difference > 0),
      
      .groups = "drop"
    ) |>
    mutate(
      predictor = "Lunar illumination",
      contrast = "New–full moon"
    )
  
  
  # ---------------------------------------------------------------------------
  # Assemble manuscript table
  # ---------------------------------------------------------------------------
  
  table_2 <- bind_rows(
    temp_tab,
    rh_tab,
    lunar_tab
  ) |>
    mutate(
      predictor = factor(
        predictor,
        levels = c(
          "Temperature",
          "Relative humidity",
          "Lunar illumination"
        )
      ),
      frequency_group = factor(
        frequency_group,
        levels = c("low", "mid", "high")
      ),
      `Frequency group` = recode(
        as.character(frequency_group),
        low = "Low",
        mid = "Mid",
        high = "High"
      ),
      `Predicted low` = sprintf("%.1f", predicted_low),
      `Predicted high` = sprintf("%.1f", predicted_high),
      `Difference (89% CrI)` = sprintf(
        "%.1f (%.1f, %.1f)",
        diff_median,
        diff_lower,
        diff_upper
      ),
      `Ratio (89% CrI)` = sprintf(
        "%.2f (%.2f, %.2f)",
        ratio_median,
        ratio_lower,
        ratio_upper
      ),
      `P(Δ > 0)` = sprintf(
        "%.3f",
        prob_increase
      )
    ) |>
    arrange(
      predictor,
      frequency_group
    ) |>
    transmute(
      Predictor = as.character(predictor),
      `Frequency group`,
      Contrast = contrast,
      `Predicted low`,
      `Predicted high`,
      `Difference (89% CrI)`,
      `Ratio (89% CrI)`,
      `P(Δ > 0)`
    )
  
  saveRDS(
    table_2,
    file.path(table_dir, "table_2.rds")
  )
  
  invisible(table_2)
}

# -----------------------------------------------------------------------------
# Generic summary-table builder used by Tables S2, S3a–c
# -----------------------------------------------------------------------------

build_summary_table <- function(data, group_col, caption_text) {
  data <- data |>
    mutate(Parameter = prettify_var(variable)) |>
    select(all_of(group_col), Parameter,
           Median, `Lower 89%`, `Upper 89%`, Rhat, `Bulk ESS`)

  flextable(data) |>
    bold(part = "header") |>
    align(j = c("Median", "Lower 89%", "Upper 89%", "Rhat", "Bulk ESS"),
          align = "center", part = "all") |>
    merge_v(j = group_col) |>
    valign(j = group_col, valign = "top") |>
    set_table_properties(layout = "autofit") |>
    set_caption(
      caption          = caption_text,
      fp_p             = fp_par(text.align = "left", padding = 3),
      align_with_table = FALSE
    )
}


# -----------------------------------------------------------------------------
# Supplementary tables: final analysis
# -----------------------------------------------------------------------------
# These functions return plain tibbles for direct use in Quarto. The driver
# saves each object as an RDS file under derived/tables/.


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

model_formula_text <- function(model) {
  f <- formula(model)
  
  if (inherits(f, "brmsformula")) {
    f <- f$formula
  }
  
  paste(deparse(f), collapse = " ")
}

model_family_text <- function(model) {
  fam <- model$family
  paste0(fam$family, " (", fam$link, " link)")
}

finite_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

finite_min <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  min(x)
}

model_diagnostics <- function(model, model_name, max_treedepth = 12) {
  
  # ---------------------------------------------------------------------------
  # Posterior draws and convergence diagnostics
  # ---------------------------------------------------------------------------
  
  draws <- posterior::as_draws_array(model)
  
  diag <- posterior::summarise_draws(
    draws,
    "rhat",
    "ess_bulk",
    "ess_tail"
  )
  
  # Some parameters can legitimately have undefined diagnostics
  # (for example, constants). Restrict extrema to finite values.
  diag_finite <- diag |>
    dplyr::filter(
      is.finite(.data$rhat),
      is.finite(.data$ess_bulk),
      is.finite(.data$ess_tail)
    )
  
  if (nrow(diag_finite) == 0) {
    stop(
      "No finite R-hat/ESS diagnostics found for ",
      model_name,
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Stan sampler diagnostics
  # ---------------------------------------------------------------------------
  
  np <- brms::nuts_params(model)
  
  divergences <- sum(
    np$Parameter == "divergent__" &
      np$Value == 1
  )
  
  max_td <- sum(
    np$Parameter == "treedepth__" &
      np$Value >= max_treedepth
  )
  
  # ---------------------------------------------------------------------------
  # Negative-binomial shape
  # ---------------------------------------------------------------------------
  
  ddf <- posterior::as_draws_df(model)
  
  shape_name <- intersect(
    c("shape", "b_shape"),
    names(ddf)
  )
  
  if (length(shape_name) == 0) {
    
    shape_med <- NA_real_
    shape_lo  <- NA_real_
    shape_hi  <- NA_real_
    
  } else {
    
    shape_draws <- ddf[[shape_name[1]]]
    
    shape_med <- median(shape_draws)
    shape_lo  <- unname(quantile(shape_draws, 0.055))
    shape_hi  <- unname(quantile(shape_draws, 0.945))
  }
  
  # ---------------------------------------------------------------------------
  # Output
  # ---------------------------------------------------------------------------
  
  tibble(
    Model = model_name,
    `Site-nights` = nobs(model),
    `Post-warmup draws` = posterior::ndraws(draws),
    `Maximum R-hat` = max(diag_finite$rhat),
    `Minimum bulk ESS` = min(diag_finite$ess_bulk),
    `Minimum tail ESS` = min(diag_finite$ess_tail),
    Divergences = divergences,
    `Max-treedepth transitions` = max_td,
    `NB shape median` = shape_med,
    `NB shape lower 89%` = shape_lo,
    `NB shape upper 89%` = shape_hi
  )
}

format_diagnostic_table <- function(x) {
  x |>
    mutate(
      `Maximum R-hat` = round(`Maximum R-hat`, 2),
      `Minimum bulk ESS` = round(`Minimum bulk ESS`),
      `Minimum tail ESS` = round(`Minimum tail ESS`),
      across(
        starts_with("NB shape"),
        ~ round(.x, 2)
      )
    )
}

# -----------------------------------------------------------------------------
# Table S1: Final model specifications
# -----------------------------------------------------------------------------

make_table_s1 <- function(models_aggregate, models_frequency) {

  model_rows <- bind_rows(
    tibble(
      Model = "Aggregate seasonal",
      Response = "Aggregate activity",
      model = list(models_aggregate$season_total)
    ),
    tibble(
      Model = "Aggregate environmental",
      Response = "Aggregate activity",
      model = list(models_aggregate$env_conditional)
    ),
    tibble(
      Model = c(
        "Low-frequency environmental",
        "Mid-frequency environmental",
        "High-frequency environmental"
      ),
      Response = c(
        "Low-frequency activity",
        "Mid-frequency activity",
        "High-frequency activity"
      ),
      model = list(
        models_frequency$low,
        models_frequency$mid,
        models_frequency$high
      )
    )
  )

  model_rows |>
    mutate(
      Formula = map_chr(model, model_formula_text),
      Likelihood = map_chr(model, model_family_text),
      `Site-nights` = map_int(model, nobs)
    ) |>
    select(Model, Response, Formula, Likelihood, `Site-nights`)
}


# -----------------------------------------------------------------------------
# Table S2: Aggregate-model sampling diagnostics
# -----------------------------------------------------------------------------

make_table_s2 <- function(models_aggregate) {
  
  bind_rows(
    model_diagnostics(
      models_aggregate$season_total,
      "Aggregate seasonal"
    ),
    model_diagnostics(
      models_aggregate$env_conditional,
      "Aggregate environmental"
    )
  ) |>
    format_diagnostic_table()
}


# -----------------------------------------------------------------------------
# Table S3: Frequency-group model sampling diagnostics
# -----------------------------------------------------------------------------

make_table_s3 <- function(models_frequency) {
  
  bind_rows(
    model_diagnostics(
      models_frequency$low,
      "Low-frequency group"
    ),
    model_diagnostics(
      models_frequency$mid,
      "Mid-frequency group"
    ),
    model_diagnostics(
      models_frequency$high,
      "High-frequency group"
    )
  ) |>
    format_diagnostic_table()
}


# -----------------------------------------------------------------------------
# Table S4: Species composition and frequency classification
# -----------------------------------------------------------------------------

make_table_s4 <- function(cdata_species) {

  species_lookup <- tribble(
    ~Species, ~scientific_name,           ~common_name,                  ~freq_group, ~char_freq_kHz,
    "MYOYUM", "Myotis yumanensis",        "Yuma myotis",                 "high", 48,
    "MYOLUC", "Myotis lucifugus",         "Little brown myotis",         "high", 42,
    "MYOCAL", "Myotis californicus",      "California myotis",           "high", 50,
    "MYOCIL", "Myotis ciliolabrum",       "Western small-footed myotis", "high", 45,
    "MYOVOL", "Myotis volans",            "Long-legged myotis",          "high", 40,
    "PARHES", "Parastrellus hesperus",    "Canyon bat",                  "high", 45,
    "LASBLO", "Lasiurus blossevillii",    "Western red bat",             "high", 42,
    "MYOEVO", "Myotis evotis",            "Long-eared myotis",           "mid",  35,
    "MYOTHY", "Myotis thysanodes",        "Fringed myotis",              "mid",  32,
    "EPTFUS", "Eptesicus fuscus",         "Big brown bat",               "mid",  27,
    "LASNOC", "Lasionycteris noctivagans","Silver-haired bat",           "mid",  26,
    "TADBRA", "Tadarida brasiliensis",    "Mexican free-tailed bat",     "mid",  25,
    "ANTPAL", "Antrozous pallidus",       "Pallid bat",                  "mid",  30,
    "CORTOW", "Corynorhinus townsendii",  "Townsend's big-eared bat",    "mid",  28,
    "LASCIN", "Lasiurus cinereus",        "Hoary bat",                   "low",  20,
    "EUMPER", "Eumops perotis",           "Western mastiff bat",         "low",  12,
    "EUDMAC", "Euderma maculatum",        "Spotted bat",                 "low",  11
  )

  cdata_species |>
    group_by(Species) |>
    summarise(
      `Total intervals` = sum(n, na.rm = TRUE),
      `Site-nights detected` = sum(n > 0, na.rm = TRUE),
      .groups = "drop"
    ) |>
    left_join(species_lookup, by = "Species") |>
    mutate(freq_group = factor(freq_group, levels = c("low", "mid", "high"))) |>
    arrange(freq_group, desc(`Total intervals`)) |>
    transmute(
      `Species code` = Species,
      `Scientific name` = scientific_name,
      `Common name` = common_name,
      `Characteristic frequency (kHz)` = char_freq_kHz,
      `Frequency group` = recode(as.character(freq_group), low = "Low", mid = "Mid", high = "High"),
      `Total intervals`,
      `Site-nights detected`
    )   |>
    mutate(
      `Frequency group` = factor(
        `Frequency group`,
        levels = c("Low", "Mid", "High")
      )
    ) |>
    arrange(
      `Frequency group`,
      `Characteristic frequency (kHz)`
    )
}


# -----------------------------------------------------------------------------
# Table S5: Retained sampling effort by site
# -----------------------------------------------------------------------------

make_table_s5 <- function(agg_data) {
  agg_data |>
    group_by(Site) |>
    summarise(
      `Retained site-nights` = n(),
      `Zero-activity nights` = sum(n == 0, na.rm = TRUE),
      `Total activity intervals` = sum(n, na.rm = TRUE),
      `Mean intervals per night` = mean(n, na.rm = TRUE),
      `Median intervals per night` = median(n, na.rm = TRUE),
      `First retained night` = min(NiteDayDate, na.rm = TRUE),
      `Last retained night` = max(NiteDayDate, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(Site)   |>
    mutate(
      `Mean intervals per night` =
        round(`Mean intervals per night`, 1),
      `Median intervals per night` =
        round(`Median intervals per night`, 1)
    )
}


# -----------------------------------------------------------------------------
# Table S6: Software and computational environment
# -----------------------------------------------------------------------------

make_table_s6 <- function() {
  pkgs <- c(
    "R", "brms", "cmdstanr", "posterior", "tidybayes",
    "mgcv", "tidyverse", "ggplot2", "patchwork", "flextable"
  )

  r_ver <- paste(R.version$major, R.version$minor, sep = ".")

  tibble(
    Software = pkgs,
    Version = map_chr(
      pkgs,
      function(x) {
        if (x == "R") return(r_ver)
        if (!requireNamespace(x, quietly = TRUE)) return(NA_character_)
        as.character(utils::packageVersion(x))
      }
    )
  )
}
