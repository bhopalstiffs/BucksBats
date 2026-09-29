# =============================================================================
# Bucks Bats: figure generators
# =============================================================================
# Figures for the final negative-binomial GAM framework.
#
# Primary conditional model:
#   n ~ s(doy) + s(doy, Site, bs = "fs") +
#       s(scale_avetemp) + s(scale_averh) + s(Lunar)
#
# Depends on summarize_functions.R, especially:
#   get_predictions(), get_epred_range(),
#   get_site_season_predictions(), summarize_site_season_predictions()
#
# Required packages (loaded in test.R):
#   tidyverse, tidybayes, patchwork, RColorBrewer, brms, lubridate
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Shared theme and labels
# -----------------------------------------------------------------------------

activity_ylab <- paste0(
  "Expected number of 5-min intervals per night\n",
  "with bat activity"
)

bat_theme <- function() {
  theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.y     = element_text(size = 9),
      axis.title.y    = element_text(size = 15, margin = margin(r = 5)),
      axis.text.x     = element_text(size = 12),
      axis.title.x    = element_text(size = 16, margin = margin(t = 5)),
      plot.title      = element_text(size = 14, face = "bold")
    )
}


# -----------------------------------------------------------------------------
# 2. Low-level effect plot
# -----------------------------------------------------------------------------

# Plot posterior expected activity for one focal predictor.
#
# For final factor-smooth models, prediction = "site_average" should normally be
# used. Predictions are generated for each observed site and averaged equally
# across sites within posterior draw.
make_effect_plot <- function(model, data, focal, focal_seq, xlab,
                             scaled_var = NULL, raw_var = NULL,
                             is_doy = FALSE, ref_date = "2015-01-01",
                             y_limits = NULL,
                             prediction = c("site_average", "population"),
                             widths = c(0.5, 0.8, 0.89)) {

  prediction <- match.arg(prediction)

  predictions <- get_predictions(
    model      = model,
    data       = data,
    focal      = focal,
    focal_seq  = focal_seq,
    scaled_var = scaled_var,
    raw_var    = raw_var,
    prediction = prediction
  )

  if (is_doy) {
    predictions <- predictions |>
      mutate(x_plot = ymd(ref_date) + days(.data[[focal]] - 1))
  } else {
    xvar <- if (!is.null(raw_var)) raw_var else focal
    predictions <- predictions |>
      mutate(x_plot = .data[[xvar]])
  }

  p <- ggplot(predictions, aes(x = x_plot, y = .epred)) +
    stat_lineribbon(
      aes(fill = after_stat(level)),
      .width = widths,
      color  = "gray20"
    ) +
    scale_fill_brewer(palette = "Greys") +
    labs(
      x = if (is_doy) NULL else xlab,
      y = activity_ylab
    ) +
    bat_theme()

  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }

  p
}


# -----------------------------------------------------------------------------
# 3. Aggregate environmental-effects figure
# -----------------------------------------------------------------------------

# Main conditional environmental effects from the final aggregate model.
# Temperature and RH are displayed on their original scales. Non-focal
# environmental covariates are held at the reference values defined in
# build_pred_grid(); predictions are equal-site averages.
make_aggregate_environment_figure <- function(model, data, fig_dir = NULL,
                                              filename = "aggregate_environment_effects.PNG") {

  temp_seq <- seq(min(data$avetemp, na.rm = TRUE),
                  max(data$avetemp, na.rm = TRUE), length.out = 100)
  rh_seq <- seq(min(data$averh, na.rm = TRUE),
                max(data$averh, na.rm = TRUE), length.out = 100)
  lunar_seq <- seq(min(data$Lunar, na.rm = TRUE),
                   max(data$Lunar, na.rm = TRUE), length.out = 100)

  p_temp <- make_effect_plot(
    model, data,
    focal = "scale_avetemp", focal_seq = temp_seq,
    scaled_var = "scale_avetemp", raw_var = "avetemp",
    xlab = "Average Nightly Temperature (\u00B0C)",
    prediction = "site_average"
  )

  p_rh <- make_effect_plot(
    model, data,
    focal = "scale_averh", focal_seq = rh_seq,
    scaled_var = "scale_averh", raw_var = "averh",
    xlab = "Average Nightly Relative Humidity (%)",
    prediction = "site_average"
  ) +
    labs(y = NULL)

  p_lunar <- make_effect_plot(
    model, data,
    focal = "Lunar", focal_seq = lunar_seq,
    xlab = "Lunar Phase (0 = new moon, 1 = full moon)",
    prediction = "site_average"
  ) +
    labs(y = NULL)

  p <- (p_temp + p_rh + p_lunar) +
    plot_annotation(tag_levels = "a")

  if (!is.null(fig_dir)) {
    ggsave(
      filename = file.path(fig_dir, filename),
      plot = p, width = 20, height = 7.5, units = "in", dpi = 300
    )
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# 4. Seasonal figures
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 4.1 Plot a posterior seasonal trajectory
# -----------------------------------------------------------------------------

plot_season_posterior <- function(post,
                                  title = NULL,
                                  ref_date = "2015-01-01",
                                  y_limits = NULL,
                                  show_y_title = TRUE) {
  
  plot_data <- post %>%
    ungroup() %>%
    mutate(
      date = lubridate::ymd(ref_date) + lubridate::days(doy - 1)
    )
  
  p <- ggplot(
    plot_data,
    aes(x = date, y = .epred)
  ) +
    tidybayes::stat_lineribbon(
      .width = c(0.50, 0.80, 0.89),
      aes(fill = after_stat(level)),
      color = "gray15",
      linewidth = 0.9
    ) +
    scale_fill_brewer(
      palette = "Greys",
      direction = -1
    ) +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%b",
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      x = NULL,
      y = if (show_y_title) activity_ylab else NULL,
      title = title
    ) +
    bat_theme() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
  
  if (!is.null(y_limits)) {
    p <- p +
      coord_cartesian(ylim = y_limits)
  }
  
  p
}


# -----------------------------------------------------------------------------
# 4.2 Plot site-specific seasonal phenology
# -----------------------------------------------------------------------------

plot_site_season_posterior <- function(post_sites,
                                       post_average,
                                       title = "Site-specific seasonal phenology",
                                       ref_date = "2015-01-01") {
  
  # Summarize each site's posterior trajectory
  site_sum <- post_sites %>%
    ungroup() %>%
    mutate(
      date = lubridate::ymd(ref_date) + lubridate::days(doy - 1)
    ) %>%
    group_by(Site, doy, date) %>%
    median_qi(
      .epred,
      .width = 0.89
    ) %>%
    ungroup()
  
  # Equal-site posterior average
  avg_data <- post_average %>%
    ungroup() %>%
    mutate(
      date = lubridate::ymd(ref_date) + lubridate::days(doy - 1)
    )
  
  ggplot() +
    
    # Individual site posterior medians
    geom_line(
      data = site_sum,
      aes(
        x = date,
        y = .epred,
        group = Site
      ),
      color = "gray55",
      linewidth = 0.45,
      alpha = 0.70
    ) +
    
    # Equal-site average + uncertainty
    tidybayes::stat_lineribbon(
      data = avg_data,
      aes(
        x = date,
        y = .epred,
        fill = after_stat(level)
      ),
      .width = c(0.50, 0.80, 0.89),
      color = "gray10",
      linewidth = 1.0
    ) +
    
    scale_fill_brewer(
      palette = "Greys",
      direction = -1
    ) +
    
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%b",
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    labs(
      x = NULL,
      y = activity_ylab,
      title = title
    ) +
    
    bat_theme() +
    
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      
      # Give the y title a little breathing room
      axis.title.y = element_text(
        margin = margin(r = 12)
      )
    )
}


# -----------------------------------------------------------------------------
# 4.3 Figure 3
# -----------------------------------------------------------------------------

make_figure3 <- function(post_total,
                         post_conditional,
                         post_sites,
                         fig_dir = NULL,
                         filename = "Figure3_seasonality.PNG") {
  
  # Common y limits for panels a and b
  upper_y <- max(
    quantile(post_total$.epred, 0.995, na.rm = TRUE),
    quantile(post_conditional$.epred, 0.995, na.rm = TRUE)
  )
  
  y_limits <- c(0, upper_y)
  
  # Panel a
  p_a <- plot_season_posterior(
    post_total,
    title = "Total seasonal pattern",
    y_limits = y_limits,
    show_y_title = TRUE
  )
  
  # Panel b
  p_b <- plot_season_posterior(
    post_conditional,
    title = "Conditional seasonal pattern",
    y_limits = y_limits,
    show_y_title = FALSE
  )
  
  # Panel c
  p_c <- plot_site_season_posterior(
    post_sites,
    post_total,
    title = "Site-specific seasonal phenology"
  )
  
  # Assemble
  p <- (
    (p_a + p_b) /
      p_c
  ) +
    patchwork::plot_layout(
      heights = c(1, 1.15)
    ) +
    patchwork::plot_annotation(
      tag_levels = "a",
      theme = theme(
        plot.tag = element_text(size = 12)
      )
    )
  
  if (!is.null(fig_dir)) {
    
    dir.create(
      fig_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    ggsave(
      filename = file.path(fig_dir, filename),
      plot = p,
      width = 16,
      height = 8.5,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  }
  
  invisible(p)
}

# -----------------------------------------------------------------------------
# Environmental response plots from saved posterior predictions
# -----------------------------------------------------------------------------

make_env_posterior_plot <- function(
    posterior_draws,
    focal,
    raw_var,
    xlab,
    title = NULL,
    ylab = activity_ylab,
    y_limits = NULL,
    .width = c(0.50, 0.80, 0.89)
) {
  
  plot_data <- posterior_draws |>
    ungroup() |>
    group_by(
      .draw,
      across(all_of(c(focal, raw_var)))
    ) |>
    summarise(
      .epred = mean(.epred),
      .groups = "drop"
    )
  
  p <- ggplot(
    plot_data,
    aes(
      x = .data[[raw_var]],
      y = .epred
    )
  ) +
    stat_lineribbon(
      .width = .width,
      aes(fill = after_stat(level)),
      color = "gray15",
      linewidth = 0.9
    ) +
    scale_fill_brewer(
      palette = "Greys",
      direction = -1
    ) +
    labs(
      x = xlab,
      y = ylab,
      title = title
    ) +
    bat_theme() +
    theme(
      legend.position = "none"
    )
  
  if (!is.null(y_limits)) {
    p <- p +
      coord_cartesian(
        ylim = y_limits
      )
  }
  
  p
}

# -----------------------------------------------------------------------------
# Figure 4: Conditional temperature responses by frequency group
# -----------------------------------------------------------------------------

make_temperature_figure <- function(
    post_low,
    post_mid,
    post_high,
    fig_dir = NULL,
    filename = "Figure4_temperature.PNG"
) {
  
  p_low <- make_env_posterior_plot(
    posterior_draws = post_low,
    focal  = "scale_avetemp",
    raw_var = "avetemp",
    xlab   = expression("Nightly mean temperature (" * degree * "C)"),
    title  = "Low-frequency group"
  )
  
  p_mid <- make_env_posterior_plot(
    posterior_draws = post_mid,
    focal  = "scale_avetemp",
    raw_var = "avetemp",
    xlab   = expression("Nightly mean temperature (" * degree * "C)"),
    title  = "Mid-frequency group",
    ylab   = NULL
  )
  
  p_high <- make_env_posterior_plot(
    posterior_draws = post_high,
    focal  = "scale_avetemp",
    raw_var = "avetemp",
    xlab   = expression("Nightly mean temperature (" * degree * "C)"),
    title  = "High-frequency group",
    ylab   = NULL
  )
  
  p <- (p_low + p_mid + p_high) +
    plot_annotation(
      tag_levels = "a"
    )
  
  if (!is.null(fig_dir)) {
    ggsave(
      filename = file.path(fig_dir, filename),
      plot = p,
      width = 16,
      height = 5.5,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  }
  
  invisible(p)
}

make_rh_figure <- function(
    post_low,
    post_mid,
    post_high,
    fig_dir = NULL,
    filename = "Figure5_humidity.PNG"
) {
  
  p_low <- make_env_posterior_plot(
    posterior_draws = post_low,
    focal   = "scale_averh",
    raw_var = "averh",
    xlab    = "Nightly mean relative humidity (%)",
    title   = "Low-frequency group"
  )
  
  p_mid <- make_env_posterior_plot(
    posterior_draws = post_mid,
    focal   = "scale_averh",
    raw_var = "averh",
    xlab    = "Nightly mean relative humidity (%)",
    title   = "Mid-frequency group",
    ylab    = NULL
  )
  
  p_high <- make_env_posterior_plot(
    posterior_draws = post_high,
    focal   = "scale_averh",
    raw_var = "averh",
    xlab    = "Nightly mean relative humidity (%)",
    title   = "High-frequency group",
    ylab    = NULL
  )
  
  p <- (p_low + p_mid + p_high) +
    plot_annotation(
      tag_levels = "a"
    )
  
  if (!is.null(fig_dir)) {
    ggsave(
      filename = file.path(fig_dir, filename),
      plot = p,
      width = 16,
      height = 5.5,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  }
  
  invisible(p)
}

make_lunar_figure <- function(
    post_low,
    post_mid,
    post_high,
    fig_dir = NULL,
    filename = "Figure6_lunar.PNG"
) {
  
  p_low <- make_env_posterior_plot(
    posterior_draws = post_low,
    focal   = "Lunar",
    raw_var = "Lunar",
    xlab    = "Lunar illumination",
    title   = "Low-frequency group"
  )
  
  p_mid <- make_env_posterior_plot(
    posterior_draws = post_mid,
    focal   = "Lunar",
    raw_var = "Lunar",
    xlab    = "Lunar illumination",
    title   = "Mid-frequency group",
    ylab    = NULL
  )
  
  p_high <- make_env_posterior_plot(
    posterior_draws = post_high,
    focal   = "Lunar",
    raw_var = "Lunar",
    xlab    = "Lunar illumination",
    title   = "High-frequency group",
    ylab    = NULL
  )
  
  p <- (p_low + p_mid + p_high) +
    plot_annotation(tag_levels = "a")
  
  if (!is.null(fig_dir)) {
    ggsave(
      filename = file.path(fig_dir, filename),
      plot = p,
      width = 16,
      height = 5.5,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  }
  
  invisible(p)
}

# -----------------------------------------------------------------------------
# 5. Posterior predictive checks for the final model
# -----------------------------------------------------------------------------

# Compact diagnostic figure corresponding to the checks used during model
# development: distribution, mean, SD, and zero proportion.
make_final_ppc_figure <- function(model, fig_dir = NULL,
                                  filename = "aggregate_final_ppc.PNG",
                                  ndraws = 100) {

  p_density <- pp_check(model, ndraws = ndraws, type = "dens_overlay") +
    ggtitle("Distribution")

  p_mean <- pp_check(model, ndraws = ndraws, type = "stat", stat = "mean") +
    ggtitle("Mean")

  p_sd <- pp_check(model, ndraws = ndraws, type = "stat", stat = "sd") +
    ggtitle("Standard deviation")

  prop_zero <- function(y) mean(y == 0)
  p_zero <- pp_check(model, ndraws = ndraws, type = "stat", stat = prop_zero) +
    ggtitle("Proportion zero")

  p <- (p_density + p_mean) / (p_sd + p_zero) +
    plot_annotation(tag_levels = "a") &
    theme_minimal(base_size = 11)

  if (!is.null(fig_dir)) {
    ggsave(
      filename = file.path(fig_dir, filename),
      plot = p, width = 13, height = 10, units = "in", dpi = 300
    )
  }

  invisible(p)
}


# Site-level posterior predictive density overlays. This is supplementary: it
# helps reveal whether pooled PPC discrepancies are concentrated at particular
# sites. Uses observed rows to define site panels.
make_site_ppc_figure <- function(model, data, fig_dir = NULL,
                                 filename = "aggregate_site_ppc.PNG",
                                 ndraws = 50) {

  site_levels <- levels(factor(data$Site))

  panels <- map(site_levels, function(s) {
    rows <- which(as.character(data$Site) == s)

    pp_check(
      model,
      ndraws = ndraws,
      type = "dens_overlay",
      newdata = data[rows, , drop = FALSE]
    ) +
      ggtitle(s) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none")
  })

  p <- wrap_plots(panels, ncol = 4) +
    plot_annotation(tag_levels = "a")

  if (!is.null(fig_dir)) {
    ggsave(
      filename = file.path(fig_dir, filename),
      plot = p, width = 16, height = 9, units = "in", dpi = 300
    )
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# 6. Residual calendar-lag ACF figure
# -----------------------------------------------------------------------------

# Plot output from calendar_acf() in the refactored test.R.
make_calendar_acf_plot <- function(acf_df, fig_dir = NULL,
                                   filename = "aggregate_calendar_acf.PNG") {

  p <- ggplot(acf_df, aes(x = lag, y = correlation)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line() +
    geom_point(size = 2) +
    scale_x_continuous(breaks = acf_df$lag) +
    labs(
      x = "Calendar lag (days)",
      y = "Residual correlation"
    ) +
    theme_minimal(base_size = 12)

  if (!is.null(fig_dir)) {
    ggsave(
      filename = file.path(fig_dir, filename),
      plot = p, width = 8, height = 5.5, units = "in", dpi = 300
    )
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# 7. Convenience wrapper for final aggregate figures
# -----------------------------------------------------------------------------

# Generate the figures that are currently justified by the final aggregate
# analysis. `season_model` can be NULL until the final total-season model has
# been fitted and accepted.
make_final_aggregate_figures <- function(env_model, data, fig_dir,
                                         season_model = NULL,
                                         acf_df = NULL) {

  env_plot <- make_aggregate_environment_figure(
    env_model, data, fig_dir = fig_dir
  )

  ppc_plot <- make_final_ppc_figure(
    env_model, fig_dir = fig_dir
  )

  site_ppc_plot <- make_site_ppc_figure(
    env_model, data, fig_dir = fig_dir
  )

  if (!is.null(acf_df)) {
    acf_plot <- make_calendar_acf_plot(
      acf_df, fig_dir = fig_dir
    )
  } else {
    acf_plot <- NULL
  }

  if (!is.null(season_model)) {
    season_plot <- make_season_figure(
      season_model, data, fig_dir = fig_dir
    )

    season_compare <- make_season_model_comparison(
      season_model = season_model,
      env_model    = env_model,
      data         = data,
      fig_dir      = fig_dir
    )
  } else {
    season_plot <- NULL
    season_compare <- NULL
  }

  invisible(list(
    environment       = env_plot,
    season            = season_plot,
    season_comparison = season_compare,
    ppc               = ppc_plot,
    site_ppc          = site_ppc_plot,
    calendar_acf      = acf_plot
  ))
}


# =============================================================================
# NOTE ON FREQUENCY-GROUP FIGURES
# =============================================================================
# The old script automatically generated low/mid/high-frequency figures from
# the old Poisson model set. Those functions are intentionally not carried into
# this production file yet. Frequency-group figures should be restored only
# after the diagnostic fits establish the appropriate likelihood/seasonal
# structure for each group. This prevents the plotting pipeline from silently
# treating unvalidated models as final results.
# =============================================================================


# =============================================================================
# Supplementary figures
# =============================================================================

# -----------------------------------------------------------------------------
# Figure S1: Site-specific seasonal phenology
# -----------------------------------------------------------------------------

make_figure_s1 <- function(post_sites, post_average, fig_dir = NULL,
                           filename = "Figure_S1_site_seasonality.PNG") {

  p <- plot_site_season_posterior(
    post_sites = post_sites,
    post_average = post_average,
    title = NULL
  )

  if (!is.null(fig_dir)) {
    ggsave(
      file.path(fig_dir, filename),
      p, width = 11, height = 6.5, units = "in", dpi = 300, bg = "white"
    )
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# Calendar-aware residual correlation
# -----------------------------------------------------------------------------

calendar_acf <- function(model, data, max_lag = 14) {

  r <- residuals(model, summary = TRUE)[, "Estimate"]

  dat <- data |>
    transmute(
      Site = as.character(Site),
      date = as.Date(NiteDayDate),
      resid = r
    )

  map_dfr(seq_len(max_lag), function(k) {

    paired <- dat |>
      transmute(
        Site,
        date2 = date + k,
        resid_lag = resid
      ) |>
      inner_join(
        dat |>
          transmute(Site, date2 = date, resid_now = resid),
        by = c("Site", "date2")
      )

    tibble(
      lag = k,
      n_pairs = nrow(paired),
      correlation = if (nrow(paired) > 2) {
        cor(paired$resid_lag, paired$resid_now, use = "complete.obs")
      } else {
        NA_real_
      }
    )
  })
}


# -----------------------------------------------------------------------------
# Figure S2: Calendar-aware residual temporal dependence
# -----------------------------------------------------------------------------

make_figure_s2 <- function(models, data_list, fig_dir = NULL,
                           filename = "Figure_S2_residual_acf.PNG",
                           max_lag = 14) {

  acf_data <- imap_dfr(
    models,
    ~ calendar_acf(.x, data_list[[.y]], max_lag = max_lag) |>
      mutate(Model = .y)
  ) |>
    mutate(
      Model = recode(
        Model,
        Aggregate = "Aggregate",
        Low = "Low-frequency group",
        Mid = "Mid-frequency group",
        High = "High-frequency group"
      )
    )

  p <- ggplot(acf_data, aes(lag, correlation)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.6) +
    facet_wrap(~ Model, ncol = 2) +
    scale_x_continuous(breaks = seq(2, max_lag, by = 2)) +
    labs(
      x = "Calendar lag (days)",
      y = "Residual correlation"
    ) +
    theme_minimal(base_size = 11)

  if (!is.null(fig_dir)) {
    ggsave(
      file.path(fig_dir, filename),
      p, width = 10, height = 7.5, units = "in", dpi = 300, bg = "white"
    )
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# Figure S3: Aggregate posterior predictive checks
# -----------------------------------------------------------------------------

make_figure_s3 <- function(model, fig_dir = NULL,
                           filename = "Figure_S3_aggregate_ppc.PNG",
                           ndraws = 100) {

  make_final_ppc_figure(
    model = model,
    fig_dir = fig_dir,
    filename = filename,
    ndraws = ndraws
  )
}


# -----------------------------------------------------------------------------
# Figure S4: Site-specific posterior predictive checks
# -----------------------------------------------------------------------------

make_figure_s4 <- function(model, data, fig_dir = NULL,
                           filename = "Figure_S4_site_ppc.PNG",
                           ndraws = 500) {

  yrep <- posterior_predict(model, ndraws = ndraws)

  site_levels <- unique(as.character(data$Site))

  site_ppc <- map_dfr(site_levels, function(s) {

    idx <- which(as.character(data$Site) == s)
    obs <- data$n[idx]

    rep_mean <- rowMeans(yrep[, idx, drop = FALSE])
    rep_sd <- apply(yrep[, idx, drop = FALSE], 1, sd)

    bind_rows(
      tibble(
        Site = s,
        Statistic = "Mean",
        Observed = mean(obs),
        Median = median(rep_mean),
        Lower = quantile(rep_mean, 0.055),
        Upper = quantile(rep_mean, 0.945)
      ),
      tibble(
        Site = s,
        Statistic = "Standard deviation",
        Observed = sd(obs),
        Median = median(rep_sd),
        Lower = quantile(rep_sd, 0.055),
        Upper = quantile(rep_sd, 0.945)
      )
    )
  })

  p <- ggplot(site_ppc, aes(x = Site)) +
    geom_linerange(
      aes(ymin = Lower, ymax = Upper),
      linewidth = 0.7
    ) +
    geom_point(aes(y = Median), shape = 21, fill = "white", size = 2.5) +
    geom_point(aes(y = Observed), shape = 4, size = 2.8, stroke = 0.9) +
    facet_wrap(~ Statistic, scales = "free_y", ncol = 1) +
    labs(
      x = NULL,
      y = "Intervals per night"
    ) +
    theme_minimal(base_size = 11)

  if (!is.null(fig_dir)) {
    ggsave(
      file.path(fig_dir, filename),
      p, width = 10, height = 7.5, units = "in", dpi = 300, bg = "white"
    )
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# Figure S5: Frequency-group posterior predictive checks
# -----------------------------------------------------------------------------

make_figure_s5 <- function(models, fig_dir = NULL,
                           filename = "Figure_S5_frequency_ppc.PNG",
                           ndraws = 100) {

  group_titles <- c(
    Low = "Low-frequency group",
    Mid = "Mid-frequency group",
    High = "High-frequency group"
  )

  make_row <- function(model, title) {

    p1 <- pp_check(model, ndraws = ndraws, type = "dens_overlay") +
      labs(title = title, subtitle = "Distribution")

    p2 <- pp_check(model, ndraws = ndraws, type = "stat", stat = "sd") +
      labs(title = NULL, subtitle = "Standard deviation")

    prop_zero <- function(y) mean(y == 0)
    p3 <- pp_check(model, ndraws = ndraws, type = "stat", stat = prop_zero) +
      labs(title = NULL, subtitle = "Proportion zero")

    p1 + p2 + p3
  }

  p <- (
    make_row(models$Low, group_titles["Low"]) /
    make_row(models$Mid, group_titles["Mid"]) /
    make_row(models$High, group_titles["High"])
  ) +
    plot_annotation(tag_levels = "a") &
    theme_minimal(base_size = 10)

  if (!is.null(fig_dir)) {
    ggsave(
      file.path(fig_dir, filename),
      p, width = 14, height = 12, units = "in", dpi = 300, bg = "white"
    )
  }

  invisible(p)
}
