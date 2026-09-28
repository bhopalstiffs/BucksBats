calendar_acf <- function(model, data, max_lag = 14) {
  
  resid_df <- data %>%
    mutate(
      resid = residuals(
        model,
        summary = TRUE
      )[, "Estimate"]
    ) %>%
    select(
      Site,
      NiteDayDate,
      resid
    )
  
  map_dfr(seq_len(max_lag), function(k) {
    
    pairs <- resid_df %>%
      select(
        Site,
        date1 = NiteDayDate,
        resid1 = resid
      ) %>%
      mutate(
        date2 = date1 + days(k)
      ) %>%
      inner_join(
        resid_df %>%
          select(
            Site,
            date2 = NiteDayDate,
            resid2 = resid
          ),
        by = c("Site", "date2")
      )
    
    tibble(
      lag = k,
      n_pairs = nrow(pairs),
      correlation = if (nrow(pairs) > 2) {
        cor(
          pairs$resid1,
          pairs$resid2,
          use = "complete.obs"
        )
      } else {
        NA_real_
      }
    )
  })
}

# -----------------------------------------------------------------------------
# Save standard posterior predictive checks
# -----------------------------------------------------------------------------

save_pp_checks <- function(
    model,
    model_name,
    output_dir = "figures/diagnostics",
    ndraws = 100,
    width = 7,
    height = 5,
    dpi = 300
) {
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Custom statistic: number of zero observations
  n_zeros <- function(x) {
    sum(x == 0)
  }
  
  # ---------------------------------------------------------------------------
  # Posterior predictive checks
  # ---------------------------------------------------------------------------
  
  p_density <- brms::pp_check(
    model,
    type = "dens_overlay",
    ndraws = ndraws
  )
  
  p_sd <- brms::pp_check(
    model,
    type = "stat",
    stat = sd,
    ndraws = ndraws
  )
  
  p_zeros <- brms::pp_check(
    model,
    type = "stat",
    stat = n_zeros,
    ndraws = ndraws
  )
  
  # ---------------------------------------------------------------------------
  # Save figures
  # ---------------------------------------------------------------------------
  
  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(model_name, "_ppc_density.png")
    ),
    plot = p_density,
    width = width,
    height = height,
    dpi = dpi
  )
  
  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(model_name, "_ppc_sd.png")
    ),
    plot = p_sd,
    width = width,
    height = height,
    dpi = dpi
  )
  
  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(model_name, "_ppc_zeros.png")
    ),
    plot = p_zeros,
    width = width,
    height = height,
    dpi = dpi
  )
  
  # ---------------------------------------------------------------------------
  # Return plots for inspection if desired
  # ---------------------------------------------------------------------------
  
  invisible(
    list(
      density = p_density,
      sd = p_sd,
      zeros = p_zeros
    )
  )
}
