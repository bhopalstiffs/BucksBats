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
# Table S1: Model specifications
# -----------------------------------------------------------------------------

make_table_s1 <- function(models_lookup, table_dir) {
  s1_data <- bind_rows(
    build_spec_row(models_lookup$agg$season_total,    "Aggregate: Season-only",       "Community total",      0.99),
    build_spec_row(models_lookup$agg$temp_total,      "Aggregate: Temperature total", "Community total",      0.95),
    build_spec_row(models_lookup$agg$env_conditional, "Aggregate: Full conditional",  "Community total",      0.95),
    build_spec_row(models_lookup$low$season_total,    "Low: Season-only",             "Low-frequency group",  0.95),
    build_spec_row(models_lookup$low$temp_total,      "Low: Temperature total",       "Low-frequency group",  0.95),
    build_spec_row(models_lookup$low$env_conditional, "Low: Full conditional",        "Low-frequency group",  0.95),
    build_spec_row(models_lookup$mid$season_total,    "Mid: Season-only",             "Mid-frequency group",  0.95),
    build_spec_row(models_lookup$mid$temp_total,      "Mid: Temperature total",       "Mid-frequency group",  0.95),
    build_spec_row(models_lookup$mid$env_conditional, "Mid: Full conditional",        "Mid-frequency group",  0.95),
    build_spec_row(models_lookup$high$season_total,   "High: Season-only",            "High-frequency group", 0.95),
    build_spec_row(models_lookup$high$temp_total,     "High: Temperature total",      "High-frequency group", 0.95),
    build_spec_row(models_lookup$high$env_conditional,"High: Full conditional",       "High-frequency group", 0.95)
  )

  caption_s1 <- paste(
    "Table S1. Specifications for all Bayesian generalized additive models fit in this study.",
    "All models used a Poisson likelihood with log link, 4 MCMC chains \u00d7 4,000 iterations",
    "(500 warmup discarded), and a normal(0, 10) prior on the intercept. Aggregate models",
    "were fit on community-total activity (intervals summed across all species) at each",
    "site-night. Frequency-stratified models were fit on group-level activity (intervals",
    "summed across species within each characteristic-frequency class). Site identity was",
    "included as a random intercept in all models. The aggregate season-only model used",
    "adapt_delta = 0.99 to eliminate divergent transitions arising from weak identification",
    "of the smoothing-spline hyperparameter; all other models used adapt_delta = 0.95."
  )

  table_s1 <- flextable(s1_data) |>
    bold(part = "header") |>
    align(j = c("N", "adapt_delta", "Post-warmup draws"),
          align = "center", part = "all") |>
    fontsize(j = "Formula", size = 9, part = "body") |>
    merge_v(j = "Dataset") |>
    valign(j = "Dataset", valign = "top") |>
    set_table_properties(layout = "autofit") |>
    set_caption(
      caption          = caption_s1,
      fp_p             = fp_par(text.align = "left", padding = 3),
      align_with_table = FALSE
    )

  save_as_docx(table_s1, path = file.path(table_dir, "Table_S1.docx"))
  invisible(table_s1)
}


# -----------------------------------------------------------------------------
# Table S2: Aggregate posterior summaries
# -----------------------------------------------------------------------------

make_table_s2 <- function(models_aggregate, table_dir) {
  s2_data <- bind_rows(
    extract_model_summary(models_aggregate$season_total,    "Season-only (total seasonal effect)"),
    extract_model_summary(models_aggregate$env_conditional, "Full conditional (direct effects)")
  )

  caption_s2 <- paste(
    "Table S2. Posterior summaries for aggregate models of community-level bat activity",
    "(N = 1,231 site-nights, 8 sites). Both models used a Poisson likelihood with site",
    "random intercepts. Posterior median, 89% credible interval, R-hat convergence",
    "diagnostic, and bulk effective sample size are shown for each parameter."
  )

  table_s2 <- build_summary_table(
    s2_data, group_col = "Model", caption_text = caption_s2
  )
  save_as_docx(table_s2, path = file.path(table_dir, "Table_S2.docx"))
  invisible(table_s2)
}


# -----------------------------------------------------------------------------
# Table S3a/b/c: Per-group posterior summaries
# -----------------------------------------------------------------------------

make_table_s3 <- function(models_lookup, table_dir) {

  # --- helper to assemble one S3 sub-table ----------------------------------
  assemble_group_summary <- function(model_field, caption_text, file_name) {
    dat <- bind_rows(
      extract_model_summary(models_lookup$low[[model_field]],  "Low"),
      extract_model_summary(models_lookup$mid[[model_field]],  "Mid"),
      extract_model_summary(models_lookup$high[[model_field]], "High")
    ) |>
      rename(Group = Model) |>
      mutate(Group = factor(Group, levels = c("Low", "Mid", "High")))

    tab <- build_summary_table(dat, group_col = "Group",
                               caption_text = caption_text)
    save_as_docx(tab, path = file.path(table_dir, file_name))
    tab
  }

  # --- S3a: env_conditional -------------------------------------------------
  caption_s3a <- paste(
    "Table S3a. Posterior summaries for the full conditional model",
    "(n ~ s(doy) + s(scale_avetemp) + s(scale_averh) + s(Lunar) + (1|Site))",
    "fit separately by frequency group (N = 1,231 site-nights per group, 8 sites).",
    "This model estimates the direct effect of each environmental covariate after",
    "conditioning on the others. Posterior median, 89% credible interval, R-hat,",
    "and bulk effective sample size shown."
  )

  # --- S3b: temp_total ------------------------------------------------------
  caption_s3b <- paste(
    "Table S3b. Posterior summaries for the temperature-total model",
    "(n ~ s(doy) + s(scale_avetemp) + s(Lunar) + (1|Site)), fit separately by",
    "frequency group. Humidity is omitted from this specification to estimate",
    "the total effect of temperature on bat activity, including any indirect",
    "effect mediated through humidity. Comparison with Table S3a quantifies",
    "the extent of humidity mediation in each group."
  )

  # --- S3c: season_total ----------------------------------------------------
  caption_s3c <- paste(
    "Table S3c. Posterior summaries for the season-only model",
    "(n ~ s(doy) + (1|Site)), fit separately by frequency group. This",
    "specification estimates the total seasonal effect of day of year on bat",
    "activity prior to weather adjustment. Comparison with the doy smooth",
    "coefficients in Table S3a quantifies the proportion of seasonal variation",
    "in each group that is mediated by nightly weather."
  )

  s3a <- assemble_group_summary("env_conditional", caption_s3a, "Table_S3a.docx")
  s3b <- assemble_group_summary("temp_total",      caption_s3b, "Table_S3b.docx")
  s3c <- assemble_group_summary("season_total",    caption_s3c, "Table_S3c.docx")

  invisible(list(s3a = s3a, s3b = s3b, s3c = s3c))
}


# -----------------------------------------------------------------------------
# Table S4: Species composition
# -----------------------------------------------------------------------------

make_table_s4 <- function(cdata_species, table_dir) {

  species_lookup <- tribble(
    ~Species, ~scientific_name,           ~common_name,                 ~freq_group, ~char_freq_kHz,
    "MYOYUM", "Myotis yumanensis",        "Yuma myotis",                "high", 48,
    "MYOLUC", "Myotis lucifugus",         "Little brown myotis",        "high", 42,
    "MYOCAL", "Myotis californicus",      "California myotis",          "high", 50,
    "MYOCIL", "Myotis ciliolabrum",       "Western small-footed myotis","high", 45,
    "MYOVOL", "Myotis volans",            "Long-legged myotis",         "high", 40,
    "PARHES", "Parastrellus hesperus",    "Canyon bat",                 "high", 45,
    "LASBLO", "Lasiurus blossevillii",    "Western red bat",            "high", 42,
    "MYOEVO", "Myotis evotis",            "Long-eared myotis",          "mid",  35,
    "MYOTHY", "Myotis thysanodes",        "Fringed myotis",             "mid",  32,
    "EPTFUS", "Eptesicus fuscus",         "Big brown bat",              "mid",  27,
    "LASNOC", "Lasionycteris noctivagans","Silver-haired bat",          "mid",  26,
    "TADBRA", "Tadarida brasiliensis",    "Mexican free-tailed bat",    "mid",  25,
    "ANTPAL", "Antrozous pallidus",       "Pallid bat",                 "mid",  30,
    "CORTOW", "Corynorhinus townsendii",  "Townsend's big-eared bat",   "mid",  28,
    "LASCIN", "Lasiurus cinereus",        "Hoary bat",                  "low",  20,
    "EUMPER", "Eumops perotis",           "Western mastiff bat",        "low",  12,
    "EUDMAC", "Euderma maculatum",        "Spotted bat",                "low",  11
  )

  species_table <- cdata_species |>
    group_by(Species) |>
    summarise(total_intervals     = sum(n),
              site_nights_present = n(),
              .groups = "drop") |>
    left_join(species_lookup, by = "Species") |>
    mutate(freq_group = factor(freq_group, levels = c("low", "mid", "high"))) |>
    arrange(freq_group, desc(total_intervals)) |>
    select(
      `Species code`                   = Species,
      `Scientific name`                = scientific_name,
      `Common name`                    = common_name,
      `Characteristic frequency (kHz)` = char_freq_kHz,
      `Frequency group`                = freq_group,
      `Total intervals`                = total_intervals,
      `Site-nights detected`           = site_nights_present
    )

  caption_s4 <- paste(
    "Table S4. Species composition of the acoustic dataset, classified by characteristic",
    "echolocation frequency group. Characteristic frequency (Fc) was drawn from species",
    "accounts in the Western Bat Working Group reference and the Kaleidoscope Pro classifier",
    "reference call library (Wildlife Acoustics, Maynard, MA, USA); group boundaries",
    "(low <25 kHz; mid 25\u201340 kHz; high >40 kHz) reflect the frequency ranges over",
    "which atmospheric absorption coefficients of ultrasound differ qualitatively in",
    "their humidity sensitivity (Lawrence & Simmons 1982). Species with plastic call",
    "frequencies (notably Tadarida brasiliensis) were classified based on their typical",
    "search-phase Fc."
  )

  species_ft <- flextable(species_table) |>
    italic(j = "Scientific name") |>
    bold(part = "header") |>
    align(j = c("Characteristic frequency (kHz)",
                "Total intervals",
                "Site-nights detected"),
          align = "center", part = "all") |>
    set_table_properties(layout = "autofit") |>
    set_caption(
      caption          = caption_s4,
      fp_p             = fp_par(text.align = "left", padding = 3),
      align_with_table = FALSE
    )

  save_as_docx(species_ft, path = file.path(table_dir, "Table_S4.docx"))
  invisible(species_ft)
}
