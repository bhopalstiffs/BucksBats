# =============================================================================
# Figure 2: Directed acyclic graph
# =============================================================================
#
# Purpose:
#   Define the causal structure used to distinguish seasonal and nightly
#   environmental associations with bat acoustic activity and generate the
#   manuscript DAG.
#
# Output:
#   figs/Figure2_DAG.PNG
#
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------

library(tidyverse)
library(dagitty)
library(ggdag)
library(ggforce)


# -----------------------------------------------------------------------------
# 2. Output directory
# -----------------------------------------------------------------------------

fig_dir <- "figs"

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------------------------------------------------------
# 3. Define causal structure
# -----------------------------------------------------------------------------

bat_dag <- dagify(
  activity    ~ doy + temperature + humidity + lunar,
  temperature ~ doy,
  humidity    ~ doy + temperature,
  
  outcome = "activity",
  
  labels = c(
    doy         = "Day of year",
    temperature = "Temperature",
    humidity    = "Relative humidity",
    lunar       = "Lunar phase",
    activity    = "Bat activity"
  ),
  
  coords = list(
    x = c(
      doy         =  2.2,
      temperature =  1.0,
      humidity    =  0.0,
      lunar       = -1.6,
      activity    =  0.0
    ),
    y = c(
      doy         =  0.0,
      temperature = -1.25,
      humidity    =  1.25,
      lunar       =  0.25,
      activity    =  0.0
    )
  )
)

# -----------------------------------------------------------------------------
# 6. Manuscript DAG plotting data
# -----------------------------------------------------------------------------

dag_nodes <- tibble(
  name = c(
    "lunar",
    "humidity",
    "activity",
    "temperature",
    "doy"
  ),
  label = c(
    "Lunar phase",
    "Relative humidity",
    "Bat activity",
    "Temperature",
    "Day of year"
  ),
  x = c(
    -2.00,
    0.0,
    -0.10,
    1.0,
    2.3
  ),
  y = c(
    0.30,
    1.30,
    0.00,
    -1.30,
    0.00
  )
)


# -----------------------------------------------------------------------------
# 7. Arrow coordinates
#
# Coordinates terminate near the ellipse boundaries rather than at node
# centers so that arrowheads remain visible.
# -----------------------------------------------------------------------------

dag_edges <- tribble(
  ~from,          ~to,
  "lunar",        "activity",
  "doy",          "activity",
  "doy",          "temperature",
  "doy",          "humidity",
  "temperature",  "activity",
  "temperature",  "humidity",
  "humidity",     "activity"
)


# -----------------------------------------------------------------------------
# 8. Calculate arrow endpoints
# -----------------------------------------------------------------------------

# Ellipse dimensions used for every node
ellipse_a <- 0.58
ellipse_b <- 0.20

edge_plot <- dag_edges |>
  left_join(
    dag_nodes |> select(from = name, x_from = x, y_from = y),
    by = "from"
  ) |>
  left_join(
    dag_nodes |> select(to = name, x_to = x, y_to = y),
    by = "to"
  ) |>
  rowwise() |>
  mutate(
    dx = x_to - x_from,
    dy = y_to - y_from,
    
    # Distance from ellipse center to its boundary in the direction
    # of the connecting edge.
    start_r = 1 / sqrt(
      (dx / ellipse_a)^2 +
        (dy / ellipse_b)^2
    ),
    
    end_r = start_r,
    
    x = x_from + start_r * dx,
    y = y_from + start_r * dy,
    
    xend = x_to - end_r * dx,
    yend = y_to - end_r * dy
  ) |>
  ungroup()

# -----------------------------------------------------------------------------
# 9. Plot Figure 2
# -----------------------------------------------------------------------------

fig2_dag <- ggplot() +
  
  # Causal arrows
  geom_segment(
    data = edge_plot,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    linewidth = 0.7,
    color = "gray20",
    arrow = arrow(
      type = "closed",
      length = unit(0.11, "in")
    )
  ) +
  
  # Nodes
  ggforce::geom_ellipse(
    data = dag_nodes,
    aes(
      x0 = x,
      y0 = y,
      a = ellipse_a,
      b = ellipse_b,
      angle = 0
    ),
    inherit.aes = FALSE,
    fill = "white",
    color = "gray30",
    linewidth = 0.8
  ) +
  
  # Node labels
  geom_text(
    data = dag_nodes,
    aes(
      x = x,
      y = y,
      label = label
    ),
    size = 4.2,
    color = "gray10"
  ) +
  
  coord_equal(
    clip = "off"
  ) +
  
  theme_void() +
  
  theme(
    plot.margin = margin(
      t = 20,
      r = 30,
      b = 20,
      l = 30
    )
  )


fig2_dag

# -----------------------------------------------------------------------------
# 8. Save Figure 2
# -----------------------------------------------------------------------------

ggsave(
  filename = file.path(
    fig_dir,
    "Figure2_DAG.PNG"
  ),
  plot = fig2_dag,
  width = 9,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)


# -----------------------------------------------------------------------------
# 9. Check implied adjustment sets
# -----------------------------------------------------------------------------

adjustmentSets(
  bat_dag,
  exposure = "temperature",
  outcome = "activity"
)

adjustmentSets(
  bat_dag,
  exposure = "humidity",
  outcome = "activity"
)

adjustmentSets(
  bat_dag,
  exposure = "lunar",
  outcome = "activity"
)

