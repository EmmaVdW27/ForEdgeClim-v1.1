############################################################################################
# ForEdgeClim - Thermal microhabitat (= niche) analysis visualisation
#
# MAIN FIGURE:
#   Daytime thermal niche characteristics
#
#   (a) Thermal niche availability
#   (b) Vertical distribution of thermal niches
#   (c) Spatial connectivity of thermal niches
#
# APPENDIX FIGURE:
#   Same analyses for nighttime conditions
#
# Vertical niche distribution:
#   point    = median height
#   interval = P10-P90 height range (central 80% of niche voxels)
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################


# ==========================================================================================
# PACKAGES
# ==========================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(cowplot)
  library(ggtext)
})


# ==========================================================================================
# INPUT
# ==========================================================================================

input_path <- "Output/thermal_niches/data/"

output_path <- "Output/thermal_niches/plots/"

dir.create(
  output_path,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==========================================================================================
# LOAD DATA
# ==========================================================================================

control <- readRDS(
  file.path(
    input_path,
    "thermal_niche_metrics_control_calibrated.rds"
  )
)

removal <- readRDS(
  file.path(
    input_path,
    "thermal_niche_metrics_removal_calibrated.rds"
  )
)


# ==========================================================================================
# ADD TREATMENT
# ==========================================================================================

control <- control %>%
  mutate(
    treatment = "Control"
  )

removal <- removal %>%
  mutate(
    treatment = "Liana removal"
  )


# ==========================================================================================
# COMBINE DATA
# ==========================================================================================

dat <- bind_rows(
  control,
  removal
)


# ==========================================================================================
# FACTOR LEVELS
# ==========================================================================================

niche_levels <- c(
  "Very cool",
  "Cool",
  "Near-average",
  "Warm",
  "Very warm"
)

season_levels <- c(
  "End dry season",
  "End wet season"
)

period_levels <- c(
  "Day",
  "Night"
)


dat <- dat %>%
  mutate(

    treatment = factor(
      treatment,
      levels = c(
        "Control",
        "Liana removal"
      )
    ),

    thermal_niche = factor(
      thermal_niche,
      levels = niche_levels
    ),

    seasonal_period = factor(
      seasonal_period,
      levels = season_levels
    ),

    period = factor(
      period,
      levels = period_levels
    )
  )


# ==========================================================================================
# COLOURS
# ==========================================================================================

# Treatment colours
treatment_colours <- c(
  "Control"       = "chocolate4",
  "Liana removal" = "darkgreen"
)

# Niche colours
# Same colour scheme as the 3D thermal niche visualisation
niche_colours <- c(
  "Very cool"    = "blue",
  "Cool"         = "lightblue",
  "Near-average" = "bisque2",
  "Warm"         = "orange",
  "Very warm"    = "red"
)


# ==========================================================================================
# COLOURED NICHE LABELS
#
# Text itself is coloured rather than trying to draw coloured backgrounds.
# ==========================================================================================

niche_labels_coloured <- c(

  "Very cool" =
    "<span style='color:blue'><b>Very cool</b></span>",

  "Cool" =
    "<span style='color:#ADD8E6'><b>Cool</b></span>",

  "Near-average" =
    "<span style='color:#EED5B7'><b>Near-average</b></span>",

  "Warm" =
    "<span style='color:orange'><b>Warm</b></span>",

  "Very warm" =
    "<span style='color:red'><b>Very warm</b></span>"
)


# ==========================================================================================
# POSITION DODGE
# ==========================================================================================

dodge_width <- 0.55


############################################################################################
############################################################################################
#
# STEP 1
#
# SUMMARISE THROUGH TIME WITHIN EACH INDIVIDUAL PLOT
#
############################################################################################
############################################################################################


# ==========================================================================================
# NICHE AVAILABILITY
#
# Percentage of the complete forest volume occupied by each thermal niche.
# ==========================================================================================

availability_plot <- dat %>%

  group_by(
    treatment,
    plot,
    seasonal_period,
    period,
    thermal_niche
  ) %>%

  summarise(

    niche_availability = mean(
      niche_availability,
      na.rm = TRUE
    ) * 100,

    .groups = "drop"
  )


# ==========================================================================================
# VERTICAL NICHE DISTRIBUTION
#
# Hourly median height, P10 height and P90 height are summarised within
# each individual plot.
# ==========================================================================================

vertical_plot <- dat %>%

  group_by(
    treatment,
    plot,
    seasonal_period,
    period,
    thermal_niche
  ) %>%

  summarise(

    median_height = median(
      median_height,
      na.rm = TRUE
    ),

    p10_height = median(
      p10_height,
      na.rm = TRUE
    ),

    p90_height = median(
      p90_height,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%

  filter(
    is.finite(median_height),
    is.finite(p10_height),
    is.finite(p90_height)
  )


# ==========================================================================================
# SPATIAL NICHE CONNECTIVITY
#
# Higher same_niche_neighbour_fraction values indicate that voxels belonging
# to a thermal niche more frequently share faces with voxels belonging to
# that same niche.
# ==========================================================================================

connectivity_plot <- dat %>%

  group_by(
    treatment,
    plot,
    seasonal_period,
    period,
    thermal_niche
  ) %>%

  summarise(

    connectivity = mean(
      same_niche_neighbour_fraction,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%

  filter(
    is.finite(connectivity)
  )


############################################################################################
############################################################################################
#
# STEP 2
#
# TREATMENT MEDIANS ACROSS THE FOUR PLOTS
#
############################################################################################
############################################################################################


# ==========================================================================================
# NICHE AVAILABILITY
# ==========================================================================================

availability_summary <- availability_plot %>%

  group_by(
    treatment,
    seasonal_period,
    period,
    thermal_niche
  ) %>%

  summarise(

    median_availability = median(
      niche_availability,
      na.rm = TRUE
    ),

    .groups = "drop"
  )


# ==========================================================================================
# VERTICAL NICHE DISTRIBUTION
# ==========================================================================================

vertical_summary <- vertical_plot %>%

  group_by(
    treatment,
    seasonal_period,
    period,
    thermal_niche
  ) %>%

  summarise(

    median_height = median(
      median_height,
      na.rm = TRUE
    ),

    p10_height = median(
      p10_height,
      na.rm = TRUE
    ),

    p90_height = median(
      p90_height,
      na.rm = TRUE
    ),

    .groups = "drop"
  )


# ==========================================================================================
# SPATIAL CONNECTIVITY
# ==========================================================================================

connectivity_summary <- connectivity_plot %>%

  group_by(
    treatment,
    seasonal_period,
    period,
    thermal_niche
  ) %>%

  summarise(

    median_connectivity = median(
      connectivity,
      na.rm = TRUE
    ),

    .groups = "drop"
  )


############################################################################################
############################################################################################
#
# COMMON THEME
#
############################################################################################
############################################################################################

common_theme <- theme_bw() +

  theme(

    plot.title = element_text(
      size = 30,
      face = "plain",
      margin = margin(
        b = 8
      )
    ),

    axis.title = element_text(
      size = 28
    ),

    axis.text.x = element_text(
      size = 24
    ),

    # Needed to render the coloured niche names
    axis.text.y = ggtext::element_markdown(
      size = 25
    ),

    axis.ticks.y = element_line(),

    strip.text = element_text(
      size = 29,
      face = "bold"
    ),

    # More space between end-dry and end-wet facets
    panel.spacing.x = unit(
      2,
      "cm"
    ),

    panel.grid.minor = element_blank(),

    # Legends are suppressed in the actual plots.
    # A separate treatment key is added above the complete figure.
    legend.position = "none",

    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )


############################################################################################
############################################################################################
#
# PANEL A
#
# THERMAL NICHE AVAILABILITY
#
############################################################################################
############################################################################################

create_availability_plot <- function(selected_period) {


  plot_individual <- availability_plot %>%

    filter(
      period == selected_period
    )


  plot_median <- availability_summary %>%

    filter(
      period == selected_period
    )


  ggplot(

    plot_individual,

    aes(
      x = niche_availability,
      y = thermal_niche,
      colour = treatment
    )

  ) +

    # ---- Individual plots ----

  geom_point(

    position = position_dodge(
      width = dodge_width
    ),

    size = 3,

    alpha = 0.25
  ) +

    # ---- Treatment median ----

  geom_point(

    data = plot_median,

    aes(
      x = median_availability,
      y = thermal_niche,
      colour = treatment
    ),

    position = position_dodge(
      width = dodge_width
    ),

    size = 5,

    inherit.aes = FALSE
  ) +

    # ---- Seasons ----

  facet_grid(
    . ~ seasonal_period
  ) +

    # ---- Treatment colours ----

  scale_colour_manual(
    values = treatment_colours,
    name = NULL
  ) +

    # ---- Niche labels ----

  scale_y_discrete(
    labels = niche_labels_coloured,
    drop = FALSE
  ) +

    # ---- X axis ----

  scale_x_continuous(
    limits = c(-5, 105),
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0)
  ) +

    # ---- Labels ----

  labs(

    x = "Microhabitat availability (% of forest volume)",

    y = NULL,

    title = "(a) Thermal microhabitat availability"
  ) +

    common_theme
}


############################################################################################
############################################################################################
#
# PANEL B
#
# VERTICAL DISTRIBUTION OF THERMAL NICHES
#
############################################################################################
############################################################################################

create_vertical_plot <- function(selected_period) {


  plot_dat <- vertical_summary %>%

    filter(
      period == selected_period
    )


  ggplot(

    plot_dat,

    aes(
      y = thermal_niche,
      colour = treatment
    )

  ) +

    # ---- P10-P90 vertical range ----
  #
  # Horizontal interval = central 80% vertical range

  geom_errorbarh(

    aes(
      xmin = p10_height,
      xmax = p90_height
    ),

    position = position_dodge(
      width = dodge_width
    ),

    height = 0,

    linewidth = 2,

    alpha = 0.30
  ) +

    # ---- Median height ----

  geom_point(

    aes(
      x = median_height
    ),

    position = position_dodge(
      width = dodge_width
    ),

    shape = 16,

    size = 5
  ) +

    # ---- Seasons ----

  facet_grid(
    . ~ seasonal_period
  ) +

    # ---- Treatment colours ----

  scale_colour_manual(
    values = treatment_colours,
    name = NULL
  ) +

    # ---- Niche labels ----

  scale_y_discrete(
    labels = niche_labels_coloured,
    drop = FALSE
  ) +

    # ---- Height axis ----

  scale_x_continuous(
    limits = c(-2, 42),
    breaks = seq(0, 40, by = 10),
    expand = c(0, 0)
  ) +

    # ---- Labels ----

  labs(

    x = "Height (m)",

    y = NULL,

    title = "(b) Vertical distribution of thermal microhabitats"
  ) +

    common_theme
}


############################################################################################
############################################################################################
#
# PANEL C
#
# SPATIAL CONNECTIVITY OF THERMAL NICHES
#
############################################################################################
############################################################################################

create_connectivity_plot <- function(selected_period) {


  plot_individual <- connectivity_plot %>%

    filter(
      period == selected_period
    )


  plot_median <- connectivity_summary %>%

    filter(
      period == selected_period
    )


  ggplot(

    plot_individual,

    aes(
      x = connectivity,
      y = thermal_niche,
      colour = treatment
    )

  ) +

    # ---- Individual plots ----

  geom_point(

    position = position_dodge(
      width = dodge_width
    ),

    size = 3,

    alpha = 0.25
  ) +

    # ---- Treatment median ----

  geom_point(

    data = plot_median,

    aes(
      x = median_connectivity,
      y = thermal_niche,
      colour = treatment
    ),

    position = position_dodge(
      width = dodge_width
    ),

    size = 5,

    inherit.aes = FALSE
  ) +

    # ---- Seasons ----

  facet_grid(
    . ~ seasonal_period
  ) +

    # ---- Treatment colours ----

  scale_colour_manual(
    values = treatment_colours,
    name = NULL
  ) +

    # ---- Niche labels ----

  scale_y_discrete(
    labels = niche_labels_coloured,
    drop = FALSE
  ) +

    # ---- Connectivity axis ----

  scale_x_continuous(
    limits = c(-0.1, 1.1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1),
    expand = c(0, 0)
  ) +

    # ---- Labels ----

  labs(

    x = "Local spatial microhabitat continuity",

    y = NULL,

    title = "(c) Local spatial continuity of thermal microhabitats"
  ) +

    common_theme
}


############################################################################################
############################################################################################
#
# MANUAL TREATMENT KEY
#
# This is deliberately drawn as a separate plot rather than using a ggplot
# legend. It therefore always appears ABOVE title (a).
#
############################################################################################
############################################################################################

treatment_legend_plot <- ggplot() +

  # ---- Control point ----

annotate(

  "point",

  x = 0.34,
  y = 0.5,

  colour = treatment_colours["Control"],

  size = 10
) +

  # ---- Control text ----

annotate(

  "text",

  x = 0.37,
  y = 0.5,

  label = "Control",

  hjust = 0,

  size = 10
) +

  # ---- Removal point ----

annotate(

  "point",

  x = 0.55,
  y = 0.5,

  colour = treatment_colours["Liana removal"],

  size = 10
) +

  # ---- Removal text ----

annotate(

  "text",

  x = 0.58,
  y = 0.5,

  label = "Liana removal",

  hjust = 0,

  size = 10
) +

  xlim(
    0,
    1
  ) +

  ylim(
    0,
    1
  ) +

  theme_void() +

  theme(

    plot.margin = margin(
      0,
      0,
      0,
      0
    )
  )


############################################################################################
############################################################################################
#
# DAYTIME MAIN FIGURE
#
############################################################################################
############################################################################################


# ==========================================================================================
# CREATE DAYTIME PANELS
# ==========================================================================================

plot_A_day <- create_availability_plot(
  "Day"
)

plot_B_day <- create_vertical_plot(
  "Day"
)

plot_C_day <- create_connectivity_plot(
  "Day"
)


# ==========================================================================================
# COMBINE DAYTIME PANELS
# ==========================================================================================

day_panels <- cowplot::plot_grid(

  plot_A_day,

  plot_B_day,

  plot_C_day,

  ncol = 1,

  rel_heights = c(
    1,
    1,
    1
  ),

  align = "v",

  axis = "lr"
)


# ==========================================================================================
# ADD TREATMENT KEY ABOVE TITLE (a)
# ==========================================================================================

main_niche_plot_day <- cowplot::plot_grid(

  treatment_legend_plot,

  day_panels,

  ncol = 1,

  rel_heights = c(
    0.055,
    1
  )
)


# ==========================================================================================
# SHOW
# ==========================================================================================

print(
  main_niche_plot_day
)


# ==========================================================================================
# SAVE
# ==========================================================================================

ggsave(

  filename = file.path(
    output_path,
    "thermal_niche_characteristics_day.png"
  ),

  plot = main_niche_plot_day,

  width = 14,

  height = 18,

  dpi = 300
)


############################################################################################
############################################################################################
#
# NIGHTTIME APPENDIX FIGURE
#
############################################################################################
############################################################################################


# ==========================================================================================
# CREATE NIGHTTIME PANELS
# ==========================================================================================

plot_A_night <- create_availability_plot(
  "Night"
)

plot_B_night <- create_vertical_plot(
  "Night"
)

plot_C_night <- create_connectivity_plot(
  "Night"
)


# ==========================================================================================
# COMBINE NIGHTTIME PANELS
# ==========================================================================================

night_panels <- cowplot::plot_grid(

  plot_A_night,

  plot_B_night,

  plot_C_night,

  ncol = 1,

  rel_heights = c(
    1,
    1,
    1
  ),

  align = "v",

  axis = "lr"
)


# ==========================================================================================
# ADD TREATMENT KEY ABOVE TITLE (a)
# ==========================================================================================

appendix_niche_plot_night <- cowplot::plot_grid(

  treatment_legend_plot,

  night_panels,

  ncol = 1,

  rel_heights = c(
    0.055,
    1
  )
)


# ==========================================================================================
# SHOW
# ==========================================================================================

print(
  appendix_niche_plot_night
)


# ==========================================================================================
# SAVE
# ==========================================================================================

ggsave(

  filename = file.path(
    output_path,
    "thermal_niche_characteristics_night.png"
  ),

  plot = appendix_niche_plot_night,

  width = 14,

  height = 18,

  dpi = 300
)


############################################################################################
############################################################################################
#
# NUMERICAL OUTPUT
#
############################################################################################
############################################################################################


# ==========================================================================================
# NICHE AVAILABILITY DIFFERENCE
#
# Removal - Control
#
# Units = percentage points of complete forest volume
# ==========================================================================================

availability_difference <- availability_summary %>%

  dplyr::select(
    treatment,
    seasonal_period,
    period,
    thermal_niche,
    median_availability
  ) %>%

  pivot_wider(

    names_from = treatment,

    values_from = median_availability
  ) %>%

  mutate(

    difference =
      `Liana removal` -
      Control
  )


cat(
  "\n\nNICHE AVAILABILITY\n"
)

cat(
  "Treatment difference = removal - control, percentage points\n\n"
)

print(
  availability_difference
)


# ==========================================================================================
# VERTICAL NICHE POSITION
# ==========================================================================================

cat(
  "\n\nVERTICAL NICHE DISTRIBUTION\n"
)

cat(
  "Median height and P10-P90 vertical range (m)\n\n"
)

print(
  vertical_summary
)


# ==========================================================================================
# CONNECTIVITY DIFFERENCE
#
# Removal - Control
# ==========================================================================================

connectivity_difference <- connectivity_summary %>%

  dplyr::select(
    treatment,
    seasonal_period,
    period,
    thermal_niche,
    median_connectivity
  ) %>%

  pivot_wider(

    names_from = treatment,

    values_from = median_connectivity
  ) %>%

  mutate(

    difference =
      `Liana removal` -
      Control
  )


cat(
  "\n\nSPATIAL NICHE CONNECTIVITY\n"
)

cat(
  "Treatment difference = removal - control\n\n"
)

print(
  connectivity_difference
)




# ==========================================================================================
# FINISHED
# ==========================================================================================

cat(
  "\nFigures saved to:\n",
  output_path,
  "\n"
)
