################################################################################
# ForEdgeClim - Vertical absolute air-temperature profiles and spatial P05-P95
# Day and night shown separately; dry and wet seasons in columns.
#
# Author: Emma Van de Walle – Q-ForestLab
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(ggh4x)
  library(patchwork)
  library(lmerTest)
  library(emmeans)
  library(purrr)
  library(tidyr)
})

control <- readRDS("Output/heterogeneity/heterogeneity_control_calibrated.rds") %>%
  mutate(treatment = "Control")
removal <- readRDS("Output/heterogeneity/heterogeneity_removal_calibrated.rds") %>%
  mutate(treatment = "Liana removal")

output_path <- "Output/heterogeneity/figures"
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Keep original seasonal and daytime definitions from the supplied script.
dat <- bind_rows(control, removal) %>%
  mutate(
    treatment = factor(treatment, levels = c("Control", "Liana removal")),
    season = if_else(doy >= 350 | doy <= 120, "Dry season", "Wet season"),
    period = if_else(hour >= 11 & hour <= 22, "Day", "Night"),
    season = factor(season, levels = c("Dry season", "Wet season")),
    period = factor(period, levels = c("Day", "Night"))
  ) %>%
  filter(is.finite(Tair_p05), is.finite(Tair_p95), is.finite(Tair_mean))

# Step 1: summarize hourly spatial percentiles through time within each plot.
# These are medians of hourly spatial percentiles, NOT percentiles of all
# voxel temperatures pooled across the entire season.
plot_summary <- dat %>%
  group_by(treatment, plot, season, period, z) %>%
  summarise(
    p05 = median(Tair_p05, na.rm = TRUE),
    p95 = median(Tair_p95, na.rm = TRUE),
    mean_T = median(Tair_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(p05), is.finite(p95), is.finite(mean_T))

# Step 2: summarize plot-specific profiles equally across plots.
vertical_summary <- plot_summary %>%
  group_by(treatment, season, period, z) %>%
  summarise(
    p05 = median(p05, na.rm = TRUE),
    p95 = median(p95, na.rm = TRUE),
    mean_T = median(mean_T, na.rm = TRUE),
    n_plots = n_distinct(plot),
    .groups = "drop"
  ) %>%
  filter(is.finite(p05), is.finite(p95), is.finite(mean_T))

# Set to NULL for automatically selected shared x-axis limits.
# Alternatively specify e.g. c(20, 34) for both day and night.
manual_x_limits <- NULL

treatment_colours <- c("Control" = "blue", "Liana removal" = "red")

make_profile <- function(selected_period) {
  df <- vertical_summary %>% filter(period == selected_period)

  p <- ggplot(df, aes(y = z, colour = treatment, fill = treatment)) +
    geom_ribbon(
      aes(xmin = p05, xmax = p95, group = treatment),
      orientation = "y", alpha = 0.18, colour = NA
    ) +
    geom_path(aes(x = mean_T, group = treatment), linewidth = 1.2) +
    facet_grid(. ~ season) +
    scale_x_continuous(
      name = expression("Air temperature (" * degree * "C)"),
      breaks = function(x) seq(
        from = ceiling(x[1]),
        to = floor(x[2]),
        by = 1
      ),
      labels = scales::label_number(accuracy = 1),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    scale_y_continuous(
      name = "Height (m)", limits = c(0, 40),
      breaks = seq(0, 40, 10), expand = expansion(mult = c(0, 0))
    ) +
    scale_colour_manual(values = treatment_colours, name = NULL) +
    scale_fill_manual(values = treatment_colours, name = NULL) +
    guides(fill = "none", colour = guide_legend(override.aes = list(linewidth = 1.5))) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 21),
      axis.text = element_text(size = 17),
      strip.text = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18),
      legend.position = "top",
      legend.justification = "left",
      legend.key.width = unit(1.5, "cm"),
      plot.margin = margin(10, 15, 10, 10)
    )

  if (!is.null(manual_x_limits)) {
    p <- p + coord_cartesian(xlim = manual_x_limits)
  }
  p
}

plot_day <- make_profile("Day")
plot_night <- make_profile("Night")

print(plot_day)
print(plot_night)

ggsave(file.path(output_path, "vertical_absolute_temperature_profiles_day.png"),
       plot_day, width = 9, height = 5, dpi = 300)
ggsave(file.path(output_path, "vertical_absolute_temperature_profiles_night.png"),
       plot_night, width = 9, height = 5, dpi = 300)




############################################################################################
#
# APPENDIX FIGURE
#
# DISTRIBUTION OF SPATIAL THERMAL HETEROGENEITY
# ACROSS VERTICAL FOREST ZONES
#
############################################################################################


# ==========================================================================================
# DEFINE VERTICAL FOREST ZONES
#
# 0-2 m
# 2-10 m
# 10-20 m
# 20-30 m
# >30 m
# ==========================================================================================

vertical_zone_data <- dat %>%

  mutate(

    height_class = cut(

      z,

      breaks = c(
        -Inf,
        2,
        10,
        20,
        30,
        Inf
      ),

      labels = c(
        "0-2 m",
        "2-10 m",
        "10-20 m",
        "20-30 m",
        ">30 m"
      ),

      right = FALSE
    )
  )

# Height classes ordered from forest canopy to forest floor
height_levels <- c(
  ">30 m",
  "20-30 m",
  "10-20 m",
  "2-10 m",
  "0-2 m"
)

vertical_zone_data <- vertical_zone_data %>%
  mutate(
    height_class = factor(
      as.character(height_class),
      levels = height_levels
    )
  )

# ---- Remove missing height classes ----

vertical_zone_data <- vertical_zone_data[
  !is.na(vertical_zone_data$height_class),
]


# ==========================================================================================
# STEP 1
#
# Average hourly heterogeneity within each:
#
# plot x day x period x 1-m layer
# ==========================================================================================

vertical_zone_data <- vertical_zone_data %>%

  group_by(
    treatment,
    plot,
    season,
    period,
    date,
    height_class,
    z
  ) %>%

  summarise(

    Tair_range90 = mean(
      Tair_range90,
      na.rm = TRUE
    ),

    .groups = "drop"
  )


# ==========================================================================================
# STEP 2
#
# Average 1-m layers within each vertical forest zone.
#
# This prevents wider zones from receiving more weight simply because
# they contain more individual 1-m layers.
# ==========================================================================================

vertical_zone_data <- vertical_zone_data %>%

  group_by(
    treatment,
    plot,
    season,
    period,
    date,
    height_class
  ) %>%

  summarise(

    Tair_range90 = mean(
      Tair_range90,
      na.rm = TRUE
    ),

    .groups = "drop"
  )


# ---- Remove non-finite observations ----

vertical_zone_data <- vertical_zone_data[
  is.finite(vertical_zone_data$Tair_range90),
]


# ==========================================================================================
# CREATE TREATMENT POSITIONS
#
# Control       = position 1
# Liana removal = position 2
#
# Treatment names are deliberately not shown on the x-axis.
# The colours and legend identify the treatments.
# ==========================================================================================

vertical_zone_data <- vertical_zone_data %>%

  mutate(

    treatment_position = case_when(

      treatment == "Control" ~ 1,

      treatment == "Liana removal" ~ 2
    )
  )



# ==========================================================================================
# APPENDIX VIOLIN PLOT
# Day:   Dry | Wet  (shared y-axis within each height class)
# Night: Dry | Wet  (shared y-axis within each height class)
# ==========================================================================================

make_zone_plot <- function(selected_period) {

  df <- vertical_zone_data %>%
    filter(period == selected_period) %>%
    mutate(
      season = factor(
        season,
        levels = c("Dry season", "Wet season")
      )
    )

  ggplot(
    df,
    aes(
      x = treatment_position,
      y = Tair_range90,
      fill = treatment,
      colour = treatment,
      group = treatment
    )
  ) +

    geom_violin(
      trim = TRUE,
      alpha = 0.25,
      scale = "width",
      linewidth = 1.3,
      width = 0.8
    ) +

    geom_boxplot(
      width = 0.16,
      outlier.shape = NA,
      alpha = 0.65,
      linewidth = 1
    ) +

    facet_grid2(
      height_class ~ season,
      scales = "free_y",
      independent = "none",
      axes = "margins",
      switch = NULL,
      labeller = labeller(
        season = function(x) {
          paste(selected_period, "-", x)
        }
      )
    ) +

    scale_fill_manual(
      values = treatment_colours,
      name = NULL
    ) +

    scale_colour_manual(
      values = treatment_colours,
      name = NULL
    ) +

    scale_x_continuous(
      limits = c(0.4, 2.6),
      breaks = NULL
    ) +

    labs(
      x = NULL,
      y = if (selected_period == "Day") {
        "Spatial thermal range (P95 - P05, °C)"
      } else {
        NULL
      }
    ) +

    theme_bw() +

    theme(
      axis.title = element_text(size = 32),
      axis.text.y = element_text(size = 28),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),

      strip.text = element_text(
        size = 30,
        face = "bold"
      ),

      strip.text.y.right = if (selected_period == "Night") {
        element_text(size = 26, face = "bold", angle = 270)
      } else {
        element_blank()
      },

      legend.text = element_text(size = 28),
      legend.position = "none",

      panel.spacing.y = unit(1.2, "cm"),
      panel.spacing.x = unit(0.3, "cm"),

      plot.margin = margin(10, 10, 10, 10)
    )
}

# Separate day and night panels
plot_zones_day <- make_zone_plot("Day")
plot_zones_night <- make_zone_plot("Night")

# Combine into four columns:
# Day Dry | Day Wet | Night Dry | Night Wet
vertical_zone_plot <-
  (plot_zones_day | plot_zones_night) +
  plot_layout(
    widths = c(1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 28),
    legend.key.width = unit(2.5, "cm"),
    legend.key.spacing.x = unit(1.5, "cm")
  )

# Save
ggsave(
  filename = file.path(
    output_path,
    "vertical_zone_thermal_heterogeneity_distributions.png"
  ),
  plot = vertical_zone_plot,
  width = 20,
  height = 19,
  dpi = 300
)



################################################################################
# ADDITIONAL FIGURE: ABSOLUTE TEMPERATURE PROFILES AT RELATIVE CANOPY HEIGHT
# Same formatting as vertical_absolute_temperature_profiles_day.png
################################################################################

canopy_heights <- tibble(
  plot = c("C8", "C13", "C12", "C16",
           "R7", "R11", "R5", "R14"),
  H_max = c(38, 33, 29, 30,
            37, 35, 36, 42)
)

# Retain existing canopy_heights table with columns plot and H_max.

# Step 1: Normalise height for each plot.
relative_dat <- dat %>%
  left_join(canopy_heights, by = "plot") %>%
  mutate(z_rel = z / H_max) %>%
  filter(
    is.finite(z_rel),
    z_rel >= 0,
    z_rel <= 1
  )

# Step 2: Summarise hourly temperature statistics within each plot.
relative_plot_summary <- relative_dat %>%
  group_by(treatment, plot, season, period, z_rel) %>%
  summarise(
    p05 = median(Tair_p05, na.rm = TRUE),
    p95 = median(Tair_p95, na.rm = TRUE),
    mean_T = median(Tair_mean, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Interpolate onto a common relative-height grid.
common_heights <- seq(0.05, 0.95, by = 0.025)

relative_interpolated <- relative_plot_summary %>%
  group_by(treatment, plot, season, period) %>%
  group_modify(~ {

    d <- arrange(.x, z_rel)

    tibble(
      z_rel = common_heights,
      p05 = approx(
        d$z_rel, d$p05,
        xout = common_heights,
        rule = 1
      )$y,
      p95 = approx(
        d$z_rel, d$p95,
        xout = common_heights,
        rule = 1
      )$y,
      mean_T = approx(
        d$z_rel, d$mean_T,
        xout = common_heights,
        rule = 1
      )$y
    )
  }) %>%
  ungroup() %>%
  filter(
    is.finite(p05),
    is.finite(p95),
    is.finite(mean_T)
  )

# Step 4: Summarise equally across plots, as in the original script.
relative_summary <- relative_interpolated %>%
  group_by(treatment, season, period, z_rel) %>%
  summarise(
    p05 = median(p05, na.rm = TRUE),
    p95 = median(p95, na.rm = TRUE),
    mean_T = median(mean_T, na.rm = TRUE),
    n_plots = n_distinct(plot),
    .groups = "drop"
  ) %>%
  filter(
    n_plots == 4,
    is.finite(p05),
    is.finite(p95),
    is.finite(mean_T)
  )

# Step 5: Create the figure using the exact original ggplot formatting.
plot_relative_day <- ggplot(
  relative_summary %>% filter(period == "Day"),
  aes(
    y = z_rel,
    colour = treatment,
    fill = treatment
  )
) +
  geom_ribbon(
    aes(
      xmin = p05,
      xmax = p95,
      group = treatment
    ),
    orientation = "y",
    alpha = 0.18,
    colour = NA
  ) +
  geom_path(
    aes(
      x = mean_T,
      group = treatment
    ),
    linewidth = 1.2
  ) +
  facet_grid(. ~ season) +
  scale_x_continuous(
    name = expression("Air temperature (" * degree * "C)"),
    breaks = function(x) seq(
      from = ceiling(x[1]),
      to = floor(x[2]),
      by = 1
    ),
    labels = scales::label_number(accuracy = 1),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_y_continuous(
    name = "Relative canopy height",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_colour_manual(
    values = treatment_colours,
    name = NULL
  ) +
  scale_fill_manual(
    values = treatment_colours,
    name = NULL
  ) +
  guides(
    fill = "none",
    colour = guide_legend(
      override.aes = list(linewidth = 1.5)
    )
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 21),
    axis.text = element_text(size = 17),
    strip.text = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "top",
    legend.justification = "left",
    legend.key.width = unit(1.5, "cm"),
    plot.margin = margin(10, 15, 10, 10)
  )

print(plot_relative_day)

ggsave(
  file.path(
    output_path,
    "vertical_absolute_temperature_profiles_relative_height_day.png"
  ),
  plot_relative_day,
  width = 9,
  height = 5,
  dpi = 300
)


################################################################################
# QUANTITATIVE COMPARISON OF VERTICAL THERMAL HETEROGENEITY
# Absolute height versus relative canopy height
################################################################################


# Common absolute canopy-height range across all plots.
H_common <- min(canopy_heights$H_max)

# ---- 1. Absolute-height treatment differences ----

absolute_ranges <- plot_summary %>%
  mutate(
    thermal_range = p95 - p05
  ) %>%
  filter(z >= 0, z <= H_common) %>%
  group_by(season, period, treatment, z) %>%
  summarise(
    thermal_range = median(thermal_range, na.rm = TRUE),
    n_plots = n_distinct(plot),
    .groups = "drop"
  ) %>%
  filter(n_plots == 4) %>%
  select(season, period, treatment, z, thermal_range) %>%
  pivot_wider(
    names_from = treatment,
    values_from = thermal_range
  ) %>%
  mutate(
    delta = `Liana removal` - Control,
    abs_delta = abs(delta),
    height_system = "Absolute height",
    height_scaled = z / H_common
  )

# ---- 2. Relative-height treatment differences ----

relative_ranges <- relative_interpolated %>%
  mutate(
    thermal_range = p95 - p05
  ) %>%
  group_by(season, period, treatment, z_rel) %>%
  summarise(
    thermal_range = median(thermal_range, na.rm = TRUE),
    n_plots = n_distinct(plot),
    .groups = "drop"
  ) %>%
  filter(n_plots == 4) %>%
  select(season, period, treatment, z_rel, thermal_range) %>%
  pivot_wider(
    names_from = treatment,
    values_from = thermal_range
  ) %>%
  mutate(
    delta = `Liana removal` - Control,
    abs_delta = abs(delta),
    height_system = "Relative canopy height",
    height_scaled = z_rel
  )

# ---- 3. Summarise treatment contrasts ----

summarise_contrasts <- function(df) {

  df %>%
    group_by(season, period) %>%
    summarise(
      mean_abs_delta = mean(abs_delta, na.rm = TRUE),
      median_abs_delta = median(abs_delta, na.rm = TRUE),
      mean_signed_delta = mean(delta, na.rm = TRUE),
      max_abs_delta = max(abs_delta, na.rm = TRUE),
      n_heights = n(),
      .groups = "drop"
    )
}

summary_absolute <- summarise_contrasts(absolute_ranges) %>%
  rename_with(
    ~ paste0(.x, "_absolute"),
    -c(season, period)
  )

summary_relative <- summarise_contrasts(relative_ranges) %>%
  rename_with(
    ~ paste0(.x, "_relative"),
    -c(season, period)
  )

comparison <- summary_absolute %>%
  left_join(
    summary_relative,
    by = c("season", "period")
  ) %>%
  mutate(
    change_abs_delta_pct =
      100 * (
        mean_abs_delta_relative /
          mean_abs_delta_absolute - 1
      )
  )

print(comparison, width = Inf)

write.csv(
  comparison,
  file.path(
    output_path,
    "vertical_heterogeneity_height_normalisation_metrics.csv"
  ),
  row.names = FALSE
)

# ---- 4. Plot treatment contrasts ----

contrast_profiles <- bind_rows(
  absolute_ranges %>%
    select(
      season, period, height_system,
      height_scaled, delta
    ),
  relative_ranges %>%
    select(
      season, period, height_system,
      height_scaled, delta
    )
) %>%
  mutate(
    height_system = factor(
      height_system,
      levels = c(
        "Absolute height",
        "Relative canopy height"
      )
    )
  )

plot_contrasts <- ggplot(
  contrast_profiles %>% filter(period == "Day"),
  aes(
    x = delta,
    y = height_scaled,
    colour = height_system
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_path(linewidth = 1.2) +
  facet_grid(. ~ season) +
  scale_y_continuous(
    name = "Normalised vertical coordinate",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = expression(
      "Treatment difference in spatial thermal range (" *
        degree * "C)"
    ),
    colour = NULL
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 21),
    axis.text = element_text(size = 17),
    strip.text = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "top",
    legend.justification = "left",
    legend.key.width = unit(1.5, "cm"),
    plot.margin = margin(10, 15, 10, 10)
  )

print(plot_contrasts)

ggsave(
  file.path(
    output_path,
    "vertical_heterogeneity_height_normalisation_comparison.png"
  ),
  plot_contrasts,
  width = 9,
  height = 5,
  dpi = 300
)

