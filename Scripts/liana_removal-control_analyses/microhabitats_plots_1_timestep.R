############################################################################################
# ForEdgeClim thermal microhabitat (= niche) visualisation - C8 and R7
#
# Date/time: 26 April 2023, 19:00 UTC
#
# Creates for each plot:
#   1. North-south averaged vertical thermal niche profile (x-z)
#   2. Ground-floor thermal niche map (x-y)
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################

# ---- Packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

# -----------------
# INPUT
# -----------------

average_over_slices = "no" # yes or no (no = single slice)

file_C8 <- "Output/2023-04-26_19/C8/model_results.rds"
file_R7 <- "Output/2023-04-26_19/R7/model_results.rds"

output_path <- "Output/thermal_niches/plots/2023-04-26_19"

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# -----------------
# THERMAL NICHE DEFINITIONS
# -----------------

# T_anomaly = voxel Tair - mean Tair of complete 3D forest volume
#
# Very cool       < -1.0 °C
# Cool            -1.0 to -0.5 °C
# Near-average    -0.5 to +0.5 °C
# Warm            +0.5 to +1.0 °C
# Very warm       > +1.0 °C

niche_levels <- c("Very cool", "Cool", "Near-average", "Warm", "Very warm")

niche_colours <- c(
  "Very cool"    = "blue",
  "Cool"         = "lightblue",
  "Near-average" = "bisque2",
  "Warm"         = "orange",
  "Very warm"    = "red"
)


classify_thermal_niche <- function(T_anomaly) {
  cut(
    T_anomaly,
    breaks = c(-Inf, -1, -0.5, 0.5, 1, Inf),
    labels = niche_levels,
    right = FALSE
  )
}

# -----------------
# LOAD MODEL RESULTS
# -----------------

model_C8 <- readRDS(file_C8)
model_R7 <- readRDS(file_R7)

# -----------------
# PREPARE 3D TEMPERATURE DATA
# -----------------

prepare_temperature_grid <- function(model_result, plot_name) {

  micro_grid <- model_result$micro_grid
  air_temp <- model_result$air_temperature

  if (nrow(micro_grid) != length(air_temp)) {
    stop(paste("micro_grid and air_temperature have different lengths for", plot_name))
  }

  temp_grid <- micro_grid %>%
    transmute(
      x = x,
      y = y,
      z = z,
      Tair = air_temp - 273.15
    )

  Tair_volume_mean <- mean(temp_grid$Tair, na.rm = TRUE)

  temp_grid <- temp_grid %>%
    mutate(
      T_anomaly = Tair - Tair_volume_mean,
      thermal_niche = classify_thermal_niche(T_anomaly)
    )

  cat(
    "\n", plot_name,
    "\nMean 3D air temperature = ", round(Tair_volume_mean, 2), " °C\n",
    sep = ""
  )

  return(temp_grid)
}

temp_C8 <- prepare_temperature_grid(model_C8, "C8")
temp_R7 <- prepare_temperature_grid(model_R7, "R7")

# -----------------
# COMMON AXIS LIMITS
# -----------------

# Same dimensions are used for C8 and R7 so that the panels can easily
# be aligned and transformed in Inkscape.

x_limits <- range(c(temp_C8$x, temp_R7$x), na.rm = TRUE)
y_limits <- range(c(temp_C8$y, temp_R7$y), na.rm = TRUE)
z_limits <- range(c(temp_C8$z, temp_R7$z), na.rm = TRUE)

cat(
  "\nAxis limits:\n",
  "x: ", paste(x_limits, collapse = " - "), "\n",
  "y: ", paste(y_limits, collapse = " - "), "\n",
  "z: ", paste(z_limits, collapse = " - "), "\n",
  sep = ""
)

# -----------------
# GROUND-FLOOR TEMPERATURE
# -----------------

prepare_ground_floor <- function(temp_grid) {

  ground_z <- min(temp_grid$z, na.rm = TRUE)

  temp_grid %>%
    filter(z == ground_z)
}

ground_C8 <- prepare_ground_floor(temp_C8)
ground_R7 <- prepare_ground_floor(temp_R7)


if(average_over_slices == "yes"){
  # -----------------
  # NORTH-SOUTH AVERAGED TEMPERATURE
  # -----------------

  # Air temperature is first averaged over the north-south (y) direction
  # for every x-z position.
  #
  # The resulting mean temperature is then compared with the mean temperature
  # of the complete original 3D forest volume.

  prepare_vertical_profile <- function(temp_grid) {

    Tair_volume_mean <- mean(temp_grid$Tair, na.rm = TRUE)

    temp_grid %>%
      group_by(x, z) %>%
      summarise(mean_temperature = mean(Tair, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        T_anomaly = mean_temperature - Tair_volume_mean,
        thermal_niche = classify_thermal_niche(T_anomaly)
      )
  }

  average_air_temperature_C8 <- prepare_vertical_profile(temp_C8)
  average_air_temperature_R7 <- prepare_vertical_profile(temp_R7)
} else{
  # -----------------
  # SINGLE VERTICAL SLICE AT Y = 30 m
  # -----------------

  # Select a single x-z slice at y = 30 m.
  # Thermal anomalies remain relative to the mean air temperature
  # of the complete original 3D forest volume.

  selected_y <- 30

  prepare_vertical_profile <- function(temp_grid, y_slice = selected_y) {

    Tair_volume_mean <- mean(temp_grid$Tair, na.rm = TRUE)

    # Select nearest available y-coordinate
    available_y <- sort(unique(temp_grid$y))
    actual_y <- available_y[which.min(abs(available_y - y_slice))]

    cat(
      "\nRequested y = ", y_slice,
      " m; selected y = ", actual_y,
      " m\n",
      sep = ""
    )

    temp_grid %>%
      filter(y == actual_y) %>%
      mutate(
        T_anomaly = Tair - Tair_volume_mean,
        thermal_niche = classify_thermal_niche(T_anomaly)
      )
  }

  slice_air_temperature_C8 <- prepare_vertical_profile(temp_C8)
  slice_air_temperature_R7 <- prepare_vertical_profile(temp_R7)
}

# -----------------
# VERTICAL NICHE PLOT
# -----------------

plot_vertical_niches <- function(dat) {

  ggplot(dat, aes(x = x, y = z, fill = thermal_niche)) +
    geom_tile() +
    scale_fill_manual(
      values = niche_colours,
      limits = niche_levels,
      drop = FALSE
    ) +
    labs(title = NULL, x = NULL, y = NULL, fill = NULL) +
    scale_y_continuous(position = "right") +
    coord_fixed(ratio = 1, clip = "off") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 20),
      legend.position = "none",
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA)
    )
}


# -----------------
# GROUND-FLOOR NICHE PLOT
# -----------------

plot_ground_niches <- function(dat) {

  ggplot(dat, aes(x = x, y = y, fill = thermal_niche)) +
    geom_tile() +
    scale_fill_manual(
      values = niche_colours,
      limits = niche_levels,
      drop = FALSE
    ) +
    labs(title = NULL, x = NULL, y = NULL, fill = NULL) +
    coord_fixed(ratio = 1, clip = "off") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.position = "none",
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA)
    )
}


# -----------------
# CREATE PLOTS
# -----------------
if(average_over_slices == "yes"){

  temp_air_niche_3D_C8 <- plot_vertical_niches(average_air_temperature_C8)
  temp_air_niche_3D_ground_C8 <- plot_ground_niches(ground_C8)

  temp_air_niche_3D_R7 <- plot_vertical_niches(average_air_temperature_R7)
  temp_air_niche_3D_ground_R7 <- plot_ground_niches(ground_R7)

} else{
  temp_air_niche_3D_C8 <- plot_vertical_niches(slice_air_temperature_C8)
  temp_air_niche_3D_ground_C8 <- plot_ground_niches(ground_C8)

  temp_air_niche_3D_R7 <- plot_vertical_niches(slice_air_temperature_R7)
  temp_air_niche_3D_ground_R7 <- plot_ground_niches(ground_R7)
}


# C8
ggsave(
  file.path(output_path, "3D_temp_air_niches_C8.pdf"),
  plot = temp_air_niche_3D_C8,
  width = 3, height = 2, dpi = 500
)

ggsave(
  file.path(output_path, "3D_temp_air_niches_ground_C8.pdf"),
  plot = temp_air_niche_3D_ground_C8,
  width = 3, height = 3, dpi = 500
)

# R7
ggsave(
  file.path(output_path, "3D_temp_air_niches_R7.pdf"),
  plot = temp_air_niche_3D_R7,
  width = 3, height = 2, dpi = 500
)

ggsave(
  file.path(output_path, "3D_temp_air_niches_ground_R7.pdf"),
  plot = temp_air_niche_3D_ground_R7,
  width = 3, height = 3, dpi = 500
)

# -----------------
# SEPARATE THERMAL NICHE LEGEND
# -----------------

legend_plot <- ggplot(
  data.frame(
    x = 1:5,
    y = 1,
    thermal_niche = factor(niche_levels, levels = niche_levels)
  ),
  aes(x = x, y = y, fill = thermal_niche)
) +
  geom_tile() +
  scale_fill_manual(
    values = niche_colours,
    limits = niche_levels,
    drop = FALSE,
    name = NULL
  ) +
  guides(
    fill = guide_legend(
      nrow = 5,
      byrow = TRUE,
      keyheight = grid::unit(0.7, "cm"),
      keywidth = grid::unit(0.7, "cm")
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 16),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )

# Extract only the legend
legend_only <- cowplot::get_legend(legend_plot)

# Save as separate PDF
ggsave(
  file.path(output_path, "thermal_niche_legend.pdf"),
  plot = legend_only,
  width = 3.2,
  height = 3,
  bg = "transparent"
)

cat("\nThermal niche plots saved to:\n", output_path, "\n")
