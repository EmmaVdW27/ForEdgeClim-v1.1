############################################################################################
# ForEdgeClim shortwave-radiation validation visualisation
#
# Compares shortwave-radiation validation performance between calibrated and
# uncalibrated ForEdgeClim simulations for one forest plot.
#
# The script visualises temporal Pearson and Spearman correlations across PAR
# sensors and spatial Pearson and Spearman correlations across days. Median
# correlations are calculated and reported for each model parameterisation.
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################

# ---- Packages ----
suppressPackageStartupMessages({
  library(ForEdgeClimOnPanama)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(patchwork)
  library(purrr)
  library(stringr)
})

# -----------------
# INPUT
# -----------------

# Set the base plot. The script will also look for "<base>_uncalibrated".
base_plot <- "R7"   # C8 or R7

# Build file map: names = model label, values = plot-string used in filenames
plot_map <- c(
  uncalibrated = paste0(base_plot, "_uncalibrated"),
  calibrated   = paste0(base_plot, "_calibrated")

)

# -----------------
# OUTPUT
# -----------------

output_path <- "Output/validation/plots/"

# -----------------
# LOAD DATA (both calibrated & uncalibrated)
# -----------------
temporal_uncal <- readRDS(
  paste0(
    "Output/validation/data/temporal_metrics_",
    plot_map[1],
    ".rds"
  )
)

temporal_cal <- readRDS(
  paste0(
    "Output/validation/data/temporal_metrics_",
    plot_map[2],
    ".rds"
  )
)

spatial_uncal <- readRDS(
  paste0(
    "Output/validation/data/spatial_metrics_",
    plot_map[1],
    ".rds"
  )
)

spatial_cal <- readRDS(
  paste0(
    "Output/validation/data/spatial_metrics_",
    plot_map[2],
    ".rds"
  )
)


# ==============================================================================
# SPATIAL CORRELATIONS — DRY AND WET SEASONS
# ==============================================================================

spatial_seasonal <- bind_rows(
  spatial_uncal %>%
    mutate(model = "uncalibrated"),

  spatial_cal %>%
    mutate(model = "calibrated")
) %>%
  mutate(
    season = if_else(
      lubridate::month(date) %in% 1:4,
      "Dry",
      "Wet"
    )
  )

spatial_seasonal_stats <- spatial_seasonal %>%
  group_by(model, season) %>%
  summarise(
    n_days = n_distinct(date),
    spatial_pearson = median(
      spatial_pearson,
      na.rm = TRUE
    ),
    spatial_spearman = median(
      spatial_spearman,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(model, season)

cat("\nSpatial correlations by season:\n")
print(spatial_seasonal_stats)

write.csv(
  spatial_seasonal_stats,
  file.path(
    output_path,
    paste0("spatial_seasonal_", base_plot, ".csv")
  ),
  row.names = FALSE
)


temporal_all <- bind_rows(

  temporal_uncal %>%
    mutate(model = "uncalibrated"),

  temporal_cal %>%
    mutate(model = "calibrated")

) %>%
  pivot_longer(
    cols = c(
      temporal_pearson,
      temporal_spearman
    ),
    names_to = "metric",
    values_to = "value"
  )

spatial_all <- bind_rows(

  spatial_uncal %>%
    mutate(model = "uncalibrated"),

  spatial_cal %>%
    mutate(model = "calibrated")

) %>%
  pivot_longer(
    cols = c(
      spatial_pearson,
      spatial_spearman
    ),
    names_to = "metric",
    values_to = "value"
  )

temporal_all$model <- factor(
  temporal_all$model,
  levels = c("uncalibrated", "calibrated")
)

spatial_all$model <- factor(
  spatial_all$model,
  levels = c("uncalibrated", "calibrated")
)

# -----------------
# STATISTICS (per model)
# -----------------

stats_df <- bind_rows(

  temporal_uncal %>%
    summarise(

      temporal_pearson =
        median(
          temporal_pearson,
          na.rm = TRUE
        ),

      temporal_spearman =
        median(
          temporal_spearman,
          na.rm = TRUE
        )

    ) %>%
    mutate(model = "uncalibrated"),

  temporal_cal %>%
    summarise(

      temporal_pearson =
        median(
          temporal_pearson,
          na.rm = TRUE
        ),

      temporal_spearman =
        median(
          temporal_spearman,
          na.rm = TRUE
        )

    ) %>%
    mutate(model = "calibrated")

) %>%
  left_join(

    bind_rows(

      spatial_uncal %>%
        summarise(

          spatial_pearson =
            median(
              spatial_pearson,
              na.rm = TRUE
            ),

          spatial_spearman =
            median(
              spatial_spearman,
              na.rm = TRUE
            )

        ) %>%
        mutate(model = "uncalibrated"),

      spatial_cal %>%
        summarise(

          spatial_pearson =
            median(
              spatial_pearson,
              na.rm = TRUE
            ),

          spatial_spearman =
            median(
              spatial_spearman,
              na.rm = TRUE
            )

        ) %>%
        mutate(model = "calibrated")

    ),

    by = "model"

  )

# plot labels:
temporal_all <- temporal_all %>%
  mutate(
    metric = recode(
      metric,
      temporal_pearson  = "Pearson",
      temporal_spearman = "Spearman"
    )
  )

spatial_all <- spatial_all %>%
  mutate(
    metric = recode(
      metric,
      spatial_pearson  = "Pearson",
      spatial_spearman = "Spearman"
    )
  )

temporal_medians <- temporal_all %>%
  group_by(model, metric) %>%
  summarise(
    median_value = median(value, na.rm = TRUE),
    .groups = "drop"
  )

spatial_medians <- spatial_all %>%
  group_by(model, metric) %>%
  summarise(
    median_value = median(value, na.rm = TRUE),
    .groups = "drop"
  )


# -----------------
# TEMPORAL PLOTS
# -----------------

temporal_plot <- ggplot(
  temporal_all,
  aes(
    x = metric,
    y = value
  )
) +
  geom_boxplot(
    fill = "white",
    colour = "blue",
    linewidth = 1
  ) +
  facet_wrap(
    ~ model,
    ncol = 2
  ) +
  theme_bw() +
  coord_cartesian(
    ylim = c(-1, 1)
  ) +
  labs(
    x = "",
    y = "Correlation",
    title = "(c) Temporal model-observation light agreement"
  ) +
  theme(
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 35),
    strip.text = element_text(size = 40)
  ) +
  geom_label(
    data = temporal_medians,
    aes(
      x = metric,
      y = median_value,
      label = sprintf("%.2f", median_value)
    ),
    inherit.aes = FALSE,
    fill = scales::alpha("white", 0.8),
    colour = "blue",
    label.size = 0,
    size = 13,
    fontface = "bold"
  )



# -----------------
# SPATIAL PLOTS
# -----------------
spatial_plot <- ggplot(
  spatial_all,
  aes(
    x = metric,
    y = value
  )
) +
  geom_boxplot(
    fill = "white",
    colour = "blue",
    linewidth = 1
  ) +
  facet_wrap(
    ~ model,
    ncol = 2
  ) +
  theme_bw() +
  coord_cartesian(
    ylim = c(-1, 1)
  ) +
  labs(
    x = "",
    y = "Correlation",
    title = "(d) Spatial model-observation light agreement"
  ) +
  theme(
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 35),
    strip.text = element_text(size = 40)
  ) +
  geom_label(
    data = spatial_medians,
    aes(
      x = metric,
      y = median_value,
      label = sprintf("%.2f", median_value)
    ),
    inherit.aes = FALSE,
    fill = scales::alpha("white", 0.8),
    colour = "blue",
    label.size = 0,
    size = 13,
    fontface = "bold"
  )

ggsave(
  file.path(
    output_path,
    paste0(
      "temporal_validation_",
      base_plot,
      "_both.png"
    )
  ),
  temporal_plot,
  width = 16,
  height = 10,
  dpi = 300
)

ggsave(
  file.path(
    output_path,
    paste0(
      "spatial_validation_",
      base_plot,
      "_both.png"
    )
  ),
  spatial_plot,
  width = 16,
  height = 10,
  dpi = 300
)


# -----------------
# CONSOLE SUMMARIES (per model)
# -----------------

stats_df %>%
  rename(
    Model = model,
    `Temporal Pearson` = temporal_pearson,
    `Temporal Spearman` = temporal_spearman,
    `Spatial Pearson` = spatial_pearson,
    `Spatial Spearman` = spatial_spearman
  ) %>%
  print(row.names = FALSE)

