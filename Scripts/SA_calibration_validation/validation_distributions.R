################################################################################
# ForEdgeClim - calibrated spatial distribution validation using HOURLY TEMPERATURE and DAILY-MEAN RADIATION
# 2023; temperature C13/R11, radiation C8/R7; both 08:00-16:00 Panama LT
#
# Author: Emma Van de Walle - Q-ForestLab
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(patchwork)
})

temperature_plots <- c(Control = "C13", Removal = "R11")
radiation_plots <- c(Control = "C8", Removal = "R7")
temperature_dir <- "Output/validation/data"
radiation_dir <- "Output/timeseries/data"
output_dir <- "Output/validation/combined_spatial_distributions"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
local_hours <- 08:16  # inclusive; Panama is UTC-5 year-round
min_hours <- 6L
min_sensors <- 10L
min_mean_PAR <- 5
kde_adjust <- 2
colours <- c(Observed = "blue", Modelled = "red")

wasserstein_1d <- function(x, y) {
  x <- sort(x[is.finite(x)]); y <- sort(y[is.finite(y)])
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  z <- sort(unique(c(x, y)))
  if (length(z) < 2) return(0)
  mid <- (head(z, -1) + tail(z, -1)) / 2
  sum(abs(findInterval(mid, x) / length(x) -
            findInterval(mid, y) / length(y)) * diff(z))
}
q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))

# TEMPERATURE: retain hourly values per sensor, only 08:00--16:00 Panama LT;
# no daily aggregation
read_temperature <- function(treatment, plot_id) {
  path <- file.path(temperature_dir, paste0("val_data_", plot_id, "_calibrated.rds"))
  if (!file.exists(path)) stop("Missing file: ", path)
  dat <- readRDS(path)
  required <- c("datetime", "obs", "mod", "TMS_x", "TMS_y")
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  dat %>%
    transmute(treatment = treatment, plot = plot_id,
              datetime = as.POSIXct(datetime, tz = "UTC"),
              sensor_x = as.numeric(TMS_x), sensor_y = as.numeric(TMS_y),
              obs = as.numeric(obs), mod = as.numeric(mod)) %>%
    filter(year(datetime) == 2023, is.finite(obs), is.finite(mod),
           is.finite(sensor_x), is.finite(sensor_y)) %>%
    group_by(treatment, plot, datetime, sensor_x, sensor_y) %>%
    summarise(obs = mean(obs), mod = mean(mod), .groups = "drop") %>%
    mutate(local_datetime = datetime - hours(5),
           date = as.Date(local_datetime, tz = "UTC"),
           local_hour = hour(local_datetime),
           period = "Day") %>%
    filter(year(local_datetime) == 2023, local_hour %in% local_hours) %>%
    # Each row represents one matched sensor-hour.
    # The >=6 hours/day filter is applied subsequently in temp_hourly.
    ungroup()
}

temp_hourly <- bind_rows(lapply(names(temperature_plots), function(trt) {
  read_temperature(trt, unname(temperature_plots[[trt]]))
})) %>%

  # Require at least 6 valid hours per sensor-day
  group_by(treatment, plot, date, sensor_x, sensor_y) %>%
  filter(n_distinct(datetime) >= min_hours) %>%
  ungroup() %>%

  # Require at least 10 valid sensors per individual hour
  group_by(treatment, plot, datetime) %>%
  filter(n() >= min_sensors) %>%

  # Calculate hourly spatial temperature anomalies
  mutate(
    obs_value = obs - mean(obs),
    mod_value = mod - mean(mod)
  ) %>%
  ungroup() %>%

  transmute(
    treatment, plot, date, datetime, period,
    variable = "Air temperature",
    sensor_id = paste(sensor_x, sensor_y, sep = "_"),
    obs_value, mod_value
  )

# RADIATION: matched hourly PAR (08:00--16:00 Panama LT)
# -> daily mean per sensor
read_radiation <- function(treatment, plot_id) {
  path <- file.path(radiation_dir, paste0("timeseries_", plot_id, "_calibrated.rds"))
  if (!file.exists(path)) stop("Missing file: ", path)
  dat <- readRDS(path)
  required <- c("datetime", "location_ID_number", "radiation_mean", "PAR_model")
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  dat %>%
    transmute(treatment = treatment, plot = plot_id,
              datetime = as.POSIXct(datetime, tz = "UTC"),
              sensor_id = as.character(location_ID_number),
              obs = as.numeric(radiation_mean), mod = as.numeric(PAR_model)) %>%
    mutate(local_datetime = datetime - hours(5),
           date = as.Date(local_datetime, tz = "UTC"),
           local_hour = hour(local_datetime)) %>%
    filter(year(local_datetime) == 2023,
           local_hour %in% local_hours,
           is.finite(obs), is.finite(mod), obs >= 0, mod >= 0) %>%
    group_by(treatment, plot, date, sensor_id) %>%
    summarise(obs = mean(obs), mod = mean(mod),
              n_hours = n_distinct(datetime), .groups = "drop") %>%
    filter(n_hours >= min_hours)
}

rad_daily <- bind_rows(lapply(names(radiation_plots), function(trt) {
  read_radiation(trt, unname(radiation_plots[[trt]]))
})) %>%
  group_by(treatment, plot, date) %>%
  filter(n_distinct(sensor_id) >= min_sensors) %>%
  # Retain original spatial-mean normalisation of daily sensor values.
  mutate(obs_mean = mean(obs), mod_mean = mean(mod)) %>%
  ungroup() %>%
  filter(obs_mean >= min_mean_PAR, mod_mean > 0) %>%
  transmute(treatment, plot, date, datetime = as.POSIXct(NA, tz = "UTC"),
            period = "Day", variable = "Radiation",
            sensor_id, obs_value = obs / obs_mean,
            mod_value = mod / mod_mean)

combined <- bind_rows(temp_hourly, rad_daily) %>%
  mutate(treatment = factor(treatment, levels = c("Control", "Removal")),
         variable = factor(variable, levels = c("Air temperature", "Radiation")),
         period = factor(period, levels = c("Day", "Night")))
if (!nrow(combined)) stop("No valid observations remain.")

# For temperature, compute a separate spatial distribution per HOUR.
# For radiation, retain the daily distribution of per-sensor daily means.
distribution_metrics <- combined %>%
  group_by(treatment, plot, variable, period, date, datetime) %>%
  summarise(n_sensors = n(),
            obs_p10 = q(obs_value, .10), mod_p10 = q(mod_value, .10),
            obs_p50 = q(obs_value, .50), mod_p50 = q(mod_value, .50),
            obs_p90 = q(obs_value, .90), mod_p90 = q(mod_value, .90),
            wasserstein = wasserstein_1d(obs_value, mod_value),
            .groups = "drop") %>%
  mutate(obs_heterogeneity = obs_p90 - obs_p10,
         mod_heterogeneity = mod_p90 - mod_p10,
         heterogeneity_bias = mod_heterogeneity - obs_heterogeneity,
         p10_bias = mod_p10 - obs_p10,
         p50_bias = mod_p50 - obs_p50,
         p90_bias = mod_p90 - obs_p90)

summary_metrics <- distribution_metrics %>%
  group_by(treatment, plot, variable, period) %>%
  summarise(n_days = n_distinct(date),
            n_hours = if (first(as.character(variable)) == "Air temperature") n() else NA_integer_,
            mean_sensors = mean(n_sensors),
            median_wasserstein = median(wasserstein),
            median_observed_P90P10 = median(obs_heterogeneity),
            median_modelled_P90P10 = median(mod_heterogeneity),
            median_heterogeneity_bias = median(heterogeneity_bias),
            median_p10_bias = median(p10_bias),
            median_p50_bias = median(p50_bias),
            median_p90_bias = median(p90_bias), .groups = "drop")

write.csv(distribution_metrics, file.path(output_dir, "distribution_metrics_hourly_temperature.csv"), row.names = FALSE)
write.csv(summary_metrics, file.path(output_dir, "distribution_summary.csv"), row.names = FALSE)
saveRDS(combined, file.path(output_dir, "spatial_distributions_hourly_temperature.rds"))

plot_data <- combined %>%
  pivot_longer(c(obs_value, mod_value), names_to = "source", values_to = "value") %>%
  mutate(source = recode(source, obs_value = "Observed", mod_value = "Modelled"),
         source = factor(source, levels = c("Observed", "Modelled")))

base_theme <- theme_bw(base_size = 17) +
  theme(strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "grey85", colour = "grey35"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.position = "bottom", legend.text = element_text(size = 17),
        panel.grid.minor = element_blank())

# Each variable gets its own x-axis label, but both columns have equal widths.
make_panel <- function(data, variable_name, geometry = "density") {
  d <- data %>% filter(variable == variable_name) %>%
    mutate(facet_title = variable_name)
  if (!nrow(d)) return(NULL)
  if (geometry == "density") {
    p <- ggplot(d, aes(x = value, colour = source, fill = source)) +
      geom_density(alpha = 0.20, linewidth = 1.1,
                   adjust = kde_adjust, na.rm = TRUE) +
      scale_fill_manual(values = colours)
  } else {
    p <- ggplot(d, aes(x = value, colour = source)) +
      stat_ecdf(linewidth = 1.1, na.rm = TRUE)
  }
  p +
    facet_grid(treatment ~ facet_title) +
    scale_colour_manual(values = colours) +
    labs(x = if (variable_name == "Air temperature")
      "Temperature anomaly relative to spatial mean (°C)"
      else "Relative light availability (spatial mean = 1)",
      y = if (variable_name == "Air temperature")
        if (geometry == "density") "Density" else "Cumulative proportion"
      else NULL,
      colour = NULL, fill = NULL) +
    base_theme +
    theme(
      strip.text.y = if (variable_name == "Air temperature")
        element_blank() else element_text(size = 18),
      strip.background.y = if (variable_name == "Air temperature")
        element_blank() else element_rect(fill = "grey85", colour = "grey35")
    ) +
    {if (variable_name == "Radiation")
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      else theme()}
}

day_data <- plot_data %>% filter(period == "Day")
for (geometry in c("density", "ecdf")) {
  p_temp <- make_panel(day_data, "Air temperature", geometry)
  p_rad <- make_panel(day_data, "Radiation", geometry)
  combined_plot <- (p_temp + p_rad +
                      plot_layout(ncol = 2, widths = c(1, 1), guides = "collect")) &
    theme(legend.position = "bottom")
  ggsave(file.path(output_dir, paste0("combined_", geometry, "_day.png")),
         combined_plot, width = 14, height = 8, dpi = 300, bg = "white")
}

print(summary_metrics, width = Inf)
cat("\nSaved hourly-temperature validation to: ", output_dir, "\n", sep = "")
