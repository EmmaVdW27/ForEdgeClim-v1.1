################################################################################
# Canopy structure-class validation: radiation (C8/R7) and air temperature (C13/R11)
# TLS / ALS / combined TLS+ALS, 2023, 08:00-16:00 Panama local time
# Source-specific pooled quintiles within each variable's two plots.
# Temperature columns: SN, value_mean, temperature_model.
#
# Author: Emma Van de Walle - Q-ForestLab
################################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readxl)
  library(RANN); library(lubridate)
})

root <- "."
timeseries_dir <- file.path(root, "Output", "timeseries", "data")
sensor_file <- file.path(root, "Data", "sensor_locations_all_plots.xlsx")
output_dir <- file.path(root, "Output", "validation", "TLS_vs_ALS_uncalibrated_canopy_classes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sources <- c("TLS", "ALS", "Combined")
class_labels <- c("Very closed", "Closed", "Intermediate", "Open", "Very open")
sensor_height <- 1
neighbourhood_radius <- 1L
min_hours <- 6L
min_sensors <- 10L
min_mean_PAR <- 5
# Set TRUE if datetime values are already Panama local wall-clock time.
datetime_is_local <- FALSE
# Choose 'auto', 'C', or 'K'. Auto converts values with median > 150 to Celsius.
temp_obs_units <- "auto"
temp_mod_units <- "auto"

temp_sensor_types <- c("TMS")

check_columns <- function(x, cols, label) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) stop(label, ": missing columns: ", paste(missing, collapse = ", "))
}
normalise_id <- function(x) {
  x <- trimws(as.character(x))
  numeric_id <- grepl("^[0-9]+(\\.0+)?$", x)
  x[numeric_id] <- as.character(as.integer(as.numeric(x[numeric_id])))
  x
}
to_celsius <- function(x, units) {
  x <- as.numeric(x)
  if (units == "K" || (units == "auto" && any(is.finite(x)) &&
                       median(x[is.finite(x)]) > 150)) x <- x - 273.15
  x
}

plot_config <- tibble(
  variable = rep(c("Radiation", "Air temperature"), each = 2),
  plot = c("C8", "R7", "C13", "R11")
)
input_files <- tidyr::crossing(plot_config, source = sources) %>%
  mutate(
    suffix = case_when(source == "TLS" ~ "", source == "ALS" ~ "_ALS",
                       TRUE ~ "_ALS_TLS"),
    ts_path = file.path(timeseries_dir,
                        paste0("timeseries_", plot, "_uncalibrated", suffix, ".rds")),
    voxel_path = file.path(root, "Data", paste0(
      case_when(source == "TLS" ~ "TLS", source == "ALS" ~ "ALS",
                TRUE ~ "ALS_TLS"), "_scaled_DTM_and_grid_", plot, "_Panama.rds"))
  )
for (path in c(sensor_file, input_files$ts_path, input_files$voxel_path)) {
  if (!file.exists(path)) stop("Missing input: ", path)
}

# Sensor metadata. Check the TOMST labels and IDs against the timeseries.
sensors <- readxl::read_excel(sensor_file)
check_columns(sensors, c("x", "y", "SN", "sensor", "plot"), "Sensor file")
sensors <- sensors %>% mutate(
  plot = case_when(
    tolower(trimws(as.character(plot))) == "cplot08" ~ "C8",
    tolower(trimws(as.character(plot))) == "rplot07" ~ "R7",
    tolower(trimws(as.character(plot))) == "cplot13" ~ "C13",
    tolower(trimws(as.character(plot))) == "rplot11" ~ "R11",
    TRUE ~ toupper(trimws(as.character(plot)))
  ),
  sensor_id = normalise_id(SN),
  sensor_type = toupper(trimws(as.character(sensor))),
  x = as.numeric(x), y = as.numeric(y)
) %>% inner_join(plot_config, by = "plot") %>%
  filter((variable == "Radiation" & sensor_type == "LS") |
           (variable == "Air temperature" & sensor_type %in% temp_sensor_types),
         !is.na(sensor_id), sensor_id != "", is.finite(x), is.finite(y)) %>%
  select(variable, plot, sensor_id, x, y)
if (anyDuplicated(sensors[c("plot", "sensor_id")])) stop("Duplicate sensor IDs")
for (pl in plot_config$plot) {
  if (!any(sensors$plot == pl)) stop("No matching sensors for ", pl,
                                     ". Check temp_sensor_types, plot labels and sensor IDs in the Excel file.")
}

calculate_structure <- function(plot_id, source_id, voxel_path) {
  obj <- readRDS(voxel_path)
  grid <- if (is.list(obj) && !is.null(obj$grid)) obj$grid else obj
  check_columns(grid, c("X", "Y", "Z", "density"), voxel_path)
  grid <- grid %>% transmute(X = as.numeric(X), Y = as.numeric(Y),
                             Z = as.numeric(Z), density = as.numeric(density))
  if (any(!is.finite(as.matrix(grid))) || any(grid$density < 0 | grid$density > 1))
    stop("Invalid voxel values: ", voxel_path)
  if (anyDuplicated(grid[c("X", "Y", "Z")])) stop("Duplicate voxel coordinates: ", voxel_path)
  xs <- sort(unique(grid$X)); ys <- sort(unique(grid$Y))
  dx <- median(diff(xs)); dy <- median(diff(ys))
  if (length(xs) < 3 || length(ys) < 3 ||
      any(abs(diff(xs) - dx) > abs(dx) * 1e-5) ||
      any(abs(diff(ys) - dy) > abs(dy) * 1e-5)) stop("Irregular voxel XY grid: ", voxel_path)
  grid <- grid %>% mutate(ix = match(X, xs), iy = match(Y, ys))
  grid_xy <- grid %>% distinct(X, Y, ix, iy)
  ss <- sensors %>% filter(plot == plot_id)
  nn <- RANN::nn2(as.matrix(grid_xy[, c("X", "Y")]),
                  query = as.matrix(ss[, c("x", "y")]), k = 1)
  ss$ix <- grid_xy$ix[nn$nn.idx[, 1]]
  ss$iy <- grid_xy$iy[nn$nn.idx[, 1]]
  ss$match_distance_m <- as.numeric(nn$nn.dists[, 1])
  bind_rows(lapply(seq_len(nrow(ss)), function(i) {
    s <- ss[i, ]
    ng <- grid %>% filter(abs(ix - s$ix) <= neighbourhood_radius,
                          abs(iy - s$iy) <= neighbourhood_radius, Z > sensor_height)
    if (!nrow(ng)) return(NULL)
    ncol <- n_distinct(paste(ng$ix, ng$iy)); nlevels <- n_distinct(ng$Z)
    if (nrow(ng) != ncol * nlevels)
      stop("Incomplete voxel neighbourhood: ", plot_id, " sensor ", s$sensor_id)
    tibble(variable = s$variable, plot = plot_id, source = source_id,
           sensor_id = s$sensor_id, openness = 1 - mean(ng$density),
           n_columns = ncol, match_distance_m = s$match_distance_m)
  }))
}
structure <- bind_rows(lapply(seq_len(nrow(input_files)), function(i) {
  f <- input_files[i, ]
  calculate_structure(f$plot, f$source, f$voxel_path)
}))
if (any(structure$n_columns < 9)) warning("Some sensors have <9 voxel columns")
if (any(structure$match_distance_m > 2)) warning("Some sensors >2 m from nearest voxel column")

read_timeseries <- function(variable_id, plot_id, source_id, path) {
  x <- readRDS(path)
  if (is.list(x) && !is.data.frame(x) && !is.null(x$data)) x <- x$data
  check_columns(x, "datetime", path)
  if (variable_id == "Radiation") {
    check_columns(x, c("location_ID_number", "radiation_mean", "PAR_model"), path)
    ids <- x$location_ID_number
    obs <- as.numeric(x$radiation_mean); mod <- as.numeric(x$PAR_model)
  } else {
    check_columns(x, c("SN", "value_mean", "temperature_model"), path)
    ids <- x$SN
    obs <- to_celsius(x$value_mean, temp_obs_units)
    mod <- to_celsius(x$temperature_model, temp_mod_units)
  }
  timestamp <- if (datetime_is_local) {
    lubridate::force_tz(as.POSIXct(x$datetime, tz = "UTC"), "America/Panama")
  } else {
    lubridate::with_tz(as.POSIXct(x$datetime, tz = "UTC"), "America/Panama")
  }
  tibble(variable = variable_id, plot = plot_id, source = source_id,
         sensor_id = normalise_id(ids),
         datetime_local = timestamp,
         date = as.Date(timestamp, tz = "America/Panama"),
         hour = lubridate::hour(timestamp), obs = obs, mod = mod) %>%
    filter(lubridate::year(date) == 2023, hour >= 8, hour <= 16,
           is.finite(obs), is.finite(mod)) %>%
    { if (variable_id == "Radiation") filter(., obs >= 0, mod >= 0) else . } %>%
    group_by(variable, plot, source, sensor_id, date, hour) %>%
    summarise(obs = mean(obs), mod = mean(mod), .groups = "drop") %>%
    group_by(variable, plot, source, sensor_id, date) %>%
    summarise(obs = median(obs), mod = median(mod), n_hours = n(), .groups = "drop") %>%
    filter(n_hours >= min_hours)
}
daily_raw <- bind_rows(lapply(seq_len(nrow(input_files)), function(i) {
  f <- input_files[i, ]
  message("Reading ", basename(f$ts_path))
  read_timeseries(f$variable, f$plot, f$source, f$ts_path)
}))
if (!nrow(daily_raw)) stop("No valid daily data")

# Pair only sensor-days shared by all three source simulations.
key <- c("variable", "plot", "sensor_id", "date")
paired <- daily_raw %>%
  select(all_of(key), source, obs, mod) %>%
  pivot_wider(names_from = source, values_from = c(obs, mod), names_sep = "_") %>%
  filter(if_all(all_of(c("obs_TLS", "obs_ALS", "obs_Combined",
                         "mod_TLS", "mod_ALS", "mod_Combined")), ~ !is.na(.x))) %>%
  inner_join(structure %>% count(variable, plot, sensor_id) %>%
               filter(n == 3L) %>% select(-n),
             by = c("variable", "plot", "sensor_id"))
if (!nrow(paired)) stop("No matched sensor-days. Check sensor IDs against timeseries.")
if (any(abs(paired$obs_TLS - paired$obs_ALS) > 1e-5 |
        abs(paired$obs_TLS - paired$obs_Combined) > 1e-5))
  warning("Observed values differ across source files; TLS observations used")
paired <- paired %>% group_by(variable, plot, date) %>%
  filter(n_distinct(sensor_id) >= min_sensors) %>%
  group_modify(~ {
    x <- .x
    if (.y$variable == "Radiation") {
      if (mean(x$obs_TLS) < min_mean_PAR ||
          any(c(mean(x$mod_TLS), mean(x$mod_ALS), mean(x$mod_Combined)) <= 0))
        return(x[0, ])
      x <- x %>% mutate(obs_value = obs_TLS / mean(obs_TLS),
                        TLS_value = mod_TLS / mean(mod_TLS),
                        ALS_value = mod_ALS / mean(mod_ALS),
                        Combined_value = mod_Combined / mean(mod_Combined))
    } else {
      x <- x %>% mutate(obs_value = obs_TLS - mean(obs_TLS),
                        TLS_value = mod_TLS - mean(mod_TLS),
                        ALS_value = mod_ALS - mean(mod_ALS),
                        Combined_value = mod_Combined - mean(mod_Combined))
    }
    x
  }) %>% ungroup()
if (!nrow(paired)) stop("No matched days after filters")

# Source-specific pooled quintiles, separately for C8/R7 and C13/R11.
structure_three <- structure %>%
  semi_join(paired %>% distinct(variable, plot, sensor_id),
            by = c("variable", "plot", "sensor_id")) %>%
  group_by(variable, source) %>%
  group_modify(~ {
    bounds <- as.numeric(quantile(.x$openness, probs = seq(0, 1, 0.2), na.rm = TRUE))
    if (anyDuplicated(bounds)) stop("Non-unique quintile boundaries: ",
                                    .y$variable, " / ", .y$source)
    .x %>% mutate(canopy_class = factor(cut(openness, breaks = bounds,
                                            labels = class_labels, include.lowest = TRUE), levels = class_labels))
  }) %>% ungroup()
class_table <- structure_three %>% group_by(variable, source) %>%
  group_modify(~ {
    b <- as.numeric(quantile(.x$openness, probs = seq(0, 1, 0.2), na.rm = TRUE))
    tibble(canopy_class = class_labels, lower = head(b, -1), upper = tail(b, -1))
  }) %>% ungroup()
paired_long <- bind_rows(lapply(sources, function(src) {
  paired %>% transmute(variable, plot, sensor_id, date, source = src,
                       Observed = obs_value, Modelled = .data[[paste0(src, "_value")]])
})) %>% left_join(structure_three %>% select(variable, plot, sensor_id, source,
                                             openness, canopy_class),
                  by = c("variable", "plot", "sensor_id", "source"))
if (anyNA(paired_long$canopy_class)) stop("Missing structure class after join")
annual <- paired_long %>% group_by(variable, plot, source, sensor_id,
                                   openness, canopy_class) %>%
  summarise(Observed = median(Observed), Modelled = median(Modelled),
            n_days = n_distinct(date), .groups = "drop")
class_summary <- annual %>% pivot_longer(c(Observed, Modelled),
                                         names_to = "series", values_to = "value") %>%
  group_by(variable, plot, source, canopy_class, series) %>%
  summarise(median_value = median(value), q25 = quantile(value, 0.25),
            q75 = quantile(value, 0.75), n_sensors = n(), .groups = "drop")

# Separate 2 x 3 figures.
colours <- c("Observed" = "blue", "Modelled" = "red")
for (variable_id in c("Radiation", "Air temperature")) {
  these_plots <- plot_config$plot[plot_config$variable == variable_id]
  fig_data <- annual %>%
    filter(variable == variable_id) %>%
    pivot_longer(
      c(Observed, Modelled),
      names_to = "series",
      values_to = "value"
    ) %>%
    mutate(
      plot = factor(plot, levels = these_plots),
      structure_input = factor(
        source,
        levels = sources,
        labels = c("TLS", "ALS", "TLS + ALS")
      ),
      series = factor(series, levels = c("Observed", "Modelled"))
    )

  p <- ggplot(
    fig_data,
    aes(x = canopy_class, y = value,
        fill = series, colour = series)
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.8),
      width = 0.65,
      alpha = 0.20,
      linewidth = 0.8,
      outlier.size = 1.5
    ) +
    facet_grid(
      plot ~ structure_input,
      labeller = labeller(
        plot = c(
          "C8"  = "Control",
          "C13" = "Control",
          "R7"  = "Removal",
          "R11" = "Removal"
        )
      )
    ) +
    scale_fill_manual(values = colours) + scale_colour_manual(values = colours) +
    scale_x_discrete(drop = FALSE, labels = c("Very\nclosed", "Closed",
                                              "Intermediate", "Open", "Very\nopen")) +
    labs(x = "Canopy structure class (source-specific pooled quintiles)",
         y = if (variable_id == "Radiation") "Relative light availability" else
           "Air temperature anomaly (°C)", fill = NULL, colour = NULL) +
    guides(colour = "none", fill = guide_legend(
      override.aes = list(alpha = 0.20, colour = c("blue", "red")))) +
    theme_bw(base_size = 17) +
    theme(strip.text = element_text(size = 23),
          axis.title = element_text(size = 24),
          axis.text.x = element_text(size = 18, colour = "black"),
          axis.text.y = element_text(size = 19, colour = "black"),
          legend.position = "bottom", legend.text = element_text(size = 21),
          panel.grid.minor = element_blank(),
          panel.spacing = grid::unit(0.9, "cm"),
          plot.margin = margin(20, 25, 20, 20))
  # Temperature anomalies can be positive or negative.
  if (variable_id == "Air temperature") p <- p + geom_hline(yintercept = 0, linewidth = 0.4)
  stem <- if (variable_id == "Radiation") "radiation" else "air_temperature"
  ggsave(file.path(output_dir, paste0("TLS_ALS_combined_", stem, ".png")),
         p, width = 23, height = 10, dpi = 300)
  ggsave(file.path(output_dir, paste0("TLS_ALS_combined_", stem, ".pdf")),
         p, width = 23, height = 10)
  print(p)
}
write.csv(paired, file.path(output_dir, "paired_daily_radiation_temperature.csv"), row.names = FALSE)
write.csv(paired_long, file.path(output_dir, "paired_daily_three_classes_both.csv"), row.names = FALSE)
write.csv(structure_three, file.path(output_dir, "sensor_structure_both.csv"), row.names = FALSE)
write.csv(class_table, file.path(output_dir, "pooled_class_boundaries_both.csv"), row.names = FALSE)
write.csv(annual, file.path(output_dir, "annual_sensor_values_both.csv"), row.names = FALSE)
write.csv(class_summary, file.path(output_dir, "class_summary_both.csv"), row.names = FALSE)
write.csv(structure_three %>% count(variable, plot, source, canopy_class,
                                    name = "n_sensors", .drop = FALSE),
          file.path(output_dir, "class_sensor_counts_both.csv"), row.names = FALSE)
cat("\nPooled class boundaries:\n"); print(class_table)
cat("\nValid matched days:\n"); print(paired %>% distinct(variable, plot, date) %>% count(variable, plot))
cat("\nCompleted: ", output_dir, "\n", sep = "")
