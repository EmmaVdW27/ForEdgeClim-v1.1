###############################################################################
#
# ForEdgeClim - Liana effects on forest microclimate
#
# Quantifies differences in near-surface microclimate and energy-balance
# variables between liana-removal and control plots throughout 2023.
#
# For each plot, daily daytime and nighttime medians are first calculated from
# the hourly data. Treatment-level daily values are then obtained by averaging
# equally across the four plots. Only dates with data from all four control and
# all four removal plots are retained.
#
# Treatment effects are expressed as:
#
#     Liana removal - Control
#
# for observed and modelled air temperature, surface temperature, latent heat
# flux, relative humidity, vapour pressure deficit, and direct and diffuse
# shortwave radiation.
#
# Daytime is defined as 11:00-22:00 UTC (06:00-17:00 local time, UTC-5);
# all remaining hours are classified as nighttime.
#
# Wet season: 1 May - 15 December
# Dry season: 16 December - 30 April
#
# Exact plot-level bootstrap confidence intervals are calculated by enumerating
# all unique bootstrap compositions of the four plots within each treatment
# and weighting them by their bootstrap probabilities.
#
# Author: Emma Van de Walle - Q-ForestLab
#
###############################################################################

library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(patchwork)
library(grid)


#########
# INPUT #
#########

temperature_from_model = "Tair_mean_1m" # Tair_mean_1m or temperature_model (only sensor voxels)

output_path <- "Output/liana-effect/"

control <- readRDS("Output/timeseries/data/timeseries_control_calibrated.rds")
removal <- readRDS("Output/timeseries/data/timeseries_removal_calibrated.rds")

###############################################################################
# CLEAN DATETIME
###############################################################################

clean_datetime <- function(df) {

  df %>%
    mutate(
      datetime = trimws(datetime),
      datetime = na_if(datetime, ""),
      datetime = na_if(datetime, "NA"),
      datetime = suppressWarnings(ymd_hms(datetime))
    ) %>%
    filter(!is.na(datetime))

}

###############################################################################
# DAY/NIGHT CLASSIFICATION (Panama, UTC)
###############################################################################

add_period <- function(df) {

  df %>%
    mutate(
      hour_utc = hour(datetime),

      period = case_when(
        hour_utc >= 11 & hour_utc < 23 ~ "Day",
        TRUE                           ~ "Night"
      )
    )

}

###############################################################################
# AGGREGATION
###############################################################################

control <- control %>%
  clean_datetime() %>%
  add_period() %>%
  mutate(date = as.Date(datetime))

removal <- removal %>%
  clean_datetime() %>%
  add_period() %>%
  mutate(date = as.Date(datetime))

control_summary <-
  control %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(plot, date, period) %>% # daily median per plot (medians to exclude outliers)
  summarise(
    obs_control    = median(value_mean, na.rm = TRUE),
    mod_control    = median(.data[[temperature_from_model]], na.rm = TRUE),
    LE_control     = median(LE_mean, na.rm = TRUE),
    RH_control     = median(RH_mean, na.rm = TRUE),
    VPD_control    = median(VPD_mean, na.rm = TRUE),
    Fd_control     = median(Fd_down_mean, na.rm = TRUE),
    Fb_control     = median(Fb_down_mean, na.rm = TRUE),
    Tsurf_control  = median(Tsurf_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(date, period) %>% # mean over plots
  summarise(
    n_control      = n(),
    obs_control    = mean(obs_control, na.rm = TRUE),
    mod_control    = mean(mod_control, na.rm = TRUE),
    LE_control     = mean(LE_control, na.rm = TRUE),
    RH_control     = mean(RH_control, na.rm = TRUE),
    VPD_control    = mean(VPD_control, na.rm = TRUE),
    Fd_control     = mean(Fd_control, na.rm = TRUE),
    Fb_control     = mean(Fb_control, na.rm = TRUE),
    Tsurf_control  = mean(Tsurf_control, na.rm = TRUE),
    .groups = "drop"
  )

removal_summary <-
  removal %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(plot, date, period) %>%
  summarise(
    obs_removal   = median(value_mean, na.rm = TRUE),
    mod_removal   = median(.data[[temperature_from_model]], na.rm = TRUE),
    LE_removal    = median(LE_mean, na.rm = TRUE),
    RH_removal    = median(RH_mean, na.rm = TRUE),
    VPD_removal   = median(VPD_mean, na.rm = TRUE),
    Fd_removal    = median(Fd_down_mean, na.rm = TRUE),
    Fb_removal    = median(Fb_down_mean, na.rm = TRUE),
    Tsurf_removal = median(Tsurf_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(date, period) %>%
  summarise(
    n_removal     = n(),
    obs_removal   = mean(obs_removal, na.rm = TRUE),
    mod_removal   = mean(mod_removal, na.rm = TRUE),
    LE_removal    = mean(LE_removal, na.rm = TRUE),
    RH_removal    = mean(RH_removal, na.rm = TRUE),
    VPD_removal   = mean(VPD_removal, na.rm = TRUE),
    Fd_removal    = mean(Fd_removal, na.rm = TRUE),
    Fb_removal    = mean(Fb_removal, na.rm = TRUE),
    Tsurf_removal = mean(Tsurf_removal, na.rm = TRUE),
    .groups = "drop"
  )

###############################################################################
# JOIN
###############################################################################

comparison <-
  full_join(
    control_summary,
    removal_summary,
    by = c(
      "date",
      "period"
    )
  )%>%
  filter(
    n_control == 4, # only use days on which we have data from the 4 control plots
    n_removal == 4 # same for removal plots
  )


###############################################################################
# DELTAS
###############################################################################

comparison <-
  comparison %>%
  mutate(

    delta_obs =
      obs_removal -
      obs_control,

    delta_mod =
      mod_removal -
      mod_control,

    delta_LE =
      LE_removal -
      LE_control,

    delta_RH =
      RH_removal -
      RH_control,

    delta_VPD =
      VPD_removal -
      VPD_control,

    delta_Fd =
      Fd_removal -
      Fd_control,

    delta_Fb =
      Fb_removal -
      Fb_control,

    delta_Tsurf =
      Tsurf_removal -
      Tsurf_control
  )


###############################################################################
# SEASON CLASSIFICATION
###############################################################################

comparison <- comparison %>%
  mutate(
    season = ifelse(
      month(date) >= 5 &
        (
          month(date) < 12 |
            (month(date) == 12 & day(date) <= 15)
        ),
      "Wet season",
      "Dry season"
    )
  )

###############################################################################
# WET SEASON BACKGROUND
###############################################################################

years <- unique(year(comparison$date))

wet_season <- data.frame(
  xmin = as.Date(paste0(years, "-05-01")),
  xmax = as.Date(paste0(years, "-12-15")),
  ymin = -Inf,
  ymax = Inf
)

###############################################################################
# SEASON RANGES
###############################################################################

years <- sort(unique(year(comparison$date)))

wet_ranges <- data.frame(
  season = "Wet season",
  xmin = as.Date(paste0(years, "-05-01")),
  xmax = as.Date(paste0(years, "-12-15"))
)

dry_ranges_1 <- data.frame(
  season = "Dry season",
  xmin = as.Date(paste0(years, "-12-16")),
  xmax = as.Date(paste0(years, "-12-31"))
)

dry_ranges_2 <- data.frame(
  season = "Dry season",
  xmin = as.Date(paste0(years, "-01-01")),
  xmax = as.Date(paste0(years, "-04-30"))
)

season_ranges <- bind_rows(
  wet_ranges,
  dry_ranges_1,
  dry_ranges_2
) %>%
  filter(
    xmin <= max(comparison$date) &
      xmax >= min(comparison$date)
  )

###############################################################################
# ABSOLUTE PROCESS VARIABLES FOR sw AND LE
###############################################################################

absolute_process_values <- comparison %>%
  mutate(
    SW_control = Fd_control + Fb_control,
    SW_removal = Fd_removal + Fb_removal
  ) %>%
  group_by(season, period) %>%
  summarise(
    SW_control = median(SW_control, na.rm = TRUE),
    SW_removal = median(SW_removal, na.rm = TRUE),

    LE_control = median(LE_control, na.rm = TRUE),
    LE_removal = median(LE_removal, na.rm = TRUE),

    .groups = "drop"
  )

print(absolute_process_values)

###############################################################################
# MACRO AIR TEMPERATURE
###############################################################################

macro_data <- readRDS("Data/TMS_macro_2023_4C_4R.rds") %>%
  mutate(
    plot = sub("^([A-Z])0+", "\\1", plot),
    datetime = as.POSIXct(datetime, tz = "UTC"),
    datetime_UTC_hour = floor_date(datetime, unit = "hour")
  ) %>%
  filter(year(datetime_UTC_hour) == 2023) %>%
  group_by(plot, datetime_UTC_hour) %>%
  summarise(
    temp_macro = mean(temp_macro, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    date = as.Date(datetime_UTC_hour),
    hour_utc = hour(datetime_UTC_hour),
    period = if_else(
      hour_utc >= 11 & hour_utc < 23,
      "Day",
      "Night"
    ),
    treatment = if_else(
      grepl("^C", plot),
      "Control",
      "Removal"
    )
  )

macro_summary <- macro_data %>%
  group_by(treatment, plot, date, period) %>%
  summarise(
    macro_temp = median(temp_macro, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(treatment, date, period) %>%
  summarise(
    n_plots = n_distinct(plot),
    macro_temp = mean(macro_temp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_plots == 4)

###############################################################################
# ABSOLUTE AIR TEMPERATURE
###############################################################################

absolute_air <- bind_rows(

  comparison %>%
    transmute(
      date,
      season,
      period,
      treatment = "Control",
      type = "Observed",
      value = obs_control
    ),

  comparison %>%
    transmute(
      date,
      season,
      period,
      treatment = "Control",
      type = "Modelled",
      value = mod_control
    ),

  comparison %>%
    transmute(
      date,
      season,
      period,
      treatment = "Removal",
      type = "Observed",
      value = obs_removal
    ),

  comparison %>%
    transmute(
      date,
      season,
      period,
      treatment = "Removal",
      type = "Modelled",
      value = mod_removal
    ),

  macro_summary %>%
    inner_join(
      comparison %>% distinct(date, period, season),
      by = c("date", "period")
    ) %>%
    transmute(
      date,
      season,
      period,
      treatment,
      type = "Macroclimate",
      value = macro_temp
    )

)


###############################################################################
# PROCESS VARIABLES LONG FORMAT
###############################################################################

comparison_long <- comparison %>%
  select(
    date,
    season,
    period,
    delta_LE,
    delta_RH,
    delta_VPD,
    delta_Fd,
    delta_Fb,
    delta_Tsurf
  ) %>%
  pivot_longer(
    cols = -c(date, season, period),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(type = "Modelled")

air_long <-
  comparison %>%
  transmute(
    date,
    season,
    period,
    variable = "c) Δ Air T (°C)",
    Observed = delta_obs,
    Modelled = delta_mod
  ) %>%
  pivot_longer(
    c(Observed, Modelled),
    names_to = "type",
    values_to = "value"
  )

comparison_long <- bind_rows(
  comparison_long,
  air_long
)

comparison_long_day <- comparison_long %>%
  mutate(
    variable = recode(
      variable,
      delta_Fd        = "a) Δ Diffuse SW (W m⁻²)",
      delta_Fb        = "b) Δ Direct SW (W m⁻²)",
      delta_Tsurf     = "c) Δ Surface T (°C)",
      "c) Δ Air T (°C)" = "d) Δ Air T (°C)",
      delta_LE        = "e) Δ LE (W m⁻²)",
      delta_RH        = "f) Δ RH (%)",
      delta_VPD       = "g) Δ VPD (kPa)"
    )
  )

comparison_long_night <- comparison_long %>%
  mutate(
    variable = recode(
      variable,
      delta_Tsurf     = "a) Δ Surface T (°C)",
      "c) Δ Air T (°C)" = "b) Δ Air T (°C)",
      delta_LE        = "c) Δ LE (W m⁻²)",
      delta_RH        = "d) Δ RH (%)",
      delta_VPD       = "e) Δ VPD (kPa)"
    )
  )


###############################################################################
# FACET FUNCTION
###############################################################################

plot_facet <- function(df, vars, period_filter, title = NULL, show_legend = FALSE) {

  df_sub <- df %>%
    filter(
      variable %in% vars,
      period == period_filter
    )

  means <- df_sub %>%
    group_by(variable, season, type) %>%
    summarise(
      mean_val = median(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    merge(season_ranges, by = "season") %>%
    mutate(
      x_pos = ifelse(
        season == "Wet season",
        as.Date("2023-08-15"),
        as.Date("2023-02-15")
      )
    ) %>%
    group_by(variable) %>%
    mutate(
      range_y = diff(range(df_sub$value[df_sub$variable == first(variable)],
                           na.rm = TRUE)),
      y_offset = case_when(
        variable == "d) Δ Air T (°C)" & type == "Observed" ~  0.04 * range_y,
        variable == "d) Δ Air T (°C)" & type == "Modelled" ~ -0.04 * range_y,
        TRUE ~ 0
      ),
      label_y = mean_val + y_offset
    ) %>%
    ungroup()

  p <- ggplot(
    df_sub,
    aes(
      x = date,
      y = value
    )
  ) +

    geom_rect(
      data = wet_season,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = "lightblue",
      alpha = 0.6
    ) +

    geom_smooth(
      aes(
        colour = type,
        fill = type
      ),
      method = "gam",
      formula = y ~ s(x, k = 20),
      size = 1,
      alpha = 0.25,
      na.rm = TRUE,
      fullrange = TRUE
    ) +

    geom_segment(
      data = means,
      aes(
        x = xmin,
        xend = xmax,
        y = mean_val,
        yend = mean_val,
        colour = type
      ),
      linewidth = 0.8,
      show.legend = FALSE
    ) +

    geom_label(
      data = means,
      aes(
        x = x_pos,
        y = label_y,
        label = paste0(
          ifelse(season == "Wet season", "Wet", "Dry"),
          " = ",
          round(mean_val, 2)
        ),
        colour = type
      ),
      fill = "white",
      label.size = 0.3,
      label.r = unit(0.15, "lines"),
      label.padding = unit(0.15, "lines"),
      size = 5,
      fontface = "bold",
      show.legend = FALSE
    ) +

    scale_x_date(
      limits = c(
        as.Date("2023-01-01"),
        as.Date("2023-12-31")
      )
    ) +

    facet_wrap(
      ~variable,
      scales = "free_y",
      ncol = 1
    ) +

    theme_bw() +

    labs(
      title = title,
      x = "Date",
      y = NULL
    ) +

    theme(
      strip.text = element_text(
        size = 18,
        face = "bold"
      ),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      plot.title = element_text(size = 18, face = "bold"),
      legend.position = if(show_legend) "bottom" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      legend.key.width = unit(2, "cm")
    ) +
    scale_colour_manual(
      values = c(
        Observed = "blue",
        Modelled = "red"
      )
    ) +

    scale_fill_manual(
      values = c(
        Observed = "blue",
        Modelled = "red"
      )
    )

  return(p)

}

###############################################################################
# DAY FACETS
###############################################################################

p_radiation_day <- plot_facet(
  comparison_long_day,
  vars = c(
    "a) Δ Diffuse SW (W m⁻²)",
    "b) Δ Direct SW (W m⁻²)"
  ),
  period_filter = "Day",
  title = "Daytime (06-18 EST)",
  show_legend = FALSE
)

p_middle_day <- plot_facet(
  comparison_long_day,
  vars = c(
    "c) Δ Surface T (°C)",
    "d) Δ Air T (°C)"
  ),
  period_filter = "Day",
  show_legend = TRUE
)

p_right_day <- plot_facet(
  comparison_long_day,
  vars = c(
    "e) Δ LE (W m⁻²)",
    "f) Δ RH (%)",
    "g) Δ VPD (kPa)"
  ),
  period_filter = "Day",
  show_legend = FALSE
)

combined_day <-
  p_radiation_day +
  p_middle_day +
  p_right_day +
  plot_layout(ncol = 3)

###############################################################################
# NIGHT FACETS
###############################################################################

p_middle_night <- plot_facet(
  comparison_long_night,
  vars = c(
    "a) Δ Surface T (°C)",
    "b) Δ Air T (°C)"
  ),
  period_filter = "Night",
  title = "Nighttime (18-06 EST)",
  show_legend = TRUE
)

p_right_night <- plot_facet(
  comparison_long_night,
  vars = c(
    "c) Δ LE (W m⁻²)",
    "d) Δ RH (%)",
    "e) Δ VPD (kPa)"
  ),
  period_filter = "Night",
  show_legend = FALSE
)

combined_night <-
  p_middle_night +
  p_right_night +
  plot_layout(ncol = 2)



###############################################################################
# absolute air plot
###############################################################################

plot_absolute_air <- function(df, period_filter, title = NULL, show_legend = FALSE){

  df_sub <- df %>%
    filter(period == period_filter)

  df_sub <- df_sub %>%
    mutate(
      type = factor(
        type,
        levels = c("Observed", "Modelled", "Macroclimate")
      )
    )

  means <- df_sub %>%
    group_by(treatment, season, type) %>%
    summarise(
      mean_val = median(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      xmin = if_else(
        season == "Wet season",
        as.Date("2023-05-01"),
        as.Date("2023-01-01")
      ),
      xmax = if_else(
        season == "Wet season",
        as.Date("2023-12-15"),
        as.Date("2023-04-30")
      ),
      x_pos = case_when(
        season == "Wet season" & type == "Macroclimate" ~ as.Date("2023-06-15"),
        season == "Wet season" & type == "Observed"     ~ as.Date("2023-08-15"),
        season == "Wet season" & type == "Modelled"     ~ as.Date("2023-10-15"),

        season == "Dry season" & type == "Macroclimate" ~ as.Date("2023-01-20"),
        season == "Dry season" & type == "Observed"     ~ as.Date("2023-02-25"),
        season == "Dry season" & type == "Modelled"     ~ as.Date("2023-04-05")
      )
    ) %>%
    group_by(treatment, season) %>%
    mutate(
      x_label = if (n() == 2 && abs(diff(mean_val)) < 1) {
        c(x_pos[1] - 20, x_pos[2] + 20)
      } else {
        x_pos
      }
    ) %>%
    ungroup()

  ggplot(
    df_sub,
    aes(
      x = date,
      y = value
    )
  ) +

    geom_rect(
      data = wet_season,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = "lightblue",
      alpha = 0.6
    ) +

    geom_smooth(
      aes(
        colour = type,
        fill = type
      ),
      method = "gam",
      formula = y ~ s(x, k = 30),
      linewidth = 1,
      alpha = 0.25,
      fullrange = TRUE
    ) +

    geom_segment(
      data = means,
      aes(
        x = xmin,
        xend = xmax,
        y = mean_val,
        yend = mean_val,
        colour = type
      ),
      linewidth = 0.8,
      show.legend = FALSE
    ) +

    geom_label(
      data = means,
      aes(
        x = x_pos,
        y = mean_val,
        label = round(mean_val, 1),
        colour = type
      ),
      fill = "white",
      size = 4,
      show.legend = FALSE
    ) +

    facet_wrap(
      ~treatment,
      ncol = 2,
      scales = "fixed"
    ) +

    scale_colour_manual(
      values = c(
        Observed = "blue",
        Modelled = "red",
        Macroclimate = "black"
      ),
      breaks = c("Observed", "Modelled", "Macroclimate"),
      limits = c("Observed", "Modelled", "Macroclimate"),
      drop = FALSE,
      name = NULL
    ) +

    scale_fill_manual(
      values = c(
        Observed = "blue",
        Modelled = "red",
        Macroclimate = "black"
      ),
      breaks = c("Observed", "Modelled", "Macroclimate"),
      limits = c("Observed", "Modelled", "Macroclimate"),
      drop = FALSE,
      name = NULL
    ) +

    labs(
      title = title,
      x = "Date",
      y = "Air temperature (°C)",
      colour = NULL
    ) +

    coord_cartesian(ylim = c(22, 29)) +

    scale_y_continuous(
      breaks = 22:29,
      labels = function(x) sprintf("%.0f", x)
    ) +

    theme_bw() +

    theme(
      strip.text = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      plot.title = element_text(
        size = 18,
        face = "bold"
      ),
      legend.position = if(show_legend) "bottom" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      legend.key.width = unit(2, "cm"),
      panel.spacing.x = unit(2, "lines")
    ) +
    scale_x_date(
      limits = c(
        as.Date("2023-01-01"),
        as.Date("2023-12-31")
      )
    )
}

###############################################################################
# CREATE FIGURES
###############################################################################

absolute_day <- plot_absolute_air(
  absolute_air,
  period_filter = "Day",
  title = "Daytime (06–18 EST)",
  show_legend = FALSE
)

absolute_night <- plot_absolute_air(
  absolute_air,
  period_filter = "Night",
  title = "Nighttime (18–06 EST)",
  show_legend = TRUE

)

absolute_combined <-
  wrap_plots(
    absolute_day,
    absolute_night,
    ncol = 1
  )

###############################################################################
# SAVE FIGURES
###############################################################################

ggsave(
  paste0(
    output_path,
    "combined_process_facets_day.png"
  ),
  combined_day,
  width = 14,
  height = 8,
  dpi = 300
)

ggsave(
  paste0(
    output_path,
    "combined_process_facets_night.png"
  ),
  combined_night,
  width = 9.33,
  height = 8,
  dpi = 300
)


ggsave(
  paste0(
    output_path,
    "absolute_air_temperature_combined.png"
  ),
  plot = absolute_combined,
  width = 12,
  height = 10,
  dpi = 300
)

###############################################################################
# METRICS
###############################################################################

metrics <- comparison %>%
  group_by(period) %>%
  summarise(
    rmse_v1.1 = sqrt(
      mean(
        (delta_mod - delta_obs)^2,
        na.rm = TRUE
      )
    ),

    correlation_v1.1 = cor(
      delta_obs,
      delta_mod,
      use = "complete.obs"
    ),
    .groups = "drop"
  )

print(metrics)


###############################################################################
# EXACT PLOT-LEVEL BOOTSTRAP
#
# 4 plots per treatment
# 35 unique bootstrap compositions per treatment
# 35 x 35 = 1225 Control x Removal combinations
#
# Each unique composition is weighted according to its probability
# under an ordinary bootstrap with replacement.
###############################################################################


###############################################################################
# 1. DAILY MEDIANS PER INDIVIDUAL PLOT
###############################################################################

control_plot_daily <- control %>%
  group_by(plot, date, period) %>%
  summarise(
    obs   = median(value_mean, na.rm = TRUE),
    mod   = median(.data[[temperature_from_model]], na.rm = TRUE),
    LE    = median(LE_mean, na.rm = TRUE),
    RH    = median(RH_mean, na.rm = TRUE),
    VPD   = median(VPD_mean, na.rm = TRUE),
    Fd    = median(Fd_down_mean, na.rm = TRUE),
    Fb    = median(Fb_down_mean, na.rm = TRUE),
    Tsurf = median(Tsurf_mean, na.rm = TRUE),
    .groups = "drop"
  )


removal_plot_daily <- removal %>%
  group_by(plot, date, period) %>%
  summarise(
    obs   = median(value_mean, na.rm = TRUE),
    mod   = median(.data[[temperature_from_model]], na.rm = TRUE),
    LE    = median(LE_mean, na.rm = TRUE),
    RH    = median(RH_mean, na.rm = TRUE),
    VPD   = median(VPD_mean, na.rm = TRUE),
    Fd    = median(Fd_down_mean, na.rm = TRUE),
    Fb    = median(Fb_down_mean, na.rm = TRUE),
    Tsurf = median(Tsurf_mean, na.rm = TRUE),
    .groups = "drop"
  )


###############################################################################
# 2. FUNCTION TO GENERATE ALL UNIQUE BOOTSTRAP COMPOSITIONS
#
# For four plots, this generates all integer vectors:
#
# (n1, n2, n3, n4)
#
# for which:
#
# n1 + n2 + n3 + n4 = 4
#
# Examples:
# 4 0 0 0
# 3 1 0 0
# 2 2 0 0
# 2 1 1 0
# 1 1 1 1
#
# Total = 35 combinations
###############################################################################

make_boot_compositions <- function(n_plots) {

  grid <- expand.grid(
    replicate(
      n_plots,
      0:n_plots,
      simplify = FALSE
    )
  )

  grid <- grid[
    rowSums(grid) == n_plots,
    ,
    drop = FALSE
  ]

  names(grid) <- paste0(
    "n",
    seq_len(n_plots)
  )


  ###########################################################################
  # BOOTSTRAP PROBABILITY
  #
  # Multinomial probability:
  #
  # n! / (n1! n2! ... nk!) * (1/k)^n
  ###########################################################################

  counts <- as.matrix(grid)

  multiplicity <- apply(
    counts,
    1,
    function(x) {
      factorial(n_plots) /
        prod(factorial(x))
    }
  )

  probability <-
    multiplicity *
    (1 / n_plots)^n_plots


  grid$weight <- probability

  grid
}


###############################################################################
# 3. CREATE THE 35 CONTROL AND 35 REMOVAL COMPOSITIONS
###############################################################################

control_plots <- unique(
  control_plot_daily$plot
)

removal_plots <- unique(
  removal_plot_daily$plot
)


control_compositions <- make_boot_compositions(
  length(control_plots)
)

removal_compositions <- make_boot_compositions(
  length(removal_plots)
)


cat(
  "Unique Control compositions:",
  nrow(control_compositions),
  "\n"
)

cat(
  "Unique Removal compositions:",
  nrow(removal_compositions),
  "\n"
)

cat(
  "Total Control x Removal combinations:",
  nrow(control_compositions) *
    nrow(removal_compositions),
  "\n\n"
)


###############################################################################
# 4. CHECK THAT THE BOOTSTRAP PROBABILITIES SUM TO 1
###############################################################################

cat(
  "Sum Control probabilities:",
  sum(control_compositions$weight),
  "\n"
)

cat(
  "Sum Removal probabilities:",
  sum(removal_compositions$weight),
  "\n\n"
)


###############################################################################
# 5. FUNCTION TO CALCULATE WEIGHTED TREATMENT MEAN
#    FOR ONE BOOTSTRAP COMPOSITION
###############################################################################

calculate_boot_treatment <- function(
    daily_data,
    plots,
    composition
) {

  counts <- as.numeric(
    composition[
      paste0(
        "n",
        seq_along(plots)
      )
    ]
  )


  ###########################################################################
  # Attach bootstrap multiplicity to each plot
  ###########################################################################

  plot_weights <- data.frame(
    plot = plots,
    boot_count = counts
  )


  dat <- daily_data %>%
    inner_join(
      plot_weights,
      by = "plot"
    ) %>%
    filter(
      boot_count > 0
    )


  ###########################################################################
  # Weighted mean across sampled plots
  #
  # A plot sampled 3 times receives weight 3.
  ###########################################################################

  dat %>%
    group_by(
      date,
      period
    ) %>%
    summarise(

      obs = weighted.mean(
        obs,
        w = boot_count,
        na.rm = TRUE
      ),

      mod = weighted.mean(
        mod,
        w = boot_count,
        na.rm = TRUE
      ),

      LE = weighted.mean(
        LE,
        w = boot_count,
        na.rm = TRUE
      ),

      RH = weighted.mean(
        RH,
        w = boot_count,
        na.rm = TRUE
      ),

      VPD = weighted.mean(
        VPD,
        w = boot_count,
        na.rm = TRUE
      ),

      Fd = weighted.mean(
        Fd,
        w = boot_count,
        na.rm = TRUE
      ),

      Fb = weighted.mean(
        Fb,
        w = boot_count,
        na.rm = TRUE
      ),

      Tsurf = weighted.mean(
        Tsurf,
        w = boot_count,
        na.rm = TRUE
      ),

      .groups = "drop"
    )
}


###############################################################################
# 6. CALCULATE ALL 35 CONTROL TREATMENT SERIES
###############################################################################

control_boot_series <- vector(
  "list",
  nrow(control_compositions)
)


for(i in seq_len(nrow(control_compositions))) {

  control_boot_series[[i]] <-
    calculate_boot_treatment(
      daily_data = control_plot_daily,
      plots = control_plots,
      composition = control_compositions[i, ]
    )

}


###############################################################################
# 7. CALCULATE ALL 35 REMOVAL TREATMENT SERIES
###############################################################################

removal_boot_series <- vector(
  "list",
  nrow(removal_compositions)
)


for(i in seq_len(nrow(removal_compositions))) {

  removal_boot_series[[i]] <-
    calculate_boot_treatment(
      daily_data = removal_plot_daily,
      plots = removal_plots,
      composition = removal_compositions[i, ]
    )

}


###############################################################################
# 8. RUN ALL 35 x 35 = 1225 COMBINATIONS
###############################################################################

exact_boot_results <- vector(
  "list",
  nrow(control_compositions) *
    nrow(removal_compositions)
)


counter <- 1


for(i in seq_len(nrow(control_compositions))) {

  for(j in seq_len(nrow(removal_compositions))) {


    #########################################################################
    # Join Control and Removal
    #########################################################################

    boot_comparison <- inner_join(

      control_boot_series[[i]],

      removal_boot_series[[j]],

      by = c(
        "date",
        "period"
      ),

      suffix = c(
        "_control",
        "_removal"
      )
    ) %>%

      mutate(

        season = ifelse(
          month(date) >= 5 &
            (
              month(date) < 12 |
                (
                  month(date) == 12 &
                    day(date) <= 15
                )
            ),
          "Wet season",
          "Dry season"
        ),


        #####################################################################
        # Removal - Control
        #####################################################################

        delta_obs =
          obs_removal -
          obs_control,

        delta_mod =
          mod_removal -
          mod_control,

        delta_LE =
          LE_removal -
          LE_control,

        delta_RH =
          RH_removal -
          RH_control,

        delta_VPD =
          VPD_removal -
          VPD_control,

        delta_Fd =
          Fd_removal -
          Fd_control,

        delta_Fb =
          Fb_removal -
          Fb_control,

        delta_SW =
          (
            Fd_removal +
              Fb_removal
          ) -
          (
            Fd_control +
              Fb_control
          ),

        delta_Tsurf =
          Tsurf_removal -
          Tsurf_control
      )


    #########################################################################
    # Calculate exactly the same median statistic as in the figures
    #########################################################################

    result <- boot_comparison %>%

      select(
        season,
        period,
        delta_obs,
        delta_mod,
        delta_LE,
        delta_RH,
        delta_VPD,
        delta_Fd,
        delta_Fb,
        delta_SW,
        delta_Tsurf
      ) %>%

      pivot_longer(
        cols = starts_with("delta_"),
        names_to = "variable",
        values_to = "effect"
      ) %>%

      group_by(
        season,
        period,
        variable
      ) %>%

      summarise(

        estimate = median(
          effect,
          na.rm = TRUE
        ),

        .groups = "drop"
      )


    #########################################################################
    # Probability of this Control x Removal combination
    #########################################################################

    result$weight <-
      control_compositions$weight[i] *
      removal_compositions$weight[j]


    result$control_combination <- i
    result$removal_combination <- j


    exact_boot_results[[counter]] <- result

    counter <- counter + 1
  }


  cat(
    "Finished Control composition",
    i,
    "of",
    nrow(control_compositions),
    "\n"
  )
}


###############################################################################
# 9. COMBINE RESULTS
###############################################################################

exact_boot_results <- bind_rows(
  exact_boot_results
)


###############################################################################
# 10. WEIGHTED QUANTILE FUNCTION
#
# Needed because the 1225 unique combinations do NOT all have the
# same bootstrap probability.
###############################################################################

weighted_quantile <- function(
    x,
    w,
    probs
) {

  ok <-
    is.finite(x) &
    is.finite(w) &
    w > 0

  x <- x[ok]
  w <- w[ok]


  ord <- order(x)

  x <- x[ord]
  w <- w[ord]


  w <- w / sum(w)

  cumulative_weight <- cumsum(w)


  sapply(
    probs,
    function(p) {

      x[
        which(
          cumulative_weight >= p
        )[1]
      ]

    }
  )
}


###############################################################################
# 11. EXACT 95% BOOTSTRAP CONFIDENCE INTERVALS
###############################################################################

exact_CI <- exact_boot_results %>%

  group_by(
    season,
    period,
    variable
  ) %>%

  summarise(

    CI_low = weighted_quantile(
      estimate,
      weight,
      0.025
    ),

    CI_high = weighted_quantile(
      estimate,
      weight,
      0.975
    ),

    .groups = "drop"
  )


###############################################################################
# 12. ORIGINAL MEDIAN EFFECTS FROM ACTUAL DATA
###############################################################################

original_effects <- comparison %>%

  mutate(

    delta_SW =
      delta_Fd +
      delta_Fb

  ) %>%

  select(
    season,
    period,
    delta_obs,
    delta_mod,
    delta_LE,
    delta_RH,
    delta_VPD,
    delta_Fd,
    delta_Fb,
    delta_SW,
    delta_Tsurf
  ) %>%

  pivot_longer(
    cols = starts_with("delta_"),
    names_to = "variable",
    values_to = "effect"
  ) %>%

  group_by(
    season,
    period,
    variable
  ) %>%

  summarise(

    median = median(
      effect,
      na.rm = TRUE
    ),

    .groups = "drop"
  )


###############################################################################
# 13. COMBINE MEDIAN + EXACT 95% BOOTSTRAP CI
###############################################################################

results_exact_CI <- original_effects %>%

  left_join(
    exact_CI,
    by = c(
      "season",
      "period",
      "variable"
    )
  ) %>%

  arrange(
    season,
    period,
    variable
  )


###############################################################################
# 14. PRINT FINAL RESULTS
###############################################################################

print(
  results_exact_CI,
  n = Inf
)


###############################################################################
# APPENDIX TABLE: ABSOLUTE VALUES + TREATMENT EFFECTS + 95% CI
###############################################################################

# Absolute seasonal medians for both treatments
absolute_summary <- comparison %>%
  mutate(
    SW_control = Fd_control + Fb_control,
    SW_removal = Fd_removal + Fb_removal
  ) %>%
  group_by(season, period) %>%
  summarise(
    obs_control     = median(obs_control, na.rm = TRUE),
    obs_removal     = median(obs_removal, na.rm = TRUE),

    mod_control     = median(mod_control, na.rm = TRUE),
    mod_removal     = median(mod_removal, na.rm = TRUE),

    Tsurf_control   = median(Tsurf_control, na.rm = TRUE),
    Tsurf_removal   = median(Tsurf_removal, na.rm = TRUE),

    LE_control      = median(LE_control, na.rm = TRUE),
    LE_removal      = median(LE_removal, na.rm = TRUE),

    RH_control      = median(RH_control, na.rm = TRUE),
    RH_removal      = median(RH_removal, na.rm = TRUE),

    VPD_control     = median(VPD_control, na.rm = TRUE),
    VPD_removal     = median(VPD_removal, na.rm = TRUE),

    Fd_control      = median(Fd_control, na.rm = TRUE),
    Fd_removal      = median(Fd_removal, na.rm = TRUE),

    Fb_control      = median(Fb_control, na.rm = TRUE),
    Fb_removal      = median(Fb_removal, na.rm = TRUE),

    SW_control      = median(SW_control, na.rm = TRUE),
    SW_removal      = median(SW_removal, na.rm = TRUE),

    .groups = "drop"
  )


# Convert to long format
absolute_long <- bind_rows(

  absolute_summary %>%
    transmute(
      season, period,
      variable = "delta_obs",
      control = obs_control,
      removal = obs_removal
    ),

  absolute_summary %>%
    transmute(
      season, period,
      variable = "delta_mod",
      control = mod_control,
      removal = mod_removal
    ),

  absolute_summary %>%
    transmute(
      season, period,
      variable = "delta_Tsurf",
      control = Tsurf_control,
      removal = Tsurf_removal
    ),

  absolute_summary %>%
    transmute(
      season, period,
      variable = "delta_LE",
      control = LE_control,
      removal = LE_removal
    ),

  absolute_summary %>%
    transmute(
      season, period,
      variable = "delta_RH",
      control = RH_control,
      removal = RH_removal
    ),

  absolute_summary %>%
    transmute(
      season, period,
      variable = "delta_VPD",
      control = VPD_control,
      removal = VPD_removal
    ),

  absolute_summary %>%
    transmute(
      season, period,
      variable = "delta_Fd",
      control = Fd_control,
      removal = Fd_removal
    ),

  absolute_summary %>%
    transmute(
      season, period,
      variable = "delta_Fb",
      control = Fb_control,
      removal = Fb_removal
    ),

  absolute_summary %>%
    transmute(
      season, period,
      variable = "delta_SW",
      control = SW_control,
      removal = SW_removal
    )
)


# Combine absolute values with treatment effects and bootstrap CIs
appendix_table <- results_exact_CI %>%

  left_join(
    absolute_long,
    by = c(
      "season",
      "period",
      "variable"
    )
  ) %>%

  mutate(

    variable_label = recode(
      variable,

      delta_obs =
        "Near-surface air temperature (observed)",

      delta_mod =
        "Near-surface air temperature (modelled)",

      delta_Tsurf =
        "Forest surface temperature",

      delta_LE =
        "Latent heat flux",

      delta_RH =
        "Relative humidity",

      delta_VPD =
        "Vapour pressure deficit",

      delta_Fd =
        "Diffuse shortwave radiation",

      delta_Fb =
        "Direct shortwave radiation",

      delta_SW =
        "Total shortwave radiation"
    ),

    unit = recode(
      variable,

      delta_obs = "°C",
      delta_mod = "°C",
      delta_Tsurf = "°C",

      delta_LE = "W m-2",
      delta_Fd = "W m-2",
      delta_Fb = "W m-2",
      delta_SW = "W m-2",

      delta_RH = "%",
      delta_VPD = "kPa"
    )

  ) %>%

  select(
    season,
    period,
    variable_label,
    unit,
    control,
    removal,
    median,
    CI_low,
    CI_high
  ) %>%
  rename(
    treatment_difference = median
  )


print(
  appendix_table,
  n = Inf
)


write.csv(
  appendix_table,
  file.path(
    output_path,
    "appendix_liana_effect_summary.csv"
  ),
  row.names = FALSE
)
