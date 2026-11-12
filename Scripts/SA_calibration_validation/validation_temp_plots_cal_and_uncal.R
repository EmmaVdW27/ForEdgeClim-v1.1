############################################################################################
# ForEdgeClim air-temperature validation visualisation
#
# Compares air-temperature validation performance between calibrated and
# uncalibrated ForEdgeClim simulations for one forest plot using TOMST TMS-4
# observations.
#
# Validation statistics (ME, residual SD, RMSE, R2, and NSE) are calculated
# across all available hourly sensor observations for each parameterisation.
# The script visualises modelled versus observed air temperature and the
# distribution of residuals by hour of day.
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
base_plot <- "C13"   # C13 or R11

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

read_val <- function(plot_string, model_label){
  fp <- file.path("Output", "validation/data", paste0("val_data_", plot_string, ".rds"))
  if(!file.exists(fp)){
    stop(sprintf("File not found: %s", fp))
  }
  readRDS(fp) %>%
    mutate(model = model_label) %>%
    # ensure expected columns exist
    relocate(model, .before = 1)
}

val_list <- imap(plot_map, ~ read_val(.x, .y))  # .y = model label
val_data_all <- bind_rows(val_list) %>%
  mutate(model = factor(model, levels = c("uncalibrated", "calibrated")))


# For convenience keep the calibrated subset as before
val_data <- val_data_all %>% filter(model == "calibrated")

# -----------------
# STATISTICS (per model)
# -----------------

stats_df <- val_data_all %>%
  group_by(model) %>%
  summarise(
    ME   = mean(resid, na.rm = TRUE), # mean error (bias)
    RMSE = sqrt(mean(resid^2, na.rm = TRUE)), # root mean square error (spread)
    SD = sd(resid, na.rm = TRUE), # standard deviation of residuals (spread)
    R2   = cor(obs, mod, use = "complete.obs")^2, # R2 (model fit)
    NSE  = 1 - sum(resid^2, na.rm = TRUE) / # Nash–Sutcliffe model efficiency coefficient (overall performance)
      sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE),
    .groups = "drop"
  )

# Labels for annotations per model (for the scatter/facets)
lab_df <- stats_df %>%
  transmute(
    model,
    lab_scatter = sprintf("R\u00B2 = %.2f\nNSE = %.2f", R2, NSE),
    lab_resid   = sprintf("RMSE = %.2f \u00B0C\nME = %.2f \u00B0C\nSD = %.2f \u00B0C", RMSE, ME, SD)
  )

# -----------------
# MODEL vs OBSERVATION (all timesteps x positions) — facet per model
# -----------------

# Position of the 1:1 label per facet
one2one_pos <- val_data_all %>%
  group_by(model) %>%
  summarise(
    x_min = min(obs, na.rm = TRUE),
    y_min = min(obs, na.rm = TRUE),
    .groups = "drop"
  )

m_o_time_and_position <- ggplot(val_data_all, aes(x = obs, y = mod)) +
  geom_point(alpha = 0.3, colour = "blue") +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 1) +
  labs(x = "Observed temperature (°C)", y = "Modelled temperature (°C)",
       title = "(a) Modelled vs observed air temperature") +
  theme_bw() +
  theme(plot.title = element_text(size = 40),
        axis.title = element_text(size = 40),
        axis.text = element_text(size = 40),
        strip.text = element_text(size = 40)) +
  facet_wrap(~ model, ncol = 2) +
  # per-facet stats label
  geom_label(
    data = lab_df,
    aes(x = -Inf, y = Inf, label = lab_scatter),
    inherit.aes = FALSE,
    hjust = -0.1, vjust = 1.2,
    label.size = 0.6,
    label.r = unit(0.15, "lines"),
    fill = "white", colour = "black",
    size = 13
  ) +
  # 1:1 annotation per facet
  geom_text(
    data = one2one_pos,
    aes(x = x_min, y = y_min, label = "1:1"),
    inherit.aes = FALSE,
    colour = "red", hjust = -2, vjust = -1, fontface = "bold", size = 13
  )

#print(m_o_time_and_position)
ggsave(file.path(output_path, paste0("model_vs_obs_over_time_and_position_", base_plot, "_both.png")),
       plot = m_o_time_and_position, width = 16, height = 10, dpi = 300)



# -----------------
# RESIDUALS: grouped by hour-of-day (24h) per model
# -----------------

val_data_all_h <- val_data_all %>%
  mutate(
    hour = lubridate::hour(datetime),
    hour_f = factor(hour, levels = 0:23, labels = sprintf("%02d", 0:23))
  )

hour_residual <- ggplot(val_data_all_h, aes(x = hour_f, y = resid)) +
  geom_boxplot(
    fill = NA,
    colour = "blue",
    linewidth = 1,
    outlier.colour = "blue",
    outlier.alpha = 0.4
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red", linewidth = 1) +
  labs(
    y = "Residual (°C)",
    x = "Hour of day",
    title = "(b) Air temperature residuals grouped by hour of day"
  ) +
  theme_bw() +
  scale_x_discrete(
    breaks = sprintf("%02d", c(0, 6, 12, 18))
  ) +
  theme(
    plot.title  = element_text(size = 40),
    axis.title  = element_text(size = 40),
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 40),
    strip.text  = element_text(size = 40)
  ) +
  facet_wrap(~ model, ncol = 2) +
  geom_label(
    data = lab_df,
    aes(x = -Inf, y = Inf, label = lab_resid),
    inherit.aes = FALSE,
    hjust = -0.1, vjust = 3.3,
    label.size = 0.6,
    label.r = unit(0.15, "lines"),
    fill = "white",
    colour = "black",
    size = 13
  )


#print(hour_residual)
ggsave(
  file.path(output_path, paste0("residual_by_hour_", base_plot, "_both.png")),
  plot = hour_residual, width = 16, height = 10, dpi = 300
)




# -----------------
# CONSOLE SUMMARIES (per model)
# -----------------

cat("\nValidation statistics per model:\n")
stats_df %>%
  mutate(
    ME   = sprintf("%.2f °C", ME),
    RMSE = sprintf("%.2f °C", RMSE),
    SD   = sprintf("%.2f °C", SD),
    R2   = sprintf("%.2f",   R2),
    NSE  = sprintf("%.2f",   NSE)
  ) %>%
  rename(
    `Model` = model,
    `ME (bias)` = ME,
    `RMSE (spread)` = RMSE,
    `SD (spread)` = SD,
    `R2 (model fit)` = R2,
    `NSE (overall performance)` = NSE
  ) %>%
  print(row.names = FALSE)

