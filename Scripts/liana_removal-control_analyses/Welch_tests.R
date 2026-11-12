################################################################################
# Four primary daytime Welch tests with Holm correction for multiple testing
#
# H1: Full-year 2023 spatial thermal heterogeneity (P95-P05) at 0-20 m canopy height.
# H2: Combined very cool + very warm availability in end-season windows.
#
# Independent replication: four plots per treatment.
# Effect direction: Liana removal minus Control.
#
# Author: Emma Van de Walle - Q-ForestLab
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

output_path <- "Output/Welch_tests"
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Read microhabitat metrics
microhabitat_raw <- bind_rows(
  readRDS("Output/thermal_niches/data/thermal_niche_metrics_control_calibrated.rds") %>%
    mutate(treatment = "Control"),
  readRDS("Output/thermal_niches/data/thermal_niche_metrics_removal_calibrated.rds") %>%
    mutate(treatment = "Liana removal")
)
required_h2 <- c("treatment", "plot", "seasonal_period", "period",
                 "thermal_niche", "niche_availability")
if (!all(required_h2 %in% names(microhabitat_raw))) {
  stop("Microhabitat RDS missing columns: ",
       paste(setdiff(required_h2, names(microhabitat_raw)), collapse = ", "))
}

# H2: calculate plot-level mean of the hourly combined extreme fractions.
# Restrict to hours where both classes are present, to avoid silently treating
# missing class observations as zero. A genuine zero should be stored as zero.
extreme_plot <- microhabitat_raw %>%
  group_by(treatment, plot, seasonal_period, period, thermal_niche) %>%
  summarise(niche_availability = mean(niche_availability, na.rm = TRUE) * 100, .groups = "drop") %>%
  filter(period == "Day", thermal_niche %in% c("Very cool", "Very warm")) %>%
  group_by(treatment, plot, seasonal_period) %>%
  summarise(
    response = if (n_distinct(thermal_niche) == 2L &&
                   all(is.finite(niche_availability)))
      sum(niche_availability) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    hypothesis = "H2: Combined extreme availability",
    seasonal_period = recode(as.character(seasonal_period),
                             "End dry season" = "Dry season",
                             "End wet season" = "Wet season")
  )

# H1: full-year 2023 heterogeneity
# Height: 0 <= z < 20 m.
heterogeneity_control <- readRDS(
  "Output/heterogeneity/heterogeneity_control_calibrated.rds"
) %>% mutate(treatment = "Control")
heterogeneity_removal <- readRDS(
  "Output/heterogeneity/heterogeneity_removal_calibrated.rds"
) %>% mutate(treatment = "Liana removal")

heterogeneity_raw <- bind_rows(heterogeneity_control, heterogeneity_removal)
required_h1 <- c("treatment", "plot", "date", "doy", "hour", "z", "Tair_range90")
if (!all(required_h1 %in% names(heterogeneity_raw))) {
  stop("Heterogeneity RDS missing columns: ",
       paste(setdiff(required_h1, names(heterogeneity_raw)), collapse = ", "))
}

heterogeneity_plot <- heterogeneity_raw %>%
  mutate(
    date = as.Date(date),
    seasonal_period = if_else(
      doy >= 350 | doy <= 120, "Dry season", "Wet season"
    )
  ) %>%
  filter(
    !is.na(date),
    format(date, "%Y") == "2023",
    hour >= 11, hour <= 22,
    z >= 0, z < 20,
    is.finite(Tair_range90)
  ) %>%
  # Average hourly P95-P05 per 1-m layer within each day, then average
  # equally across layers and days, retaining plot as the replication unit.
  group_by(treatment, plot, seasonal_period, date, z) %>%
  summarise(layer_day_range = mean(Tair_range90), .groups = "drop") %>%
  group_by(treatment, plot, seasonal_period, date) %>%
  summarise(day_range = mean(layer_day_range),
            n_layers = n_distinct(z), .groups = "drop") %>%
  group_by(treatment, plot, seasonal_period) %>%
  summarise(
    response = mean(day_range),
    n_days = n_distinct(date),
    min_layers_per_day = min(n_layers),
    .groups = "drop"
  ) %>%
  mutate(hypothesis = "H1: Thermal heterogeneity 0-20 m")

# Diagnostic output: inspect temporal coverage and layer completeness.
print(heterogeneity_plot, n = Inf)
write.csv(heterogeneity_plot,
          file.path(output_path, "H1_heterogeneity_plot_values.csv"),
          row.names = FALSE)
if (any(heterogeneity_plot$min_layers_per_day < 20)) {
  warning("Some H1 plot-days have fewer than 20 height layers; inspect coverage.")
}

all_plot_data <- bind_rows(
  heterogeneity_plot %>% select(treatment, plot, seasonal_period, response, hypothesis),
  extreme_plot %>% select(treatment, plot, seasonal_period, response, hypothesis)
) %>%
  mutate(treatment = as.character(treatment),
         seasonal_period = as.character(seasonal_period),
         plot = as.character(plot))

if (any(!is.finite(all_plot_data$response)))
  stop("Non-finite plot-level response detected. Check the input summaries.")
if (anyDuplicated(all_plot_data[c("hypothesis", "seasonal_period", "treatment", "plot")]))
  stop("Duplicate plot-level entries detected.")

expected_seasons <- c("Dry season", "Wet season")
expected_hypotheses <- c("H1: Thermal heterogeneity 0-20 m",
                         "H2: Combined extreme availability")
if (!setequal(unique(all_plot_data$seasonal_period), expected_seasons) ||
    !setequal(unique(all_plot_data$hypothesis), expected_hypotheses) ||
    !setequal(unique(all_plot_data$treatment), c("Control", "Liana removal")))
  stop("Unexpected or missing hypothesis, season or treatment labels.")

# One Welch test per hypothesis x season. Effect: removal minus control.
welch_one <- function(df) {
  x <- df$response[df$treatment == "Liana removal"]
  y <- df$response[df$treatment == "Control"]
  if (length(x) != 4L || length(y) != 4L)
    stop("Expected exactly four independent plots per treatment and test.")
  fit <- t.test(x, y, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  tibble(
    n_control = length(y), n_removal = length(x),
    mean_control = mean(y), sd_control = sd(y),
    mean_removal = mean(x), sd_removal = sd(x),
    mean_difference = mean(x) - mean(y),
    ci_low = unname(fit$conf.int[1]),
    ci_high = unname(fit$conf.int[2]),
    p_unadjusted = fit$p.value
  )
}

stats_primary <- all_plot_data %>%
  group_by(hypothesis, seasonal_period) %>%
  group_modify(~ welch_one(.x)) %>%
  ungroup() %>%
  mutate(
    p_holm = p.adjust(p_unadjusted, method = "holm"),
    significant_holm_05 = p_holm < 0.05,
    units = if_else(grepl("^H1", hypothesis), "deg C", "percentage points")
  ) %>%
  arrange(hypothesis, factor(seasonal_period, levels = expected_seasons))

if (nrow(stats_primary) != 4L) stop("Expected exactly four tests.")

cat("\nFOUR PRIMARY DAYTIME WELCH TESTS (HOLM-CORRECTED)\n")
cat("Effect = Liana removal minus Control; 95% CIs are unadjusted Welch CIs.\n\n")
print(stats_primary, n = Inf, width = Inf)

write.csv(all_plot_data,
          file.path(output_path, "primary_welch_plot_values.csv"),
          row.names = FALSE)
write.csv(stats_primary,
          file.path(output_path, "primary_welch_four_tests_holm.csv"),
          row.names = FALSE)

cat("\nStatistical results saved to: ", output_path, "\n")
