################################################################################
# ForEdgeClim annual simulation and shortwave-radiation validation
#
# Runs ForEdgeClim for one forest plot for every hourly timestep of 2023 using
# the calibrated shortwave-radiation parameters.
#
# For each timestep, the script:
#   - runs ForEdgeClim using the plot-specific TLS forest structure;
#   - extracts mean microclimate and energy-balance variables at 1 m height
#     (Tair, Tsurf, LE, RH, VPD, and direct/diffuse shortwave radiation);
#   - matches simulated air temperature to TOMST sensor locations;
#   - matches simulated PAR to PAR sensor locations when observations are
#     available;
#   - stores the resulting annual model timeseries.
#
# Shortwave-radiation performance is evaluated using daily median PAR per
# sensor, retaining sensor-days with at least 8 valid hourly observations.
# Temporal Pearson and Spearman correlations are calculated per sensor and
# spatial Pearson and Spearman correlations per day.
#
# Author: Emma Van de Walle – Q-ForestLab
################################################################################

# ---- Packages ----
suppressPackageStartupMessages({
  library(ForEdgeClimOnPanama)
  library(dplyr)
  library(Metrics)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(future)
  library(future.apply)
  library(RANN)
  library(progressr)
})

# ---- Parallel ----
handlers(global = TRUE)
handlers("txtprogressbar")
options(progressr.enable = TRUE)

cores <- 16
plan(multicore, workers = cores)
print(paste("Using", cores, "cores."))

start_analysis <- Sys.time()

# -----------------
# INPUT
# -----------------

plotnr = "C8" # C8 or R7
season_string = "C8_calibrated"

params = list(
  # uncalibrated values:
  # omega = 0.52, Kb_v = 1.25,

  # calibrated values:
  omega = 0.5266289, Kb_v = 0.5058042,

  # fixed values:
  beta0 = 0.325, betad = 0.325,
  omega_g_v = 0.13, Kd_v = 0.775,

  e_forest = 0.965, beta_lw = 0.325, omega_lw = 0.035,
  Kd_lw_v = 0.3, omega_g_lw_v = 0.055,

  h = 10,
  p_ground = 0.225,
  k_soil = 1.225,

  infl_macro_ref = 32.5,
  infl_soil_ref = 5,
  infl_forest_ref = 5,
  g_macro_ref = 25,
  g_soil_ref = 10,
  g_forest_ref = 12.5,

  gs_max = 0.004,
  stomata_closure = 0.75,

  k_u = 0.7,
  u_ref = 1.6
)

# ---- Time ----
all_datetimes <- seq(
  from = as.POSIXct("2023-01-01 00:00:00", tz = "UTC"),
  to   = as.POSIXct("2023-12-31 23:00:00", tz = "UTC"),
  by   = "hour"
)

## only daytime hours
#all_datetimes <- all_datetimes[
#  as.integer(format(all_datetimes, "%H")) >= 11 &
#  as.integer(format(all_datetimes, "%H")) <= 22
#]


# plot specifics
if(plotnr == "C8"){
  plot_map = "C8"
  height_canopy = 38
  scale_factor = 5.81/6.22
  voxel_TLS <- readRDS("Data/TLS_scaled_DTM_and_grid_C8_Panama.rds")
} else { # plotnr == "R7"
  plot_map = "R7"
  height_canopy = 37
  scale_factor = 6.19/6.22
  voxel_TLS <- readRDS("Data/TLS_scaled_DTM_and_grid_R7_Panama.rds")
}

# -----------------
# MAIN LOOP
# -----------------
process_datetime <- function(dt, i) {

  list2env(params, envir = .GlobalEnv)
  create_input_drivers()
  create_physical_constants()

  datetime <- as.POSIXct(dt, tz = "UTC")

  out <- tryCatch({

    # ---- Imports ----
    import_temp_observations(datetime)
    import_rad_observations(datetime)
    import_RH_u_observations(datetime)

    # ---- Model ----
    grid_local <- voxel_TLS$grid
    grid_local$density <- grid_local$density * scale_factor

    res <- run_foredgeclim(grid_local, datetime)

    # air temp
    micro_grid <- res$micro_grid
    air_temp   <- res$air_temperature

    if (any(!is.finite(air_temp)))
      stop("NaN/Inf in air_temp")

    # ---- Extract z = 1 ----
    micro_grid <- res$micro_grid
    sel <- micro_grid$z == 1

    surf_temp_mean <- mean(micro_grid$temperature[sel] - 273.15, na.rm = TRUE)
    LE_mean <- mean(res$latent_flux[sel], na.rm = TRUE)
    RH_mean <- mean(res$RH[sel], na.rm = TRUE)
    VPD_mean <- mean(res$VPD[sel], na.rm = TRUE)
    Fd_mean <- mean(res$sw_rad$F_d_down[sel], na.rm = TRUE)
    Fb_mean <- mean(res$sw_rad$F_b_down[sel], na.rm = TRUE)

    # ---- Observations ----
    key <- format(datetime, "%Y%m%d_%H")

    # PAR observations
    file_par <- paste0(
      "Data/PAR_obs_",
      key,
      "_",
      plotnr,
      ".csv"
    )

    if (!file.exists(file_par))
      return(NULL)

    PAR_obs <- read.csv(file_par)

    # TOMST observations
    TOMST <- read.csv(
      paste0(
        "Data/microtemp_air_obs_",
        key,
        "_",
        plotnr,
        ".csv"
      )
    )

    # ---- Model grids ----
    par_grid <- res$sw_rad %>%
      filter(Z == 1) %>%
      mutate(
        PAR_model = (F_d_down + F_b_down) * 2.02
      ) %>%
      select(
        X,
        Y,
        PAR_model
      )

    temp_air_grid <- micro_grid %>%
      mutate(
        temperature = air_temp - 273.15
      ) %>%
      filter(z == 1) %>%
      select(x, y, temperature)

    # ---- Mean model temperature at 1 m height ----
    Tair_mean_1m <- mean(temp_air_grid$temperature, na.rm = TRUE)

    # ---- Air temperature nearest neighbour ----
    nn_temp <- nn2(
      data = as.matrix(temp_air_grid[, c("x","y")]),
      query = as.matrix(TOMST[, c("x","y")]),
      k = 1
    )

    idx_temp <- nn_temp$nn.idx[,1]

    TOMST <- TOMST %>%
      mutate(
        temperature_model = temp_air_grid$temperature[idx_temp],
        Tair_mean_1m = Tair_mean_1m,
      )

    validation <- NULL

    # ============================================================
    # CASE 1: PAR observations available
    # ============================================================

    if (nrow(PAR_obs) > 0) {

      nn <- nn2(
        data = as.matrix(par_grid[, c("X","Y")]),
        query = as.matrix(PAR_obs[, c("x","y")]),
        k = 1
      )

      idx <- nn$nn.idx[,1]

      PAR_obs <- PAR_obs %>%
        left_join(
          TOMST %>%
            select(
              x,
              y,
              value_mean,
              temperature_model,
              Tair_mean_1m
            ),
          by = c("x","y")
        ) %>%
        mutate(
          datetime = datetime,

          PAR_model = par_grid$PAR_model[idx],
          resid = PAR_model - radiation_mean,

          LE_mean = LE_mean,
          RH_mean = RH_mean,
          VPD_mean = VPD_mean,
          Fd_down_mean = Fd_mean,
          Fb_down_mean = Fb_mean,
          Tsurf_mean = surf_temp_mean
        )

      valid <-
        is.finite(PAR_obs$radiation_mean) &
        is.finite(PAR_obs$PAR_model)

      if (any(valid)) {

        validation <- tibble(
          datetime = datetime,
          sensor_id = PAR_obs$location_ID_number[valid],
          obs = PAR_obs$radiation_mean[valid],
          mod = PAR_obs$PAR_model[valid],
          resid = PAR_obs$resid[valid],
          sensor_x = PAR_obs$x[valid],
          sensor_y = PAR_obs$y[valid]
        )

      }

      timeseries <- PAR_obs

    } else {

      # ============================================================
      # CASE 2: No PAR observations (night)
      # ============================================================

      timeseries <- TOMST %>%
        transmute(
          datetime,

          x,
          y,

          value_mean,
          temperature_model,
          Tair_mean_1m,

          radiation_mean = NA_real_,
          PAR_model      = NA_real_,
          resid          = NA_real_,

          LE_mean = LE_mean,
          RH_mean = RH_mean,
          VPD_mean = VPD_mean,
          Fd_down_mean = Fd_mean,
          Fb_down_mean = Fb_mean,
          Tsurf_mean = surf_temp_mean
        )

    }

    list(
      timeseries = timeseries,
      validation = validation
    )

  }, error = function(e) {
    warning(paste("Crash:", format(datetime), "-", e$message))
    NULL
  })

  out
}

# run parallel, but keep index
with_progress({
  p <- progressor(along = all_datetimes)

  parts <- future_lapply(
    seq_along(all_datetimes),
    FUN = function(i) {
      p()  # <-- update progress
      process_datetime(all_datetimes[i], i)
    },
    future.seed = TRUE
  )
})

# -----------------
# COMBINE OUTPUT
# -----------------

timeseries_all <- bind_rows(lapply(parts, `[[`, "timeseries"))
val_data       <- bind_rows(lapply(parts, `[[`, "validation"))

val_data_daily <- val_data %>%
  mutate(
    date = as.Date(datetime)
  ) %>%
  group_by(date, sensor_id) %>%
  summarise(
    obs = median(obs, na.rm = TRUE),
    mod = median(mod, na.rm = TRUE),
    n_hours = n(),
    sensor_x = first(sensor_x),
    sensor_y = first(sensor_y),
    .groups = "drop"
  ) %>%
  filter(n_hours >= 8) %>% # only keep days in which more than 8 hourly timesteps are taken
  mutate(
    resid = mod - obs
  )

if(nrow(val_data_daily) == 0){
  stop("No daily validation data remaining after n_hours filter.")
}

cat("\nDistribution of valid hourly observations per sensor-day:\n")
print(table(val_data_daily$n_hours))

cat(
  sprintf(
    "\nDaily validation dataset contains %d sensor-days\n",
    nrow(val_data_daily)
  )
)

# ---- Save ----
dir.create("Output/timeseries", recursive = TRUE, showWarnings = FALSE)
saveRDS(timeseries_all, paste0("Output/timeseries/timeseries_", season_string, ".rds"))

dir.create("Output/validation/data", recursive = TRUE, showWarnings = FALSE)
saveRDS(val_data_daily, paste0("Output/validation/data/val_data_", season_string, ".rds"))

# -----------------
# VALIDATION STATS
# -----------------

# -----------------
# TEMPORAL VALIDATION
# -----------------

temporal_metrics <- val_data_daily %>%
  group_by(sensor_id) %>%
  summarise(

    temporal_pearson =
      ifelse(
        n() > 2,
        cor(
          obs,
          mod,
          method = "pearson",
          use = "complete.obs"
        ),
        NA_real_
      ),

    temporal_spearman =
      ifelse(
        n() > 2,
        cor(
          obs,
          mod,
          method = "spearman",
          use = "complete.obs"
        ),
        NA_real_
      ),

    n_days = n(),

    .groups = "drop"
  )

saveRDS(
  temporal_metrics,
  paste0(
    "Output/validation/data/temporal_metrics_",
    season_string,
    ".rds"
  )
)

# -----------------
# SPATIAL VALIDATION
# -----------------

spatial_metrics <- val_data_daily %>%
  group_by(date) %>%
  summarise(

    n_sensors = n(),

    spatial_pearson =
      ifelse(
        n() > 2,
        cor(
          obs,
          mod,
          method = "pearson",
          use = "complete.obs"
        ),
        NA_real_
      ),

    spatial_spearman =
      ifelse(
        n() > 2,
        cor(
          obs,
          mod,
          method = "spearman",
          use = "complete.obs"
        ),
        NA_real_
      ),

    .groups = "drop"
  )

saveRDS(
  spatial_metrics,
  paste0(
    "Output/validation/data/spatial_metrics_",
    season_string,
    ".rds"
  )
)

# -----------------
# SUMMARY
# -----------------

cat("\nLight validation statistics:\n")

cat(sprintf(
  "\nTemporal Pearson r (median across sensors): %.2f\n",
  median(
    temporal_metrics$temporal_pearson,
    na.rm = TRUE
  )
))

cat(sprintf(
  "Temporal Spearman rho (median across sensors): %.2f\n",
  median(
    temporal_metrics$temporal_spearman,
    na.rm = TRUE
  )
))

cat(sprintf(
  "\nSpatial Pearson r (median across days): %.2f\n",
  median(
    spatial_metrics$spatial_pearson,
    na.rm = TRUE
  )
))

cat(sprintf(
  "Spatial Spearman rho (median across days): %.2f\n",
  median(
    spatial_metrics$spatial_spearman,
    na.rm = TRUE
  )
))

# ---- Cleanup ----
plan(sequential)

end_analysis <- Sys.time()
cat(sprintf('\nTotal running time = %.2f min\n',
            as.numeric(end_analysis - start_analysis, units = "mins")))
