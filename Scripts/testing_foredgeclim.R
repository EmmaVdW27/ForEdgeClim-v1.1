###############################################################################
#
# Example script demonstrating how to run the ForEdgeClim microclimate model.
#
# The script defines the required model inputs, imports meteorological forcing
# data, loads the TLS-derived 3D forest structure, runs ForEdgeClim for the
# specified time period, and saves and visualises the model output.
#
# Author: Emma Van de Walle - Q-ForestLab
#
###############################################################################


library(ForEdgeClimOnPanama)
start_script = Sys.time()

#########
# INPUT #
#########

plotnr <- "R11" # C13, R11, R7, C8, R5, R14, C12 or C16

# First hour to simulate:
start_time <- as.POSIXct("2023-04-26 19:00:00", tz = "UTC")
# Last hour to simulate:
end_time <- start_time
# ForEdgeClim is run for each hour within this period:
datetime_series <- seq(start_time, end_time, by = "hour")

#########

if(plotnr == "C13"){
  plot_map = "C13"
  height_canopy = 33
  # original TLS las file:
  TLS_input_file <- 'Data/Data_Kasper/plot13_trimmed_5cm.las'
  ALS_input_file <- 'Data/Data_Kasper/plot13_ALS.las'
  # If the structural voxel grid has already been generated using
  # generate_DTM_grid_TLS(), the processed TLS file can be loaded directly:
  TLS_filtered_file <- 'Data/TLS_scaled_DTM_and_grid_C13_Panama.rds'
} else if(plotnr == "R11"){
  plot_map = "R11"
  height_canopy = 35
  TLS_input_file <- 'Data/Data_Kasper/plot11_trimmed_5cm.las'
  ALS_input_file <- 'Data/Data_Kasper/plot11_ALS.las'
  TLS_filtered_file <- 'Data/TLS_scaled_DTM_and_grid_R11_Panama.rds'
} else if(plotnr == "C8"){
  plot_map = "C8"
  height_canopy = 38
  TLS_input_file <- 'Data/Data_Kasper/plot08_trimmed_5cm.las'
  ALS_input_file <- 'Data/Data_Kasper/plot08_ALS.las'
  TLS_filtered_file <- 'Data/TLS_scaled_DTM_and_grid_C8_Panama.rds'
} else if(plotnr == "R7"){
  plot_map = "R7"
  height_canopy = 37
  TLS_input_file <- 'Data/Data_Kasper/plot07_trimmed_5cm.las'
  ALS_input_file <- 'Data/Data_Kasper/plot07_ALS.las'
  TLS_filtered_file <- 'Data/TLS_scaled_DTM_and_grid_R7_Panama.rds'
} else if(plotnr == "R5"){
  plot_map = "R5"
  height_canopy = 36
  TLS_input_file <- 'Data/Data_Kasper/plot05_trimmed_5cm.las'
  TLS_filtered_file <- 'Data/TLS_scaled_DTM_and_grid_R5_Panama.rds'
} else if(plotnr == "R14"){
  plot_map = "R14"
  height_canopy = 42
  TLS_input_file <- 'Data/Data_Kasper/plot14_trimmed_5cm.las'
  TLS_filtered_file <- 'Data/TLS_scaled_DTM_and_grid_R14_Panama.rds'
} else if(plotnr == "C12"){
  plot_map = "C12"
  height_canopy = 29
  TLS_input_file <- 'Data/Data_Kasper/plot12_trimmed_5cm.las'
  TLS_filtered_file <- 'Data/TLS_scaled_DTM_and_grid_C12_Panama.rds'
} else {
  plot_map = "C16"
  height_canopy = 30
  TLS_input_file <- 'Data/Data_Kasper/plot16_trimmed_5cm.las'
  TLS_filtered_file <- 'Data/TLS_scaled_DTM_and_grid_C16_Panama.rds'
}

for (current_datetime in datetime_series) {

  # Make sure the current datetime is a POSIXct object:
  current_datetime <- as.POSIXct(current_datetime, origin = "1970-01-01", tz = "UTC")
  cat(paste0('Running model for plot ', plotnr, ' at ', current_datetime, ' UTC...\n'))
  cat("\n")

  #################################################
  # INPUT DRIVERS, CONSTANTS AND MODEL PARAMETERS #
  #################################################

  create_input_drivers()
  create_physical_constants()
  create_model_parameters()

  # Create output path:
  output_path <- file.path("Output", format(current_datetime, "%Y-%m-%d_%H"), plot_map)
  # Make sure the output path exists:
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  }

  ########################################################
  # EXTRACT THE OBSERVATIONS FROM THE INPUT DRIVER FILES #
  ########################################################

  # Import meteorological observations used as model drivers and for validation:
  import_temp_observations(current_datetime)
  import_rad_observations(current_datetime)
  import_RH_u_observations(current_datetime)

  #######################
  # 3D FOREST STRUCTURE #
  #######################

  # Generate the structural voxel grid from the TLS LAS file if needed:
  # voxel_TLS = generate_DTM_grid_TLS(TLS_input_file, ALS_input_file, voxel_length)
  # saveRDS(voxel_TLS, TLS_filtered_file)
  # Load the preprocessed TLS-derived voxel grid:
  voxel_TLS = readRDS(TLS_filtered_file)

  #############
  # RUN MODEL #
  #############

  res = run_foredgeclim(voxel_TLS$grid, current_datetime)

  saveRDS(res, paste0(output_path, '/model_results.rds'))

  ############
  # PLOTTING #
  ############

  # # Digital terrain model and structural grid plots:
  # plots_dtm_struct(dtm = voxel_TLS$dtm, grid = voxel_TLS$grid, output_path)
  #
  # # Shortwave radiation plots:
  # plots_sw(sw_rad = res$sw_rad, output_path)
  #
  # # Longwave radiation plots:
  # plots_lw(lw_rad = res$lw_rad, output_path)
  #
  # # Flux plots:
  # plots_flux(res$micro_grid, res$net_radiation, res$sensible_flux, res$latent_flux, res$ground_flux, output_path, current_datetime)
  #
  # # Temperature plots:
  # plots_temp(res$micro_grid, res$air_temperature, output_path, current_datetime)


}
end_script = Sys.time()
cat(paste0('Total run time = ', round(as.numeric(end_script - start_script, units = "secs"), 2), ' s'))
