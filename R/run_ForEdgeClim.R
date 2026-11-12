###############################################################################
#
#                               ForEdgeClim
#
# A microclimate model that solves the energy balance equation in fragmented
# forest areas: Rn = LE + H + G, where,
# Rn = Net radiation (shortwave & longwave)
# LE = Latent heat flux (~ evapotranspiration)
# H = sensible heat flux
# G = Ground heat flux
###############################################################################


#' Execute the microclimate model ForEdgeClim
#'
#' @param structure_grid Dataframe with structural density values
#' @param datetime Datetime object representing the current simulation time
#' @return Dataframe with simulated microclimate temperatures and fluxes
#' @importFrom dplyr group_by ungroup mutate left_join
#' @importFrom RANN nn2
#' @export
run_foredgeclim <- function(structure_grid, datetime, trace_idx = NULL) {

  #############################################
  # 2D SHORTWAVE RADIATIVE TRANSFER MODELLING #
  #############################################
  voxel_grid = structure_grid

  final_results <<- shortwave_two_stream_RTM(datetime, lat, lon, voxel_grid,
                                             F_sky_diff_init, F_sky_dir_init,
                                             omega_g_v,
                                             Kd_v, Kb_v,
                                             omega, betad, beta0) # A lot of these parameters are global parameters and actually do not need to be stated as arguments.
  # saveRDS(final_results, 'Output/SW_DATA.rds')

  den <<- final_results$density
  net_SW <- final_results$net_sw
  SW_down <- final_results$F_d_down + final_results$F_b_down
  SW_up <- final_results$F_d_up
  x_coords <- voxel_grid$X
  y_coords <- voxel_grid$Y
  z_coords <<- voxel_grid$Z

  # dimensions grid
  x_dim <<- max(x_coords)  # Length in x-direction
  y_dim <<- max(y_coords)  # Length in y-direction
  z_dim <<- max(z_coords)  # Length in z-direction


  ###########################################
  # WIND MODULE
  ###########################################

  # Convert wind speed from km/h to m/s
  u_top_ms <- u_top / 3.6

  # Create a grouping ID for each (x, y) column
  # This replaces the need for looping over unique x-y combinations
  xy_id <- interaction(x_coords, y_coords, drop = TRUE)

  # Sort all voxels:
  # - first by column (xy_id)
  # - then by height (descending z from canopy top downwards)
  ord <- order(xy_id, -z_coords)

  # Apply sorting
  xy_id_sorted <- xy_id[ord]
  den_sorted <- den[ord]

  # Replace NA density values with 0
  den_sorted[is.na(den_sorted)] <- 0

  # Compute cumulative LAI per column (from top to bottom)
  # This replaces the inner loop that builds LAI_cum
  LAI_cum <- ave(den_sorted, xy_id_sorted, FUN = cumsum)

  # Determine where canopy starts in each column:
  # TRUE from the first voxel with structure downward
  has_canopy <- ave(den_sorted > 0, xy_id_sorted, FUN = function(x) cumsum(x) > 0)

  # Compute wind speed:
  # - above canopy: u_top_ms
  # - below canopy: attenuated exponentially based on cumulative LAI
  # - enforce minimum wind speed of 0.2 m/s
  u_sorted <- ifelse(
    !has_canopy,
    u_top_ms,
    pmax(u_top_ms * exp(-k_u * LAI_cum), 0.2)
  )

  # Reorder results back to original voxel order
  u <- numeric(length(den))
  u[ord] <- u_sorted

  # Conductances:
  g_wind <- pmax(0.1, sqrt(u / u_ref))

  g_forest <<- g_forest_ref * g_wind
  g_macro  <<- g_macro_ref  * g_wind

  # no wind dependence for soil conductance
  g_soil   <<- g_soil_ref

  # Distances of influence stay unadjusted:
  infl_soil   <<- infl_soil_ref
  infl_forest <<- infl_forest_ref
  infl_macro  <<- infl_macro_ref


  #############################
  # INITIALIsING TEMPERATURES #
  #############################

  # First longwave RTM for initializing temperatures
  micro_grid <- data.frame(
    x = x_coords,
    y = y_coords,
    z = z_coords,
    temperature = macro_temp,
    T_soil = macro_temp
  )

  # Inverse optical depth and auxiliary variables for longwave RTM
  mu_lw <- ifelse(den == 0, 1, 1 / (Kd_lw_v * den))
  n = length(den)
  gamma_plus_lw <<- rep(0.5 * (1 + sqrt((1 - omega_lw) / (1 - (1 - 2 * beta_lw) * omega_lw))) , n)
  gamma_plus_lw[n+1] <<- 1
  gamma_minus_lw <<- rep(0.5 * (1 - sqrt((1 - omega_lw) / (1 - (1 - 2 * beta_lw) * omega_lw))) , n)
  gamma_minus_lw[n+1] <<- 0
  lambda_lw <<- sqrt((1 - (1 - 2 * beta_lw) * omega_lw) * (1 - omega_lw) / (mu_lw)^2)
  lambda_lw[n+1] <<- 0

  final_results_lw = longwave_two_stream_RTM(voxel_grid, micro_grid, lw_two_stream,
                                             F_sky_lw, omega_g_lw_v, macro_temp,
                                             Kd_lw_v, omega_lw, beta_lw)
  net_LW <- final_results_lw$net_lw
  LW_down <- final_results_lw$F_d_down
  LW_up <- final_results_lw$F_d_up
  net_rad_ground <- ifelse(z_coords == 1, SW_down - SW_up + LW_down - LW_up, NA)

  # Net radiation forest surface
  net_radiation_init <- net_SW + net_LW

  # Define ground flux as a % of net radiation at ground
  G = calculate_G(net_rad_ground)

  # Define T_soil at ground surface from ground flux with ground conductance k_soil
  micro_grid$T_soil = G * stable_soil_depth / k_soil + T_soil_deep

  # Define micro_grid$temperature based on macro temp and initial net radiation + random noise
  micro_grid$temperature =  macro_temp + 0.01*net_radiation_init + rnorm(nrow(micro_grid), mean = 0, sd = 0.1)

  # First approximation of air temperature, based on forest surface temperature and random noise
  T_air_vec = micro_grid$temperature + rnorm(nrow(micro_grid), mean = 0, sd = 0.1)

  # Define T_ground for every voxel in a vertical xy-column as T_soil
  micro_grid <- micro_grid |>
    group_by(x, y) |>
    mutate(T_ground = T_soil[z == 1][1])  |>
    ungroup()

  # Determine weights for each component in the linearization via exponential weighting
  #x_threshold = length_transect # x-value from where macro influence is being felt, ie, forest edge value
  z_threshold = height_canopy  # z-value from where macro influence is being felt, ie, canopy top value
  #dis_macro_x = ifelse(x_coords > x_threshold, 0, x_threshold - x_coords + 1)
  dis_macro_z = ifelse(z_coords > z_threshold, 0, z_threshold - z_coords + 1)
  dis_soil = z_coords
  dis_forest = 0.5 # size of a voxel edge = 1m, so distance from structure is on average 0.5m
  alpha_macro = log(0.5) / infl_macro
  alpha_soil = log(0.5)/infl_soil
  alpha_forest = log(0.5)/infl_forest
  #w_macro_x = exp(alpha_macro*dis_macro_x)
  w_macro_z = exp(alpha_macro*dis_macro_z)
  w_soil = exp(alpha_soil*dis_soil)
  w_forest = exp(alpha_forest*dis_forest)

  # Air temperature update based on a linearization as is done in microclimc (Maclean and Klinges, 2021) and its successor microclimf.
  # The updated air temperature is weighted by the convection coefficient and the distance to each boundary.
  # This linearization simplifies the complex nonlinear interactions between leaves and air, making the calculations more efficient.
  # This linearization is mostly valid when there are small temperature differences between T air en T leaf, when there is
  # sufficient air streaming, when the net radiation or the humidity is not extreme, when the forest structure is quite homogeneous.
  # The linearization calculates the air temperature as influenced by the macro and soil boundary temperature and the surface temperature.
  T_air_vec = ( w_macro_z *g_macro * macro_temp  + w_soil*g_soil * micro_grid$T_ground + w_forest*g_forest * micro_grid$temperature ) / ( w_macro_z*g_macro + w_soil*g_soil + w_forest*g_forest)

  #########################################
  # ITERATION TO CLOSE THE ENERGY BALANCE #
  #########################################

  # Start iteration
  max_iter = 1000
  W = 1  # Weighting step, start value
  W_min = 0.01  # Weighting step, min value
  error_prev = Inf

  # Optional tracing of Newton convergence for selected voxels
  if (!is.null(trace_idx)) {
    trace_list <- vector("list", max_iter)  # one entry per iteration
  }


  for (iter in 1:max_iter) {

    #################
    # NET RADIATION #
    #################

    final_results_lw = longwave_two_stream_RTM(voxel_grid, micro_grid, lw_two_stream,
                                               F_sky_lw, omega_g_lw_v, macro_temp,
                                               Kd_lw_v, omega_lw, beta_lw) # A lot of these parameters are global parameters and actually do not need to be stated as arguments.


    net_LW = final_results_lw$net_lw
    LW_down = final_results_lw$F_d_down
    LW_up = final_results_lw$F_d_up
    net_rad_ground = ifelse(z_coords == 1, SW_down - SW_up + LW_down - LW_up, NA)

    # Net radiation for structure
    net_radiation = net_SW + net_LW


    ###############
    # Ground flux #
    ###############

    # G is positive if E is entering the soil
    G = calculate_G(net_rad_ground)
    micro_grid$T_soil = G*stable_soil_depth/k_soil + T_soil_deep
    micro_grid <- micro_grid |>
      group_by(x, y) |>
      mutate(T_ground = T_soil[z == 1][1])  |>
      ungroup()

    G_voxel <- rep(0, length(net_radiation))
    G_voxel[z_coords == 1] <- G[z_coords == 1]

    #############################
    # AIR TO AIR HEAT DIFFUSION #
    #############################

    dt <- 1
    T_air_vec = T_air_vec - calculate_D(T_air_vec, micro_grid$T_ground) * dt / (Cp * V * rho)


    ######################
    # SENSIBLE HEAT FLUX #
    ######################

    sensible_flux = calculate_H(micro_grid$temperature,T_air_vec)


    ####################
    # LATENT HEAT FLUX #
    ####################

    #latent_flux <- calculate_LE(micro_grid$temperature-273.15,net_radiation, G_voxel) # temp in calculate_LE must be in °C
    latent_flux_res <- calculate_LE_PM_C(micro_grid$temperature-273.15, T_air_vec-273.15, SW_down, net_radiation, G_voxel)
    latent_flux <- latent_flux_res$LE


    ##################
    # ENERGY BALANCE #
    ##################

    energy_balance_surf <- net_radiation - sensible_flux - latent_flux
    error_current = max(abs(energy_balance_surf))


    # Check oscillations, if so, reduce W
    if (iter > 1 && error_current > error_prev) {
      W = max(W * 0.8, W_min)  # Reduce W by 20%, but never smaller than W_min
    }
    error_prev = error_current  # Update previous error

    # Check convergence
    if (max(abs(energy_balance_surf)) < energy_balance_tolerance){
      #print(paste0('Convergence is reached after ', iter, ' iterations.'))
      break
    }

    # Use Newton's method to update surface temperature as is done in the SCOPE 2.0 model (Yang et al. 2021)
    #micro_grid$temperature <- micro_grid$temperature - ifelse(energy_balance_surf == 0, 0, W*energy_balance_surf/de_dT(micro_grid$temperature,net_radiation) )
    micro_grid$temperature <- micro_grid$temperature - ifelse(
      energy_balance_surf == 0,
      0,
      W*energy_balance_surf/de_dT_PM(micro_grid$temperature,
                                     SW_down,
                                     net_radiation,
                                     G_voxel)
    )

    # Determine planes-averaged surface temperature (to be used in T_air_Vec update for zero-density voxels).
    temp <- micro_grid$temperature
    mask<- as.numeric(den > 0)

    # --- Z direction ---
    sum_z   <- rowsum(temp * mask, z_coords)
    count_z <- rowsum(mask, z_coords)
    mean_z  <- sum_z / count_z
    mean_temp_z <- mean_z[z_coords]

    # --- Y direction ---
    sum_y   <- rowsum(temp * mask, y_coords)
    count_y <- rowsum(mask, y_coords)
    mean_y  <- sum_y / count_y
    mean_temp_y <- mean_y[y_coords]

    # --- X direction ---
    sum_x   <- rowsum(temp * mask, x_coords)
    count_x <- rowsum(mask, x_coords)
    mean_x  <- sum_x / count_x
    mean_temp_x <- mean_x[x_coords]

    # --- Combine ---
    mean_temp <- (mean_temp_x + mean_temp_y + mean_temp_z) / 3

    # Air temperature update based on a linearization as is done in microclimc (Maclean & Klinges, 2021).
    T_air_vec <- ( w_macro_z*g_macro * macro_temp + w_soil*g_soil * micro_grid$T_ground + ifelse(den == 0, w_forest*g_forest * mean_temp, w_forest*g_forest *micro_grid$temperature) )/ ( w_macro_z*g_macro + w_soil*g_soil + w_forest*g_forest)

  }

  ###########
  # RESULTS #
  ###########

  res <- list(
    sw_rad = final_results, # SW RTM output
    lw_rad = final_results_lw, # LW RTM output
    micro_grid = micro_grid,
    air_temperature = T_air_vec,
    net_radiation = net_radiation,
    sensible_flux = sensible_flux,
    latent_flux = latent_flux,
    ground_flux = G,
    RH = latent_flux_res$RH,
    VPD = latent_flux_res$VPD
  )

  return(res)

}



