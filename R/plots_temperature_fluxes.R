#' Function to create several temperature plots
#'
#' @param micro_grid Grid with micro variables
#' @param T_air_vec Vector with air temperatures
#' @param output_path The output_path to store the plots
#' @param datetime yyyy-mm-dd hh:mm:ss in UTC
#' @return Several temperature plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @importFrom scales rescale
#' @import ggplot2 png grid
#' @importFrom ggtext element_markdown
#' @export
#'
plots_temp <- function(micro_grid, T_air_vec, output_path, datetime){

  # Define grid with surface temperatures
  temp_surf_grid = micro_grid
  # Convert surface temperature values from Kelvin to degrees Celsius
  temp_surf_grid$temperature = micro_grid$temperature - 273.15
  temp_surf_grid$T_soil = ifelse(micro_grid$z == 1, micro_grid$T_soil - 273.15, NA)
  temp_surf_grid$density = den

  # Define grid with air temperatures
  temp_air_grid = micro_grid
  # Convert air temperature values from Kelvin to degrees Celsius
  temp_air_grid$temperature = T_air_vec - 273.15

  # Select the available voxel slice closest to y = 30 m.
  selected_y <- 30
  available_y <- sort(unique(temp_air_grid$y[is.finite(temp_air_grid$y)]))
  if (length(available_y) == 0L) stop("No valid y-coordinates found")
  actual_y <- available_y[which.min(abs(available_y - selected_y))]

  # Forest surface temperatures in one x-z slice, without y-averaging.
  # Keep only occupied voxels, as in the original surface plot.
  surface_temperature_slice <- temp_surf_grid |>
    dplyr::filter(y == actual_y, !is.na(den), den != 0) |>
    dplyr::transmute(x, z, temperature)

  # Air temperatures in the same x-z slice.
  air_temperature_slice <- temp_air_grid |>
    dplyr::filter(y == actual_y) |>
    dplyr::transmute(x, z, temperature)

  # Calculate ground surface temperature for forest surface
  surface_temperature_ground <- temp_surf_grid |>
    filter(den != 0) |>
    filter(z == 1)

  # Calculate central transect temperature for forest surface
  surface_temperature_central_transect <- temp_surf_grid |>
    filter(den != 0) |>
    filter(y == actual_y)

  # Colors to be used in the plots
  colors = c("blue", "lightblue", "orange", "red")

  # 2D tile plot:
  temp_surf_plot = ggplot(surface_temperature_slice, aes(x = x, y = z, fill = temperature)) +
    geom_tile()+
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      #limits = c(26, 44),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Forest surface temperature"),
         x = "Easting (m)", y = "Height (m)", fill = "Temperature (°C)")+
    coord_fixed(ratio = 1, clip = "off") +
    theme_bw() #+
  # theme(
  #   plot.title = element_text(size = 20),
  #   plot.caption = element_text(hjust = 0, color = "black"),
  #   axis.title = element_text(size = 20),
  #   axis.text = element_text(size = 16),
  #   strip.text = element_text(size = 18),
  #   legend.title = element_text(size = 18),
  #   legend.text  = element_text(size = 16),
  #   legend.position = "right"
  #   )


  ggsave(paste0(output_path, '/temp_surface.png'), plot = temp_surf_plot, width = 9, height = 3, dpi = 500)


  # Calculate the average of the soil temperature over all y-slices for each combination of x and z
  average_soil_temperature <- temp_surf_grid |>
    filter(T_soil != 0) |>
    group_by(x, z) |>
    summarise(mean_temperature = mean(T_soil, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  temp_soil_plot = ggplot(average_soil_temperature, aes(x = x, y = z, fill = mean_temperature)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Soil temperature"),
         x = "Easting (m)", y = NULL, fill = "Temperature (°C)") +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      # plot.title = element_text(size = 20),
      # plot.caption = element_text(hjust = 0, color = "black"),
      # axis.title = element_text(size = 20),
      # axis.text.x = element_text(size = 16),
      axis.text.y = element_blank(),
      # strip.text = element_text(size = 18),
      # legend.title = element_text(size = 18),
      # legend.text  = element_text(size = 16)
    )

  ggsave(paste0(output_path, '/temp_soil.png'), plot = temp_soil_plot, width = 9, height = 3, dpi = 500)

  # Calculate the average of the air temperature over all y-slices for each combination of x and z
  average_air_temperature <- temp_air_grid |>
    group_by(x, z) |>
    summarise(mean_temperature = mean(temperature, na.rm = TRUE), .groups = 'keep')

  # Calculate the air temperature at ground level
  ground_air_temperature <- temp_air_grid |>
    filter(z ==1)

  # Calculate the air temperature at central transect
  central_transect_air_temperature <- temp_air_grid |>
    filter(y == actual_y)

  # 2D tile plot:
  temp_air_plot = ggplot(air_temperature_slice, aes(x = x, y = z, fill = temperature)) +
    geom_tile()+
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      #limits = c(26, 44),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Air temperature"),
         y = "Height (m)", x = "Easting (m)", fill = "Temperature (°C)")+
    coord_fixed(ratio = 1) +
    theme_bw() #+
  # theme(
  #   plot.title = element_text(size = 20),
  #   plot.caption = element_text(hjust = 0, color = "black"),
  #   axis.title = element_text(size = 20),
  #   axis.text = element_text(size = 16),
  #   strip.text = element_text(size = 18),
  #   legend.title = element_text(size = 18),
  #   legend.text  = element_text(size = 16)
  # )

  ggsave(paste0(output_path, '/temp_air.png'), plot = temp_air_plot, width = 9, height = 3, dpi = 500)

  # 2D tile plot for ground level surface temperature
  temp_surf_plot_2D_ground = ggplot(surface_temperature_ground, aes(x = x, y = y, fill = temperature)) +
    geom_tile()+
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      #limits = c(25.5, 38),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = "Forest surface temperature at ground level",
         x = 'Easting (m)', y = "Norhting (m)", fill = "Temperature (°C)")+
    coord_fixed(ratio = 1) +
    theme_bw() #+
  # theme(
  #   plot.title = element_text(size = 20),
  #   plot.caption = element_text(hjust = 0, color = "black"),
  #   axis.title = element_text(size = 20),
  #   axis.text = element_text(size = 16),
  #   strip.text = element_text(size = 18),
  #   legend.title = element_text(size = 18),
  #   legend.text  = element_text(size = 16)
  # )

  ggsave(paste0(output_path, '/temp_surface_ground.png'), plot = temp_surf_plot_2D_ground, width = 6, height = 6, dpi = 300)

  # 2D tile plot for ground level air temperature
  temp_air_plot_2D_ground = ggplot(ground_air_temperature, aes(x = x, y = y, fill = temperature)) +
    geom_tile()+
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      #limits = c(25.5, 38),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = "Air temperature at ground level",
         x = 'Easting (m)', y = "Norhting (m)", fill = "Temperature (°C)")+
    coord_fixed(ratio = 1) +
    theme_bw() #+
  # theme(
  #   plot.title = element_text(size = 20),
  #   plot.caption = element_text(hjust = 0, color = "black"),
  #   axis.title = element_text(size = 20),
  #   axis.text = element_text(size = 16),
  #   strip.text = element_text(size = 18),
  #   legend.title = element_text(size = 18),
  #   legend.text  = element_text(size = 16)
  # )

  ggsave(paste0(output_path, '/temp_air_ground.png'), plot = temp_air_plot_2D_ground, width = 6, height = 6, dpi = 300)

  # 3D plots #
  ############

  # 3D plot averaged air temperature

  temp_air_plot_3D = ggplot(air_temperature_slice, aes(x = x, y = z, fill = temperature)) +
    geom_tile()+
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      limits = c(26, 33),
      breaks = seq(26, 33, by = 2),
      oob = scales::squish,
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = NULL,
         x = NULL, y = NULL, fill = NULL)+
    scale_y_continuous(
      position = "right"
    ) +
    coord_fixed(ratio = 1, clip = "off") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 20),
      legend.title = element_text(size = 18),
      legend.text  = element_text(size = 16),
      legend.position = "right",
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    )


  ggsave(paste0(output_path, '/3D_temp_air.pdf'), plot = temp_air_plot_3D, width = 4, height = 3, dpi = 500)

  # 3D plot air temperature at ground surface

  temp_air_plot_3D_ground = ggplot(ground_air_temperature, aes(x = x, y = y, fill = temperature)) +
    geom_tile()+
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      limits = c(26, 33),
      breaks = seq(26, 33, by = 2),
      oob = scales::squish,
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = NULL,
         x = NULL, y = NULL, fill = NULL)+
    coord_fixed(ratio = 1, clip = "off") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.position = 'none',
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    )


  ggsave(paste0(output_path, '/3D_temp_air_ground.pdf'), plot = temp_air_plot_3D_ground, width = 3, height = 3, dpi = 500)

  # 3D plot averaged forest surface temerpature

  temp_surf_plot_3D = ggplot(surface_temperature_slice, aes(x = x, y = z, fill = temperature)) +
    geom_tile()+
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      limits = c(25, 39),
      breaks = seq(25, 39, by = 3),
      oob = scales::squish,
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = NULL,
         x = NULL, y = NULL, fill = NULL)+
    scale_y_continuous(
      position = "right"
    ) +
    coord_fixed(ratio = 1, clip = "off") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 20),
      legend.title = element_text(size = 18),
      legend.text  = element_text(size = 16),
      legend.position = "right",
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    )


  ggsave(paste0(output_path, '/3D_temp_surface.pdf'), plot = temp_surf_plot_3D, width = 4, height = 3, dpi = 500)

  # 3D plot forest surface temperature at ground surface

  temp_surf_plot_3D_ground = ggplot(surface_temperature_ground, aes(x = x, y = y, fill = temperature)) +
    geom_tile()+
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      limits = c(25, 39),
      breaks = seq(25, 39, by = 3),
      oob = scales::squish,
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = NULL,
         x = NULL, y = NULL, fill = NULL)+
    coord_fixed(ratio = 1, clip = "off") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.position = 'none',
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    )


  ggsave(paste0(output_path, '/3D_temp_surface_ground.pdf'), plot = temp_surf_plot_3D_ground, width = 3, height = 3, dpi = 500)
}

#' Function to create several flux plots
#'
#' @param micro_grid Grid with micro variables
#' @param net_radiaiton Vector with net radiation
#' @param sensible_flux Vector with sensible heat fluxes
#' @param latent_flux Vector with latent heat fluxes
#' @param G Vector with ground heat fluxes
#' @param output_path The output_path to store the plots
#' @param datetime yyyy-mm-dd hh:mm:ss in UTC
#' @return Several flux plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @importFrom scales rescale
#' @import ggplot2
#' @export
#'
plots_flux <- function(micro_grid, net_radiation, sensible_flux, latent_flux, G, output_path, datetime){

  # Define grid with fluxes
  fluxes_grid = micro_grid
  fluxes_grid$net_radiation = net_radiation
  fluxes_grid$sensible_flux = sensible_flux
  fluxes_grid$latent_flux = latent_flux
  fluxes_grid$ground_flux = G

  # Calculate the average of the net radiation over all y-slices for each combination of x and z
  average_radiation <- fluxes_grid |>
    filter(net_radiation != 0) |>
    group_by(x, z) |>
    summarise(mean_radiation = mean(net_radiation, na.rm = TRUE), .groups = 'keep')


  # Colors to be used in the plots
  colors = c("blue", "lightblue", "orange", "red")

  # 2D tile plot:
  rad_plot = ggplot(average_radiation, aes(x = x, y = z, fill = mean_radiation)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Net radiation"),
         subtitle = "+ value ⮕ flux entering surface",
         x = "Easting (m)", y = "Height (m)", fill = "Rₙ (W/m²)")+
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "black")
    )

  ggsave(paste0(output_path, '/flux_net_radiation.png'), plot = rad_plot, width = 9, height = 3, dpi = 300)

  # Calculate the average of the sensible heat flux over all y-slices for each combination of x and z
  average_H <- fluxes_grid |>
    filter(sensible_flux != 0) |>
    group_by(x, z) |>
    summarise(mean_H = mean(sensible_flux, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  H_plot = ggplot(average_H, aes(x = x, y = z, fill = mean_H)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Sensible heat surface flux"),
         subtitle = "+ value ⮕ flux lost to air",
         x = "Easting (m)", y = "Height (m)", fill = "H (W/m²)")+
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "black")
    )

  ggsave(paste0(output_path, '/flux_sensible.png'), plot = H_plot, width = 9, height = 3, dpi = 300)

  # Calculate the average of the latent heat flux over all y-slices for each combination of x and z
  average_LE <- fluxes_grid |>
    filter(latent_flux != 0) |>
    group_by(x, z) |>
    summarise(mean_LE = mean(latent_flux, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  LE_plot = ggplot(average_LE, aes(x = x, y = z, fill = mean_LE)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Latent heat surface flux"),
         subtitle = "+ value ⮕ flux lost to air",
         x = "Easting (m)", y = "Height (m)", fill = "LE (W/m²)")+
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "black")
    )

  ggsave(paste0(output_path, '/flux_latent.png'), plot = LE_plot, width = 9, height = 3, dpi = 300)

  # Calculate the average of the ground heat flux over all y-slices for each combination of x and z
  average_G <- fluxes_grid |>
    filter(ground_flux != 0) |>
    group_by(x, z) |>
    summarise(mean_G = mean(ground_flux, na.rm = TRUE), .groups = 'keep')

  # 2D tile plot:
  ground_plot = ggplot(average_G, aes(x = x, y = z, fill = mean_G)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colors,
      values = rescale(c(0.25, 0.5, 0.75, 1)),
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Ground heat flux"),
         subtitle = "+ value ⮕ flux entering soil",
         x = "Easting (m)", y = NULL, fill = "G (W/m²)")+
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0, color = "black"),
      axis.text.y = element_blank()
    )

  ggsave(paste0(output_path, '/flux_ground.png'), plot = ground_plot, width = 9, height = 3, dpi = 300)


}
