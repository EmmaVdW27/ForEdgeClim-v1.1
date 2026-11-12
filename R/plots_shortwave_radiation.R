#' Function to create several shortwave radiation plots
#'
#' @param sw_rad Dataframe with shortwave radiative values
#' @param output_path The output_path to store the plots
#' @return Several radiation plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @importFrom scales rescale
#' @import ggplot2
#' @export
#'
plots_sw <- function(sw_rad, output_path){

  # 2D RTM plots
  ##############

  # Calculate averages across Y-slices
  final_avg_results <- sw_rad |>
    group_by(X, Z) |>
    summarise(avg_F_d_down = mean(F_d_down, na.rm = TRUE),
              avg_F_d_up = mean(F_d_up, na.rm = TRUE),
              avg_F_b_down = mean(F_b_down, na.rm = TRUE), .groups = 'keep')

  # Calculate values at ground surface, i.e., Z = 1
  final_ground_surface_results <- sw_rad |>
    filter(Z == 1, na.rm = TRUE)

  # Calculate values at central transect, i.e., Y = 30
  final_central_transect_results <- sw_rad |>
    filter(Y == 30, na.rm = TRUE) |>
    mutate(F_down_total = F_d_down + F_b_down)

  # F_d_down plot:
  F_d_down = ggplot(final_avg_results, aes(x = X, y = Z, fill = avg_F_d_down)) +
    geom_tile()+
    scale_fill_viridis_c(option = "inferno", #limits = c(0, 700), breaks = c(0, 50, 100, 250, 600),
                        # trans = scales::pseudo_log_trans(sigma = 5), oob = scales::squish,
                         guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Shortwave diffuse downward radiation"),
         x = "Easting (m)", y = "Height (m)", fill = bquote("Flux"~(W~m^{-2})) )+
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

  ggsave(paste0(output_path, '/SW_F_d_down.png'), plot = F_d_down, width = 9, height = 3, dpi = 500)


  # F_d_up plot:
  F_d_up = ggplot(final_avg_results, aes(x = X, y = Z, fill = avg_F_d_up)) +
    geom_tile()+
    scale_fill_viridis_c(option = "inferno", #limits = c(0, 700), breaks = c(0, 50, 100, 250, 600),
                        # trans = scales::pseudo_log_trans(sigma = 5), oob = scales::squish,
                         guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Shortwave diffuse upward radiation"),
         x = "Easting (m)", y = "Height (m)", fill = bquote("Flux"~(W~m^{-2})) )+
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

  ggsave(paste0(output_path, '/SW_F_d_up.png'), plot = F_d_up, width = 9, height = 3, dpi = 500)

  # F_b_down plot:
  F_b_down = ggplot(final_avg_results, aes(x = X, y = Z, fill = avg_F_b_down)) +
    geom_tile()+
    scale_fill_viridis_c(option = "inferno", #limits = c(0, 700), breaks = c(0, 50, 100, 250, 600),
                         #trans = scales::pseudo_log_trans(sigma = 5), oob = scales::squish,
                         guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Shortwave direct-beam radiation"),
         fill  = bquote("Flux"~(W~m^{-2})), x = "Easting (m)", y = "Height (m)")+
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

  ggsave(paste0(output_path, '/SW_F_b_down.png'), plot = F_b_down, width = 9, height = 3, dpi = 500)


  # F_d_down + F_b_down: single vertical slice at y = 30 m

  F_b_d_down <- ggplot(
    final_central_transect_results,
    aes(x = X, y = Z, fill = F_down_total)
  ) +
    geom_tile() +
    scale_fill_viridis_c(
      option = "inferno",
      guide = guide_colorbar(
        barwidth = 1,
        barheight = 10,
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    labs(
      title = "Shortwave direct-beam & diffuse downward radiation",
      x = "Easting (m)",
      y = "Height (m)",
      fill = bquote(Flux~(W~m^{-2}))
    ) +
    coord_fixed(ratio = 1) +
    theme_bw()

  ggsave(
    paste0(output_path, "/SW_F_b_d_down.png"),
    plot = F_b_d_down,
    width = 9,
    height = 3,
    dpi = 500
  )

  # F_down ground surface plot:
  F_down_surface_2D = ggplot(final_ground_surface_results, aes(x = X, y = Y, fill = F_b_down + F_d_down)) +
    geom_tile()+
    scale_fill_viridis_c(
      option = "inferno",
      #limits = c(0, 650),
      guide = guide_colorbar(
        barwidth = 1,
        barheight = 10,
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    labs(title = "Shortwave downward radiation at ground level",
         x = "Easting (m)", y = "Norhting (m)", fill = bquote(Flux~(W~m^{-2})))+
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

  ggsave(paste0(output_path, '/SW_F_b_d_down_ground.png'), plot = F_down_surface_2D, width = 6, height = 6, dpi = 300)

  # 3D plots #
  ############

  # F_down plot for 3D visualisation
  # Single vertical slice at y = 30 m

  F_down <- ggplot(
    final_central_transect_results,
    aes(x = X, y = Z, fill = F_down_total)
  ) +
    geom_tile() +
    scale_fill_viridis_c(
      option = "inferno",
      limits = c(0, 650),
      trans = "sqrt",
      breaks = c(0, 50, 100, 200, 400, 600),
      guide = guide_colorbar(
        barwidth = 1,
        barheight = 10,
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = NULL,
      fill = NULL
    ) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    scale_y_continuous(position = "right") +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 20),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.position = "right",
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    )

  ggsave(
    paste0(output_path, "/3D_SW_down.pdf"),
    plot = F_down,
    width = 4,
    height = 3,
    dpi = 500,
    bg = "transparent"
  )

  # F_down ground surface plot:
  F_down_surface = ggplot(final_ground_surface_results, aes(x = X, y = Y, fill = F_b_down + F_d_down)) +
    geom_tile()+
    scale_fill_viridis_c(
      option = "inferno",
      limits = c(0, 650),
      trans = "sqrt",
      breaks = c(0, 50, 100, 200, 400, 600),
      guide = guide_colorbar(
        barwidth = 1,
        barheight = 10,
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    labs(title = NULL,
         x = NULL, y = NULL, fill = NULL)+
    coord_fixed(ratio = 1) +
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

  ggsave(paste0(output_path, '/3D_SW_down_ground.pdf'), plot = F_down_surface, width = 3, height = 3, dpi = 500, bg = "transparent")

}
