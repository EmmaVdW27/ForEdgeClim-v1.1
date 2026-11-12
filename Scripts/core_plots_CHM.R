###############################################################################
# Canopy height and sensor location maps
#
# Creates canopy-height maps for the four intensively monitored plots in
# Panama (C8, C13, R7, and R11) from 1 m TLS-derived structural grids.
# Air temperature (TMS) sensor locations are shown for all plots, while PAR
# sensor locations are additionally shown for C8 and R7. The four maps are
# combined into a single figure with shared axes and legends.
#
# Author: Emma Van de Walle - Q-ForestLab
###############################################################################


library(patchwork)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(readxl)

# ------------------------------------------------------------
# Load sensor locations
# ------------------------------------------------------------

sensor_locations <- read_excel(
  "Data/sensor_locations_all_plots.xlsx"
) %>%
  mutate(
    plot_short = case_when(
      plot == "Cplot08" ~ "C8",
      plot == "Cplot13" ~ "C13",
      plot == "Rplot07" ~ "R7",
      plot == "Rplot11" ~ "R11",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(plot_short),

    # C8 and R7: TMS + LS, but exclude vertical_1 light sensor
    (plot_short %in% c("C8", "R7") &
       (sensor == "TMS" |
          (sensor == "LS" & SN != "vertical_1"))) |

      # C13 and R11: TMS only
      (plot_short %in% c("C13", "R11") &
         sensor == "TMS")
  )



# ------------------------------------------------------------
# Load TLS data
# ------------------------------------------------------------

tls_C8  <- readRDS("Data/TLS_scaled_DTM_and_grid_C8_Panama.rds")
tls_C13 <- readRDS("Data/TLS_scaled_DTM_and_grid_C13_Panama.rds")
tls_R7  <- readRDS("Data/TLS_scaled_DTM_and_grid_R7_Panama.rds")
tls_R11 <- readRDS("Data/TLS_scaled_DTM_and_grid_R11_Panama.rds")

# ------------------------------------------------------------
# Function to create canopy-height map with sensor locations
# ------------------------------------------------------------

create_canopy_height_plot <- function(grid,
                                      site_name,
                                      sensors,
                                      show_x = FALSE,
                                      show_y = FALSE,
                                      show_sensor_legend = FALSE) {

  canopy_height <- grid %>%
    filter(!is.na(density), density > 0) %>%
    group_by(X, Y) %>%
    summarise(
      canopy_height = max(Z),
      .groups = "drop"
    ) %>%
    complete(
      X = seq(min(grid$X), max(grid$X), by = 1),
      Y = seq(min(grid$Y), max(grid$Y), by = 1),
      fill = list(canopy_height = NA_real_)
    )

  sensors_plot <- sensors %>%
    filter(plot_short == site_name)

  p <- ggplot(
    canopy_height,
    aes(x = X, y = Y, fill = canopy_height)
  ) +
    geom_tile() +

    geom_point(
      data = sensors_plot,
      aes(
        x = x,
        y = y,
        shape = sensor
      ),
      inherit.aes = FALSE,
      size = 3.5,
      stroke = 1,
      colour = "black",
      fill = "white"
    ) +

    scale_shape_manual(
      values = c(
        "TMS" = 21,
        "LS"  = 24
      ),
      labels = c(
        "TMS" = "Air temperature sensor",
        "LS"  = "PAR sensor"
      ),
      name = NULL,
      guide = if (show_sensor_legend) {
        guide_legend(
          override.aes = list(
            size = 4,
            fill = "white",
            colour = "black"
          )
        )
      } else {
        "none"
      }
    ) +

    scale_fill_viridis_c(
      option = "viridis",
      direction = 1,
      limits = c(0, 40),
      na.value = "white",
      name = "Canopy height (m)",
      guide = guide_colorbar(
        barwidth = 1,
        barheight = 15,
        frame.colour = "black",
        ticks.colour = "black",
        title.position = "top",
        title.hjust = 0.5
      )
    ) +

    coord_fixed() +

    scale_x_continuous(
      limits = c(0, 60),
      breaks = seq(0, 60, 20),
      expand = c(0, 0)
    ) +

    scale_y_continuous(
      limits = c(0, 60),
      breaks = seq(0, 60, 20),
      expand = c(0, 0)
    ) +

    labs(
      title = site_name,
      x = NULL,
      y = NULL
    ) +

    theme_bw(base_size = 18) +

    theme(
      plot.title = element_text(
        size = 22,
        face = "bold",
        hjust = 0.5
      ),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      legend.title = element_text(
        size = 20,
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      legend.text = element_text(size = 20),
      legend.position = "right"
    )

  if (!show_x) {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  if (!show_y) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }

  return(p)
}


# ------------------------------------------------------------
# Create canopy height maps
# ------------------------------------------------------------

p1 <- create_canopy_height_plot(
  tls_C8$grid,
  "C8",
  sensor_locations,
  show_y = TRUE,
  show_sensor_legend = TRUE
)

p2 <- create_canopy_height_plot(
  tls_C13$grid,
  "C13",
  sensor_locations
)

p3 <- create_canopy_height_plot(
  tls_R7$grid,
  "R7",
  sensor_locations,
  show_x = TRUE,
  show_y = TRUE
)

p4 <- create_canopy_height_plot(
  tls_R11$grid,
  "R11",
  sensor_locations,
  show_x = TRUE
)

# ------------------------------------------------------------
# Combine the four panels
# ------------------------------------------------------------

p1 <- p1 + theme(plot.margin = margin(r = 20))
p3 <- p3 + theme(plot.margin = margin(r = 20))

final_plot <-
  (p1 + p2) /
  (p3 + p4) +
  plot_layout(
    guides = "collect"
  ) +
  plot_annotation(
    title = NULL,
    theme = theme(
      plot.title = element_text(
        size = 24,
        face = "bold",
        hjust = 0.5
      ),
      plot.margin = margin(
        t = 20,
        r = 20,
        b = 30,
        l = 10
      )
    )
  ) &
  theme(
    legend.position = "right",
    legend.box.spacing = unit(-100, "pt"),
    legend.title = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 18)
    ),
    legend.text = element_text(
      size = 20
    )
  )


# ------------------------------------------------------------
# Add common x- and y-axis titles
# ------------------------------------------------------------

final_plot_axes <- ggdraw() +
  draw_plot(final_plot) +

  draw_label(
    "Northing (m)",
    x = 0.09,
    y = 0.50,
    angle = 90,
    size = 22
  ) +

  draw_label(
    "Easting (m)",
    x = 0.44,
    y = 0.02,
    size = 22
  )


# Display figure
final_plot_axes


# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------

ggsave(
  "Output/Figure_canopy_height.png",
  final_plot_axes,
  width = 14,
  height = 10,
  dpi = 600
)

