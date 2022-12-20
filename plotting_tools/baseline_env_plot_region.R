library(maps)
library(ggplot2)

#Baseline environmental plot
baseline_plot_region <- function(region = "full"){
  #' Baseline map plot
  #'
  #' This function creates a baseline map plot to show the predictions on.
  #'
  #' @param region      character, "full", "SO", "NA", "NP" depending on whether we want a focus on one of the CPR regions
  #'
  #' @returns           ggplot basic map object of the chosen region
  #' @export

  #Plotting preparation
  world <- ggplot2::map_data("world")

  #Make baseline plot with map etc. -> make these before std and MESS plots
  if (region == "full"){
    p_base <- ggplot2::ggplot() +
      ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, group = group),
                   data = world[world$long <= 180,], fill = "grey85", colour = "grey50", size = 0.3) +
      ggplot2::scale_x_continuous(name = "", breaks = c(-180,-120,-60,0,60,120,180),
                         labels = c("180°","120°W","60°W","GM","60°E","120°E","180°"),
                         expand = c(0,0)) +
      ggplot2::scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
                         labels = c("-90°N","-60°N","-30°N","Eq","30°N","60°N","90°N"), expand = c(0,0)) +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),legend.key = ggplot2::element_rect(fill = "grey50"),
            panel.grid.major = ggplot2::element_line(colour = "white",linetype = "dashed")) +
      ggplot2::coord_quickmap() +
      ggplot2::theme(legend.position = "bottom",
            axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
  } else if (region == "SO"){
    #For polar projections
    # Defines the x axes required
    x_lines <- seq(-120,180, by = 60)
    #Create extra data.frames for lines
    lines_df <- data.frame(x_1 = 180,
                           y_1 = seq(-45, -85, by = -10),
                           label_1 = paste0(seq(-45, -85, by = -10), "°N"))

    xlines_df <- data.frame(x_lines = x_lines,
                            y = -32,
                            label =  c("120°W", "60°W", "0°", "60°E", "120°E", "180°W") )

    p_base <- ggplot2::ggplot() +
      ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, group = group),
                   data = world[(world$long <= 180) & (world$lat <= -40),], fill = "grey85", colour = "grey50", size = 0.3) +
      # Convert to polar coordinates
      ggplot2::coord_map("ortho", orientation = c(-90, 0, 0)) +
      ggplot2::scale_y_continuous(breaks = seq(-40, -90, by = -5), labels = NULL) +
      # Removes Axes and labels
      ggplot2::scale_x_continuous(breaks = NULL) +
      ggplot2::labs(x = "", y = "") +
      # Adds labels
      ggplot2::geom_text(ggplot2::aes(x = x_1, y = y_1, hjust = -0.2, label = label_1), data = lines_df, size = 2.5) +
      ggplot2::geom_text(ggplot2::aes(x = x_lines, y = y, label = label), data = xlines_df, size = 2.5) +
      # Adds axes
      ggplot2::geom_hline(ggplot2::aes(yintercept = -40), size = 1)  +
      ggplot2::geom_segment(ggplot2::aes(y = -40, yend = -90, x = x_lines, xend = x_lines), data = xlines_df, linetype = "dashed") +
      # Change theme to remove axes and ticks
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
            panel.grid.major = ggplot2::element_line(size = 0.25, linetype = 'dashed',
                                            colour = "black"),
            axis.ticks=ggplot2::element_blank(),
            legend.key = ggplot2::element_rect(fill = "grey50"),
            legend.position = "bottom")

  } else if (region == "NA"){
    #"NA_CPR", ifelse((((basin_name == "PacificOcean") | (basin_name == "ArcticOcean")) & ((Latitude > 40) & (Latitude < 65)) & ((Longitude < - 100) | (Longitude > 150))),
    #Make north Atlantic map
    p_base <- ggplot2::ggplot() +
      ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, group = group),
                   data = world,
                   fill = "grey85", colour = "grey50", size = 0.3) +
      ggplot2::coord_quickmap(xlim = c(-100, 30), ylim = c(40, 70)) +
      ggplot2::scale_x_continuous(name = "", breaks = c(-90,-60, -30, 0, 30),
                         labels = c("90°","60°W","30°W","GM","30°E")) +
      ggplot2::scale_y_continuous(name = "", breaks = seq(40, 70, by = 5),
                         labels = paste0(seq(40, 70, by = 5), "°N")) +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),legend.key = ggplot2::element_rect(fill = "grey50"),
            panel.grid.major = ggplot2::element_line(colour = "white",linetype = "dashed")) +
      ggplot2::theme(legend.position = "bottom",
            axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))

  } else if (region == "NP"){
    world_pacific <- world[(world$long <= -120) | (world$long >= 120), ]
    world_pacific$long <- ifelse(world_pacific$long > 0, world_pacific$long - 360, world_pacific$long)
    #Make north Pacific map
    p_base <- ggplot2::ggplot() +
      ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, group = group),
                   data = world_pacific,
                   fill = "grey85", colour = "grey50", size = 0.3) +
      ggplot2::coord_quickmap(xlim = c(-240, -120), ylim = c(38, 67)) +
      ggplot2::scale_x_continuous(name = "", breaks = c(-240, -210, -180, -150, -120),
                         labels = c("120°E", "150°E", "180°E/W","150°W","120°W")) +
      ggplot2::scale_y_continuous(name = "", breaks = seq(40, 65, by = 5),
                         labels = paste0(seq(40, 65, by = 5), "°N")) +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", size = 5),legend.key = ggplot2::element_rect(fill = "grey50"),
            panel.grid.major = ggplot2::element_line(colour = "white",linetype = "dashed")) +
      ggplot2::theme(legend.position = "bottom",
            axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  return(p_base)
}