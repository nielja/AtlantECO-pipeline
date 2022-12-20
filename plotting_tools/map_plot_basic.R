library(ggplot2)
library(maps)
library(ggExtra)
library(dplyr)
library(scales)

map_plot <- function(df, target, target_name, grid_coords = FALSE,
                     cmap = "viridis", marginal = FALSE,
                     monthly = FALSE, alpha = NA,
                    legend.position.marginal = "bottom",
                     max_col = NaN){
  #' Create basic map plot
  #'
  #' This function creates a simple map plot for scatter data.
  #'
  #' @param df                          data frame with data
  #' @param target                      character name of target variable
  #' @param target_name                 character name of target variable to put in legend
  #' @param grid_coords                 Bool, if this is True, use Latitude and Longitude, otherwise lon_gridded and lat_gridded
  #' @param cmap                        Color map to be used for continuous target variables
  #' @param marginal                    Bool, whether to add marginal density histogram
  #' @param monthly                     Bool, if True, wrap the plots by monthly facet_wrap
  #' @param alpha                       Numeric, transparency value, if NA, will be set automatically based on the number of data points
  #' @param legend.position.marginal    position of the legend if marginal plot is shown
  #' @param max_col                     Numeric value of the maximum value of the colorbar, if NA, will be set automatically
  #'
  #' @returns                           Returns a ggplot showing the target scatter variable
  #' @export


  #Alpha as function of the number of points
  if (is.na(alpha)){
    alpha = 1000/(nrow(df))}
  if (alpha > 1){alpha = 1}

  #Get coordinate names
  if (grid_coords == TRUE){
    lon_name <- "lon_gridded"
    lat_name <- "lat_gridded"
  } else{
    lon_name <- "Longitude"
    lat_name <- "Latitude"
  }

  world <- ggplot2::map_data("world")
  #Plot
  p <- df %>%
    #sample_n(10000) %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x=get(lon_name), y=get(lat_name),
                   color = get(target)), alpha=alpha, size = 1) +
    ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, group = group),
                 data = world[world$long <= 180,], fill = "grey85", colour = "grey50", size = 0.3) +
    ggplot2::scale_x_continuous(name = "", breaks = c(-180,-120,-60,0,60,120,180),
                       labels = c("180°","120°W","60°W","GM","60°E","120°E","180°"), expand = c(0,0)) +
    ggplot2::scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
                       labels = c("-90°N","-60°N","-30°N","Eq","30°N","60°N","90°N"), expand = c(0,0)) +

    #scale_x_continuous(name = "", breaks = c(-180,-120,-60,0,60,120,180),
   #                   labels = c(bquote(180*degree),bquote(120*degree*W),bquote(60*degree*W),"GM",bquote(60*degree*E),bquote(120*degree*E),bquote(180*degree*E)), expand = c(0,0)) +
   #scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
   #                   labels = c(bquote(-90*degree*N),bquote(-60*degree*N),bquote(-30*degree*N),"Eq",bquote(30*degree*N),bquote(60*degree*N),bquote(90*degree*N)), expand = c(0,0)) +
   ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
          legend.key = ggplot2::element_rect(fill = "grey50"),
          panel.grid.major = ggplot2::element_line(colour = "white",
          linetype = "dashed"),
          legend.position = "bottom",
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::labs(color = target_name) +
    ggplot2::coord_quickmap()


  #Color map
  if (is.na(max_col)){
    max_col = quantile(df[target], 0.99, na.rm = TRUE)
  }

  #Add color map
  p <- p + ggplot2::scale_color_gradientn(colours = rev(hcl.colors(10, palette = cmap)),
                         limits = c(min(df[target], na.rm = TRUE),max_col)
                                                    ,
  oob = scales::squish)

  #Monthly facet wrap
  if (monthly == T){
    p <- p + ggplot2::facet_wrap(factor(Month)~., nrow = 4) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 4))
  }


  #Marginal density plot
  if (marginal == TRUE){
    p <- p + ggplot2::theme(legend.position = legend.position.marginal) +
      ggplot2::guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                    barheight = 0.7, barwidth = 10))
    if (categorical == TRUE){
      p <- ggExtra::ggMarginal(p, type = "density", groupFill = TRUE)

    } else {
      p <- ggExtra::ggMarginal(p, type = "density", fill = "grey")
    }
    #return(ggMarginal(p, type = "density", fill = "grey"))

  }

  return(p)

}
