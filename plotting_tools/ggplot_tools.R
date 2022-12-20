library(ggplot2)


`-.gg` <- function(plot, layer) {
  #' Add elements in background of existing ggplot
  #'
  #' @param plot        ggplot object
  #' @param layer       additional layer to ggplot
  #'
  #' @returns           ggplot object with the new layer added below the existing ones
  #' @export

  #Insert layer in ggplot below existing layers
  if (missing(layer)) {
    stop("Cannot use `-.gg()` with a single argument. Did you accidentally put - on a new line?")
  }
  if (!ggplot2::is.ggplot(plot)) {
    stop('Need a plot on the left side')
  }
  plot$layers = c(layer, plot$layers)
  plot
}