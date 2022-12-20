library(ggplot2)
library(cowplot)
library(gtable)

shift_legend <- function(p){
  #' Reposition legend into empty facet
  #'
  #' This function repositions the legend of a faceted ggplot into an empty facet
  #'
  #' @param p       faceted ggplot
  #'
  #' @returns       gtable object of the ggplot with shifted legend
  #' @export

  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
    gp <- ggplot2::ggplotGrob(p) # convert to grob
    } else {
    message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
    return(p)
    }
  } else {
    gp <- p
  }

  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }

  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
  max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")

  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable::gtable_add_grob(x = gp,
  grobs = gp[["grobs"]][[guide.grob]],
  t = empty.facet.panels[["t"]],
  l = empty.facet.panels[["l"]],
  b = empty.facet.panels[["b"]],
  r = empty.facet.panels[["r"]],
  name = "new-guide-box")

  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- cowplot::gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- cowplot::gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- cowplot::gtable_remove_grobs(gp, "guide-box")

return(gp)
}

