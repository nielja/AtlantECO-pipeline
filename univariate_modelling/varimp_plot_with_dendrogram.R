#Load all the deviance explained files for one type of plankton and make variable importance plots from this
library(ggplot2)
library(here)
library(plyr)
library(dplyr)
library(scales)
library(cowplot)
library(ggdendro)

varimp_plot_univariate_dendro <- function(df,
                                          df_path,
                                          predictors,
                                          remove_outliers = TRUE,
                                          d_nam = Sys.Date()){

  #' Function to analyse the results of the univaraite analysis
  #'
  #' This function shows and saves heatmaps of the deviance explained by each predictor as calculated by
  #' the univariate analyses. It also computes correlation between the predictors and shows these as a
  #' dendrogram.
  #'
  #' @param df              plankton dataframe
  #' @param df_path         path to the folder where the deviance explained files are stored
  #' @param predictors      vector of predictors
  #' @param remove_outliers Bool, if True, we load the datasets calculated without outliers
  #' @param d_nam           Date on wich the pdp files were stored, defaults to today
  #'
  #' @returns               returns and saves a plot in the df_path folder with deviance explained heatmap and a dendrogram


  #Set dataset choice
  outlier_name <- ifelse(remove_outliers == TRUE, "wo_outliers", "w_outliers")

  #Get list of files
  lf <- list.files(path = df_path)

  #Load all those relevant ones
  dev_files <- sort(lf[grepl(paste("_dev_expained_", sep = ""), lf, fixed = TRUE)])
  dev_df <- lapply(paste(df_path, dev_files, sep = "/"), read.csv) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(aggregation_level = rep(c("01°", "10°", "05°"), each = n()/3),
                    latitudinal = rep(c("latitudinal", "cell-wise"), each = n()/6, times = 3),
                    Version = paste(aggregation_level, latitudinal, sep = " - "))

  #Aspects for plotting
  max_dev <- max(dev_df$DevExplained)
  version_levels = sort(unique(dev_df$Version))

  #Make a dendrogram for the dataset
  #Correlations between predictors
  corloads = stats::cor(df[,predictors], use = "pairwise.complete.obs")
  #Replace all 1 by NaN
  corloads[corloads == 1] = NaN
  #dendrogram
  dissimilarity = 1 - abs(corloads)
  distance = as.dist(dissimilarity)
  dhc <- as.dendrogram(hclust(distance))

  # Rectangular lines
  ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
  ddata_labels = ggdendro::label(ddata)
  #Plot
  p_dendro <- ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend),(segment(ddata))) +
    ggplot2::coord_flip() +
    ggdendro::theme_dendro() +
    ggplot2::scale_x_continuous(expand = c(0, 0.5)) +
    ggplot2::scale_y_reverse(expand = c(0, 0)) +
    ggplot2::theme(plot.margin = unit(c(0, 0, 0, 0), "mm")) +

    ggplot2::geom_hline(yintercept = 0.3, color = "red", lty = "dashed")

  #Plot heatmap
  p <- dev_df %>%
    #Filter only relevnat predictors
    dplyr::filter(Variable %in% predictors) %>%
    dplyr::group_by(Version, Variable) %>%
    dplyr::summarise(DevExplained = mean(DevExplained)) %>%
    ggplot2::ggplot(ggplot2::aes(x = factor(Version, levels = version_levels,
                             labels = version_levels),
                  y = factor(Variable, levels = predictors,
                    labels = predictors), fill = DevExplained)) +
    ggplot2::geom_tile(color = "white") +
    #ggplot2::facet_wrap(Model~.) +
    ggplot2::scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlBu")),
                         breaks = seq(0, 1, by = 1/20),
    limits = c(0, max_dev), oob = scales::squish)  +
    ggplot2::theme_light() +
    ggplot2::scale_y_discrete(limits = ddata_labels$label, position = "right") +
    ggplot2::geom_text(ggplot2::aes(label = round(DevExplained, 2)), size = 2.5) +
    ggplot2::labs(x = "Dataset Version", y = "", fill = "Deviance\nExplained") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1),
          legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), "mm")) +
    ggplot2::scale_size_continuous(guide = "none")

  #Arrange them properly
  p_grid <- cowplot::plot_grid(p_dendro, p, align = "h", rel_widths = c(0.6, 1))

  #Save this plot
  ggplot2::ggsave(filename = paste0(df_path, "/varimp_plot_univariate_with_dendro_", outlier_name, Sys.Date(), ".png"),
         plot = p_grid, width = 200, height = 100, dpi = 300, unit = "mm")

  #Return the plot
  return(p_grid)
}



