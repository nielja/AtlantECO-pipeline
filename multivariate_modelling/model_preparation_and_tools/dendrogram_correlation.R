library(ggplot2)
library(ggdendro)
library(dplyr)
library(tidyr)

dendro_correlation <- function(model_setup,
                               full_predictor_set = FALSE,
                               dendrogram = FALSE){
  #' Function to return the maximum correlation between two variables in the chosen set as well as a dendrogram.
  #'
  #' @param model_setup         dataframe of model_setup parameters from h2o_model_pipeline
  #' @param full_predictor_set  Bool, if True, the dendrogram will be created with all possible environmental predictors
  #' @param dendrogram          Bool, if True, saves a dendrogram for the different datasets and variable choices
  #'
  #' @return                    model_setup with two new columns: maximum correlation between variables, and whether the model is okay to be run, also if chosen, saves a dendrogram


  # ====== Get unique combinations of dataframes and variable choices
  unique_df <- unique(model_setup$data_path)
  unique_combination <- unique(model_setup[,c("data_path", "predictors")])

  #Get all predictors used
  all_predictors <- unique(unlist(model_setup$predictors))

  if (full_predictor_set == TRUE){
    all_predictors <- c("Temperature_surface","Temperature_MLD",
                        "Salinity_surface", "Salinity_MLD",
                        "Nitrate_surface", "Nitrate_MLD", "Nitrate_MLD_gradient", "log_Nitrate_surface", "log_Nitrate_MLD",
                        "Phosphate_surface", "Phosphate_MLD", "Phosphate_MLD_gradient", "P_star_surface", "P_star_MLD",
                        "log_Phosphate_surface", "log_Phosphate_MLD", "log_P_star_surface", "log_P_star_MLD",
                        "Silicate_surface", "Silicate_MLD", "Si_star_surface", "Si_star_MLD",
                        "log_Silicate_surface", "log_Silicate_MLD", "log_Si_star_surface", "log_Si_star_MLD",
                        "Oxygen_surface", "Oxygen_MLD", "Oxygen_depth_200",
                        "MLD_MLD", "log_MLD","tAlk", "log_tAlk", "DIC", "pCO2", "RevelleFactor", "Omega_Ca", "Omega_Ar",
                        "Chl_a","log_Chl_a", "PAR", "PIC", "log_PIC", "BBP443", "log_BBP443", "z_eu", "Wind",
                        "EKE", "log_EKE", "Dist2Coast", "log_Dist2Coast", "Kd_490",  "log_Kd_490")
  }

  #Add a column to model_setup for max correlation and whether this is okay or not
  model_setup["max_correlation"] = NaN
  model_setup["correlation_below_threshold"] = TRUE
  model_setup["collapsed_predictors"] = paste(model_setup$predictors, sep = "_")

  # ====== Loop over dataframes ============
  for (d in unique_df){
    #Get dataset name and taxon name
    d_nam <- unique(model_setup[model_setup$data_path == d, "df_name"])
    tax_name <- unique(model_setup[model_setup$data_path == d, "tax_name"])
    #Load dataset
    df <- read.csv(paste(here::here(),d, sep = "/"))

    #Correlations between predictors
    corloads = cor(df[,all_predictors], use = "pairwise.complete.obs")
    #Replace all 1 by NaN
    corloads[corloads == 1] = NaN

    #For each row of the dataset, get maximum correlation for given parameters
    for (p in (model_setup[model_setup$data_path == d,]$predictors)){

      #Get maximum correlation
      if (length(p) == 1){
        max_cor = 0
        p_collapse = p
      } else {
        max_cor = max(abs(corloads[p, p]), na.rm = TRUE)
        #If p is more than one element, collapse it
        p_collapse = paste('c(\"', paste(p, collapse = '\", \"'), '\")', sep = "")
      }
      model_setup[(model_setup$data_path == d) & (model_setup$collapsed_predictors == p_collapse),"max_correlation"] = max_cor
    }


    #If a dendrogram should be produced, do so...
    if (dendrogram == TRUE){
      dissimilarity = 1 - abs(corloads)
      distance = stats::as.dist(dissimilarity)
      dhc <- stats::as.dendrogram(hclust(distance))

      # Rectangular lines
      ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
      ddata_labels = ggdendro::label(ddata)
      #Plot
      p_dendro <- ggplot2::ggplot() +
        ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend),(ggdendro::segment(ddata))) +
        ggplot2::coord_flip() +
        ggplot2::ylim(1.2, -1) +
        ggdendro::theme_dendro() +
        ggplot2::geom_text(ggplot2::aes(x = x, y = y, label = label, angle = 0, hjust = 0),
                  data= ddata_labels, size = 4) +
        ggplot2::geom_hline(yintercept = 0.3, color = "red", lty = "dashed") +
        ggtitle(paste(tax_name, d_nam))

      ggplot2::ggsave(paste(unique(model_setup$out_path), paste(tax_name, d_nam, Sys.Date(), ".png", sep = "_"),
                   sep = "/"), p_dendro, width = 200, height = 280, dpi = 300,
             unit = "mm")
    }

  }
  #All those rows where max_correlation > .7 are set to FALSE
  model_setup[model_setup$max_correlation > 0.7, "correlation_below_threshold"] = FALSE
  model_setup <- model_setup %>%
    dplyr::select(-collapsed_predictors)

  return(model_setup)
}
