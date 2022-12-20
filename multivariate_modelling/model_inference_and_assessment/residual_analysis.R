library(plyr)
library(dplyr)
library(ggplot2)
library(here)
library(scales)
source(here::here("./General_Pipeline_PCA/Code/plotting_tools/baseline_env_plot_region.R"))
source(here::here("./General_Pipeline_PCA/Code/plotting_tools/shift_legend_ggplot.R"))
source(here::here("./General_Pipeline_PCA/Code/plotting_tools/geo_tools.R"))

residual_analysis <- function(h2o_train,
                              h2o_test,
                              predictors,
                              target,
                              model_id,
                              model_list,
                              model_name_list,
                              region,
                              out_path,
                              unit_small){
  #' Analyse model residuals
  #'
  #' This function analyses model residuals on the training and testing set in the geographic and environmental parameter
  #' space.
  #'
  #' @param h2o_train               H2O training dataframe
  #' @param h2o_test                H2O testing dataframe
  #' @param predictors              vector of names of predictors to use
  #' @param target                  vector of name of target variable
  #' @param model_id                Integer model id, only relevant in case of multiple models trained at the same time
  #' @param model_list              list of the models
  #' @param model_name_list         list of model names for plotting
  #' @param region                  Character, "full", "SO", "NA", "NP" depending on whether we want a focus on one of the CPR regions
  #' @param out_path                path to folder where the plots should be saved
  #' @param unit_small              bquote version of the unit the original target variable is in
  #'
  #' @returns                       Saves three plots, residuals vs predictions, map of residuals and residuals vs environmental variables
  #'
  #' @export



  #Merge the two to have full dataset to evaluate residuals on
  pred_df <- as.data.frame(h2o::h2o.rbind(h2o_train, h2o_test))
  #Just predictor dataset
  pred_df_predictors <- h2o::as.h2o(pred_df[predictors])
  #Loop over models and make predictions
  #Predict target variable for entire var space
  for (i in seq_along(model_list)){
    m <- model_list[i][[1]]
    m_nam <- model_name_list[i]

    #For GLM, add polynomial terms
    if (m_nam == "GLM"){
      #Adding polynomial terms
      h2o_pred_glm <- pred_df_predictors
      predictors_glm <- predictors
      for (p in predictors){
        p_nam = paste(p, "^2", sep = "")
        predictors_glm <- append(predictors_glm, p_nam)
        h2o_pred_glm[p_nam] <- h2o_pred_glm[p]^2
      }
      #Prediction for GLM
      pr_model <- h2o::h2o.predict(object = m, newdata = h2o_pred_glm)
    } else {
      #For all other model types, regular prediction for H2O models
      pr_model <- h2o::h2o.predict(object = m, newdata = pred_df_predictors)
    }
    #Append prediction column to env_df
    pred_df[m_nam] = as.data.frame(pr_model)$predict
  }

  #Pivot to long form
  pred_df <- pred_df %>%
    tidyr::pivot_longer(cols = model_name_list, names_to = "Model", values_to = "predictions") %>%
    dplyr::mutate(residuals = predictions - get(target))

  #Unit for axis labels
  unit_xlabel = bquote(Model~predictions~(.(unit_small)))
  unit_ylabel = bquote(Model~residuals~(.(unit_small)))

  #Plot residuals vs predictions for this
  resid_plot <- pred_df %>%
    ggplot2::ggplot(ggplot2::aes(x = predictions, y = residuals)) +
    ggplot2::geom_point(alpha = 0.05) +
    ggplot2::geom_smooth(method = "lm", color = "red") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
    ggplot2::facet_wrap(factor(Model, levels = model_name_list, labels = model_name_list)~.) +
    ggplot2::theme_light() +
    ggplot2::labs(x = unit_xlabel, y = unit_ylabel)

  #Plot residual map for yearly averages on a 5 degree scale
  pred_df_yearly <- pred_df %>%
    dplyr::mutate(lon_gridded = plyr::round_any(lon_gridded, 5),
           lat_gridded = plyr::round_any(lat_gridded, 5)) %>%
    dplyr::group_by(lon_gridded, lat_gridded, Model) %>%
    dplyr::summarise(residuals = mean(residuals, na.rm = TRUE))
  col_limit <- max(abs(c(quantile(pred_df_yearly$residuals, 0.01), quantile(pred_df_yearly$residuals, 0.99))))

  resid_map <- baseline_plot_region(region = region) +
    ggplot2::geom_tile(ggplot2::aes(x = lon_gridded, y = lat_gridded, fill = residuals), size = 1,
              data = pred_df_yearly) +
    ggplot2::scale_fill_gradient2() +
    ggplot2::facet_wrap(factor(Model, levels = model_name_list, labels = model_name_list)~., ncol = 2) +
    ggplot2::labs(fill = unit_ylabel) +
    ggplot2::guides(fill = guide_colorbar(title.position = "top", hjust = 0.5))

  #Shift legend into empty frame
  resid_map <- shift_legend(resid_map)

  #Residuals in enviromental space
  pred_df_env <- pred_df %>%
    tidyr::pivot_longer(cols = predictors, names_to = "PredictorName", values_to = "PredictorValue")

  #Plot for this
  env_resid_plot <- pred_df_env %>%
    ggplot2::ggplot(ggplot2::aes(x = PredictorValue, y = residuals)) +
    ggplot2::geom_point(alpha = 0.05) +
    ggplot2::geom_smooth(method = "lm", color = "red") +
    facet_grid(rows = vars(factor(Model, levels = model_name_list, labels = model_name_list)),
               cols = vars(factor(PredictorName, levels = predictors, labels = predictors)),
    scales = "free_x") +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0))

  #Create plot names
  resid_plot_name <- paste(out_path, paste0(model_id, "_residuals_vs_predictions_", Sys.Date(), ".png"), sep = "/")
  resid_map_name <- paste(out_path, paste0(model_id, "_residual_maps_", Sys.Date(), ".png"), sep = "/")
  resid_env_name <- paste(out_path, paste0(model_id, "_residual_env_space_", Sys.Date(), ".png"), sep = "/")

  #Save the plots
  ggplot2::ggsave(resid_plot_name, resid_plot, width = 200, height = 150, dpi = 300, unit = "mm")
  ggplot2::ggsave(resid_map_name, resid_map, width = 200, height = 150,  dpi = 300, unit = "mm")
  ggplot2::ggsave(resid_env_name, env_resid_plot, width = 200, dpi = 300, height = 200, unit = "mm")
}
