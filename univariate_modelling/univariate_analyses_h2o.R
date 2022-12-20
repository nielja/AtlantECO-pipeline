#Run a univariate GLM and GAM to identify the most important predictors
library(tidyr)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidymv)
library(gridExtra)
library(here)
library(h2o)
library(GGally)
library(FactoMineR)
library(factoextra)

#Set source for other code
source(here::here("./General_Pipeline_PCA/Code/multivariate_modelling/model_inference_and_assessment/pdp_plot_h2o.R"))
source(here::here("./General_Pipeline_PCA/Code/plotting_tools/load_env_data.R"))
source(here::here("./General_Pipeline_PCA/Code/plotting_tools/geo_tools.R"))


outlier_iqr_col <- function(x) {
  #' Outlier removal function
  #'
  #' Remove outliers according to the z-score, but to apply to one column only
  #'
  #' @param x   values of which to identify outliers, usually a column of a dataframe
  #'
  #' @returns   returns the vector of values with outliers set to NA

  z <- (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
  x[z > 3] <- NA

  return(x)
}


#Function to run this
univariate_analysis_env_space <- function(df, aggregation_level, predictors, target,
                                          out_path = here::here("./General_Pipeline_PCA/Data/intermediate_data/"),
                                          pdp_calculation = TRUE,
                                          env_background = TRUE,
                                          latitudinal = FALSE,
                                          remove_outliers = TRUE,
                                          env_path =  here::here("./General_Pipeline_PCA/Data/environmental_climatologies/")){
  #' Function to run the full univariate analysis and create plots of environmental background coverage.
  #'
  #' This function runs a univariate analysis for a list of predictors at a given level of spatial aggregation.
  #' Two GLMs (one linear, one quadratic & linear) and a GAM are trained to predict the target variable
  #' as a function of each predictor variable. Additionally, the distribution of the environmental parameters
  #' in the plankton data set can be plotted in comparison to the distribution of the entire globe.
  #'
  #' @param df                      plankton dataFrame
  #' @param aggregation_level       Numeric, degree at which we want to aggregate (1°, 5°, 10°)
  #' @param predictors              list of predictors for which to run the analysis
  #' @param target                  characer name of target variable
  #' @param out_path                path to save the results to
  #' @param pdp_calculation         Bool, if True, partial dependence plots will be calculated
  #' @param env_background          Bool, if True, plots showing the distribution of environmental predictors in the
  #'                                plankton data and the entire globe will be made
  #' @param latitudinal             Bool, if True, conduct spatial aggregation over latitudinal bands, if False do it cell-wise
  #' @param remove_outliers         Bool, if True, remove outlier values based on the z-score > 3
  #' @param env_path                path to the environmental climatologies
  #'
  #' @returns                       Saves a number of partial depence plot calculations, deviance explained etc.
  #'                                The results can be visualized with the varimp_plot_univariate_dendro function.
  #' @export

  ### ========= Set basic characteristics ====================
  latitudinal_name <- ifelse(latitudinal == TRUE, "latitudinal_", "")
  outlier_name <- ifelse(remove_outliers == TRUE, "wo_outliers_", "w_outliers_")

  #### ===== Variable selection ===== ####
  #Load environmental data
  env_df <- load_env_data(predictor_list = predictors,
                          decimalCoords = FALSE,
                          env_path = env_path)

  #Apply aggregate function to environmental data
  env_df_agg <- aggregate_function(env_df, aggregation_level,
                                   predictors, latitudinal = latitudinal)

  #Pivot environmental data to long form
  env_df_all <- env_df_agg %>%
    tidyr::pivot_longer(cols = all_of(predictors),
                 names_to = "Variable", values_to = "VariableValue")


  #Apply aggregation function to plankton dataframes
  df_plankton_agg <- aggregate_function(df,
                             aggregation_level =  aggregation_level,
                             mean_cols = append(target, predictors),
                                        latitudinal = latitudinal)

  #Pivot to long form as well
  df_agg_all <- df_plankton_agg %>%
    tidyr::pivot_longer(cols = all_of(predictors),
                 names_to = "Variable", values_to = "VariableValue")

  #Check if there is a folder for UnivariateModels already, if not, create it
  if (ifelse(!dir.exists(file.path(here::here(), str_remove(out_path, here::here()))),
         dir.create(file.path(here::here(), str_remove(out_path, here::here()))), FALSE) == TRUE &
  ifelse(!dir.exists(file.path(out_path, "UnivariateModels")),
         dir.create(file.path(out_path, "UnivariateModels")), FALSE) == TRUE){
    print("output folders created")
  } else {"problem with the output folders"}


  ### === Calculate and save covariance matrix ====

  #Correlations between predictors
  cor_func <- function(df){
    #' Function to calculate correlation matrix for dataframe df for all predictors
    #'
    #' @param df        dataframe to calculate correlation matrix for
    #'
    #' @returns         correlation matrix for all predictors

    corloads <- stats::cor(df[, predictors], use = "pairwise.complete.obs")

    #Replace all 1 by NaN
    corloads[corloads == 1] <- NaN

    #To long form and append DataFrame name
    corloads <- corloads %>%
      as.data.frame() %>%
      dplyr::mutate(cor_from = rownames(.)) %>%
      tidyr::pivot_longer(cols = all_of(predictors), names_to = "cor_to", values_to = "cor_value") %>%
      dplyr::mutate(aggregation_level = aggregation_level)
    #Return corloads
    return(corloads)
  }
  #Apply this function to each dataframe
  corloads_df <- cor_func(df_plankton_agg)

  #Save this modified covariance matrix
  out_csv_name <- paste(out_path, paste(aggregation_level, "_covariance_matrix_",
                                         latitudinal_name, outlier_name,
                                        Sys.Date(), ".csv", sep =""), sep = "/")
  write.csv(corloads_df, out_csv_name)


  #Different types of models
  model_name_list <- c("Linear GLM", "Quadratic GLM", "GAM")

  if (pdp_calculation == TRUE){
    ### ======== Modelling =======================
    #Start H2O instance
    h2o::h2o.init()
    h2o::h2o.no_progress()

    #Set nbins for pdp analysis
    nbins = 25

    get_pdp_df <- function(p, df){
      #' Function to get the partial dependence dataframe for a variable p and a dataframe df.
      #'
      #' Trains a linear GLM, a quadratic GLM and a GAM and evaluates these at nbins points of the
      #' variable space.
      #'
      #' @param p       name of a predictor variable
      #' @param df      dataframe on which to calculate the pdp
      #'
      #' @returns       dataframe with deviance explained values for the given predictor


      print(paste(p))

      #Remove outliers
      if (remove_outliers == TRUE){
        df <- df %>%
          dplyr::select(all_of(c(all_of(target), p))) %>%
          dplyr::mutate_at(vars(target, p), list(outlier_iqr_col)) %>%
          tidyr::drop_na(any_of(c(target, p)))
      }

      #H2O training dataframe
      h2o_train <- h2o::as.h2o(df)

      #Train a linear GLM
      glm_linear <- h2o::h2o.glm(x = p,
                            y = target,
                            training = h2o_train,
                            nfolds = 5,
                            model_id = "GLM_linear")

      #Get deviance explained
      dev_exp_glm_l <- (glm_linear@model$cross_validation_metrics_summary["null_deviance", "mean"] -
        glm_linear@model$cross_validation_metrics_summary["residual_deviance", "mean"])/glm_linear@model$cross_validation_metrics_summary["null_deviance", "mean"]


      #Train a quadratic GLM
      h2o_train_glm <- h2o_train
      h2o_train_glm[paste(p, "^2", sep = "")] <- h2o_train_glm[p]^2
      p_glm <- c(p, paste(p, "^2", sep = ""))

      glm_quadratic <- h2o::h2o.glm(x = p_glm,
                               y = target,
                               training = h2o_train_glm,
                               nfolds = 5,
                               model_id = "GLM_quadratic")
      #Get deviance explained
      dev_exp_glm_q <- (glm_quadratic@model$cross_validation_metrics_summary["null_deviance", "mean"] -
        glm_quadratic@model$cross_validation_metrics_summary["residual_deviance", "mean"])/glm_quadratic@model$cross_validation_metrics_summary["null_deviance", "mean"]


      #Train a GAM
      gam_model <- h2o::h2o.gam(x = p,
                           y = target,
                           training = h2o_train,
                           nfolds = 5,
                           gam_columns = p,
                           model_id = "GAM_Model")
      #Get deviance explained
      dev_exp_gam <- (gam_model@model$cross_validation_metrics_summary[gam_model@model$cross_validation_metrics_summary$names == "null_deviance", "mean"] -
        gam_model@model$cross_validation_metrics_summary[gam_model@model$cross_validation_metrics_summary$names == "residual_deviance", "mean"])/
        gam_model@model$cross_validation_metrics_summary[gam_model@model$cross_validation_metrics_summary$names == "null_deviance", "mean"]

      #PDP plot for the three models
      p_plot <- pdp_plot_multi(model_list = c(glm_linear, glm_quadratic, gam_model),
                     model_name_list = model_name_list,
                     model_id = p,
                     predictors = p,
                     target = target,
                     predictor_name_list = p,
                     h2o_train = h2o_train,
                     point_density = FALSE,
                     rug = FALSE, nbins = nbins,
                     alpha = 0.9,
                     scatter = TRUE,
                     out_path = here::here(),
                     return_plot = TRUE,
                     return_df = TRUE)

      #Save this pdp dataframe along with the dataframe name and the deviance explained
      pdp_df <- data.frame(p_plot[[1]]) %>%
        dplyr::mutate(DevExplained = rep(c(dev_exp_glm_l, dev_exp_glm_q, dev_exp_gam), each = nbins))

      return(pdp_df)
    }

    #Run this function for each predictor
    pdp_out <- lapply(predictors, get_pdp_df, df = df_plankton_agg)
    #Join all datasets into one big one
    pdp_df <- dplyr::bind_rows(pdp_out)

    #Save the partial dependence relations as csv file
    out_csv_name <- paste(out_path, paste(aggregation_level, "_univariate_PDP_table_",
                                          latitudinal_name, outlier_name,
                                          Sys.Date(), ".csv", sep =""), sep = "/")
    write.csv(pdp_df, out_csv_name)

    #Plot and save the pdps with the curves
    if (remove_outliers == TRUE){
      df_agg_all <- df_agg_all %>%
        #Remove outliers
        dplyr::mutate_at(vars(target, VariableValue), list(outlier_iqr_col)) %>%
        tidyr::drop_na(any_of(c(target, "VariableValue")))
    }

    #Subdataset for deviance explained text
    plot_text <- df_agg_all %>%
      dplyr::group_by(Variable) %>%
      dplyr::summarise(x = max(VariableValue, na.rm = TRUE)) %>%
      merge(pdp_df) %>%
      dplyr::group_by(Variable, Model) %>%
      dplyr::summarise(DevExplained = mean(DevExplained),
                       x = mean(x)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(DevExpText = paste(round(DevExplained*100, 2), "%"),
             y = rep(c(max(df_agg_all[target], na.rm = TRUE),
                        mean(c(max(df_agg_all[target], na.rm = TRUE),quantile(df_agg_all[target][[1]], 0.75, na.rm = TRUE))),
                       quantile(df_agg_all[target][[1]], 0.75, na.rm = TRUE)), times = length(predictors)))

    #Actual plot
    p_out <- ggplot2::ggplot(ggplot2::aes(x = VariableValue, y = ResponseValue), data = pdp_df) +
      ggplot2::facet_wrap(.~factor(Variable, levels = predictors,
                                          labels = predictors), scale = "free_x",
      nrow = 2) +
      #Rug of training data
      ggplot2::geom_rug(ggplot2::aes(x = VariableValue, y = get(target)),
               data = df_agg_all, alpha = 0.05, outside = TRUE) +
      #Scatter points as well because why not
      ggplot2::geom_point(ggplot2::aes(x = VariableValue, y = get(target)), data = df_agg_all, alpha = 0.3) +
      ggplot2::theme_light() +
      ggplot2::coord_cartesian(clip = "off") +
      #PDP lines
      ggplot2::geom_line(ggplot2::aes(color = factor(Model,
                                   levels = model_name_list,
                                   labels = model_name_list)), lwd = 1) +
      #Deviance explained
      ggplot2::geom_label(ggplot2::aes(x = x, y = y, label = DevExpText,
                     color = factor(Model,
                                    levels = model_name_list,
                                    labels = model_name_list)),
                 data = plot_text, size = 3, hjust = 1, alpha = 0.8) +
      ggplot2::theme(strip.text.x= ggplot2::element_text(color ="black", angle = 0),
            strip.text.y = ggplot2::element_text(color = "black", angle = 0),
      legend.position = "bottom") +
      ggplot2::labs(x = "Variable Value", y = target, color = "Model")


    print(p_out)
    out_name <- paste(out_path, paste(aggregation_level, "_univariate_PDP_",
                                      latitudinal_name, outlier_name,
                                      Sys.Date(), ".png", sep =""), sep = "/")
    ggplot2::ggsave(out_name, p_out, width = 200, height = 150, dpi = 300, unit = "mm")

    #Save explained deviance as extra dataframe
    devexp <- pdp_df %>%
      dplyr::group_by(Model, Variable) %>%
      dplyr::summarise(DevExplained = mean(DevExplained))
    out_dev <- paste(out_path, paste(aggregation_level, "_dev_expained_",
                                     latitudinal_name, outlier_name,
                                     Sys.Date(), ".csv", sep =""),
                     sep = "/")
    write.csv(devexp, out_dev)

  }

  if (env_background == TRUE){
    ### ==== ENV BACKGROUND PLOTS ===================
    #Add all environmental data to our predictor dataset
    plot_pairs <- df_agg_all %>%
      dplyr::mutate(point_source = "Observation")
    plot_pairs <- env_df_all %>%
      dplyr::mutate(point_source = "Environmental background") %>%
      dplyr::bind_rows(plot_pairs)

    #Create colormap
    myColors <- c("gray40","orange")
    names(myColors) <- c("Environmental background", "Observation")
    colScale <- ggplot2::scale_colour_manual(name = "",values = myColors)
    colScaleFill <- ggplot2::scale_fill_manual(name = "", values = myColors)

    env_background_plot <- function(){
      g <- ggplot2::ggplot(ggplot2::aes(x = VariableValue, ..scaled..,
                      fill = point_source,
                      color = point_source), data = plot_pairs) +
        ggplot2::geom_density(alpha = 0.4, na.rm = TRUE) +
        ggplot2::facet_wrap(factor(Variable, levels = predictors,
                          labels = predictors)~., scale = "free_x") +
          colScale +
          colScaleFill +
        ggplot2::theme(legend.position = "bottom",
              strip.text = ggplot2::element_text(size = 6),
              axis.text = ggplot2::element_text(size = 5.5)) +
        ggplot2::labs(x = "Variable Value", y = "Scaled density")

      print(g)

      #Return plot
      return(g)
    }
    g <- env_background_plot()
    out_name <- paste(out_path, paste(aggregation_level, "_environmental_background_",
                                      latitudinal_name, outlier_name,
                                      Sys.Date(), ".png", sep =""), sep = "/")
    ggplot2::ggsave(out_name, g, width = 200, height = 100, dpi = 300, unit = "mm")

  }

}