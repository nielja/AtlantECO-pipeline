library(h2o)
library(ggplot2)
library(dplyr)
library(here)
library(pryr)
library(envnames)

source(here::here("./General_Pipeline_PCA/Code/plotting_tools/ggplot_tools.R"))

pdp_plot_multi <- function(model_list,
                           model_name_list,
                           model_id,
                           predictors,
                           target,
                           target_name = "",
                           predictor_name_list = c(),
                           h2o_train = data.frame(),
                           nbins = 25,
                           rug = FALSE,
                           scatter = FALSE,
                           point_density = TRUE,
                           alpha = 0.7,
                           out_path = here::here("Plots"),
                           return_plot = FALSE,
                           return_df = FALSE){
    #' Create partial dependence plots for a list of models
    #'
    #' This function creates partial depence plots (PDP) for a list of models, it's similar to the h2o function
    #' but also works for GAMs with smoothing terms and GLMs with polynomial terms
    #'
    #' @param model_list            list of the models
    #' @param model_name_list       list of model names for plotting
    #' @param model_id              Integer model id, only relevant in case of multiple models trained at the same time
    #' @param predictors            vector of names of predictors to use
    #' @param target                vector of name of target variable
    #' @param target_name           character name of target variable to show in plot
    #' @param predictor_name_list   vector of bquote predictor names for the PDP plots
    #' @param h2o_train             H2O training dataframe
    #' @param nbins                 Integer number of bins to calculate pdp for
    #' @param rug                   Bool, if True, show a rug of observation values in the PDP
    #' @param scatter               Bool, if True, show scatter points of the observation data in the PDP
    #' @param point_density         Bool, if True, show a point density plot for the PDP
    #' @param alpha                 Numeric alpha value for point transparency
    #' @param out_path              Path to save the plots to
    #' @param return_plot           Bool, if True, returns the plot as well, just for diagnostic testing
    #' @param return_df             Bool, if True, returns the pdp dataframe, just for diagnostic testing
    #'
    #' @return                      Returns a list of 1) the path to the pdp dataframe or if return_df == TRUE 1) the pdp dataframe
    #'                              and depending on return_plot also 2) the plot element
    #' @export


  #Set target name if this isn't done yet
  if (target_name == ""){
    target_name = "target"
  }

  pdp_df <- data.frame(Model = character(),
                  Variable = character(),
                  VariableValue = numeric(),
                  ResponseValue = numeric()
                  )

  #Loop over models to calculate importance
  for (n in seq_along(model_list)){
      m <- model_list[[n]]
      m_nam <- model_name_list[[n]]
      #print(m_nam)

      #Create empty df for this model
      model_df <- data.frame(Variable = character(),
                              VariableValue = numeric(),
                              ResponseValue = numeric())

      #Calculate pdp values per variable
      if (m@algorithm != "glm" & m@algorithm != "gam"){
        pdp_vals <- h2o::h2o.partialPlot(m,
                        data = h2o_train,
                        nbins = nbins,
                        plot = FALSE)
      } else if (m@algorithm == "glm"){
        #print("GLM Model")
        #If there are no quadratic terms, use regular calculation as well
          if (length(m@model$coefficients_table$names[-1]) == length(predictors)){
              pdp_vals <- h2o::h2o.partialPlot(m,
                            data = h2o_train,
                            nbins = nbins,
                            plot = FALSE)
              pdp_vals <- list(pdp_vals)
          } else {
              #print("Quadratic terms, calculating pdp manually...")
              #Make calculation by hand
              pdp_vals = list()
              #Calculate the pdp manually
              #Get average values for all variables in training set
              mean_vals_train <- as.data.frame(h2o_train) %>%
                          tidyr::pivot_longer(all_of(predictors), values_to = "VariableValue", names_to = "Variable") %>%
                          dplyr::group_by(Variable) %>%
                          dplyr::summarise(VariableMean = as.numeric(mean(VariableValue, na.rm = TRUE)),
                          VariableMin = as.numeric(min(VariableValue, na.rm = TRUE)),
                          VariableMax = as.numeric(max(VariableValue, na.rm = TRUE)))
              if (nrow(mean_vals_train) == 1){
                mean_vals_train["Variable"] = predictors
              }

              #Loop over variables and create pdp for each
              for (p in predictors){
                  #print(p)
                  #Create predictor dataset
                  pred_vals <- seq(mean_vals_train[mean_vals_train$Variable == p, "VariableMin"][[1]],
                                  mean_vals_train[mean_vals_train$Variable == p, "VariableMax"][[1]],
                                  length.out = nbins)
                  glm_train_df <- data.frame(p = pred_vals)
                  colnames(glm_train_df) <- p
                  for (p_other in predictors[-which(predictors == p)]){
                    glm_train_df[p_other] <- rep(mean_vals_train[mean_vals_train$Variable == p_other, "VariableMean"][[1]],
                                              times = nbins)
                  }
                  #For all predictors, add quadratic terms
                  for (p in predictors){
                      p_squared_nam <- paste0(p, "^2")
                      glm_train_df[p_squared_nam] <- glm_train_df[p]^2
                  }

                  #Predict mean responses
                  glm_train_df["mean_response"] = as.data.frame(h2o::h2o.predict(object = m, newdata = h2o::as.h2o(glm_train_df)))$predict

                  pdp_vals <- append(pdp_vals, list(glm_train_df))
                  }
          }
      } else if (m@algorithm == "gam"){
        #print("GAM Model, calculating pdp manually...")
        pdp_vals = list()
        #Calculate the pdp manually
        #Get average values for all variables in training set
          mean_vals_train <- as.data.frame(h2o_train) %>%
                  tidyr::pivot_longer(all_of(predictors), values_to = "VariableValue", names_to = "Variable") %>%
                  dplyr::group_by(Variable) %>%
                  dplyr::summarise(VariableMean = as.numeric(mean(VariableValue, na.rm = TRUE)),
                            VariableMin = as.numeric(min(VariableValue, na.rm = TRUE)),
                            VariableMax = as.numeric(max(VariableValue, na.rm = TRUE)))
          if (nrow(mean_vals_train) == 1){
            mean_vals_train["Variable"] = predictors
          }

          #Loop over variables and create pdp for each
          for (p in predictors){
              #Create predictor dataset
              pred_vals <- seq(mean_vals_train[mean_vals_train$Variable == p, "VariableMin"][[1]],
                               mean_vals_train[mean_vals_train$Variable == p, "VariableMax"][[1]],
                               length.out = nbins)
              gam_train_df <- data.frame(p = pred_vals)
              colnames(gam_train_df) <- p
              for (p_other in predictors[-which(predictors == p)]){
                  gam_train_df[p_other] <- rep(mean_vals_train[mean_vals_train$Variable == p_other, "VariableMean"][[1]],
                                                times = nbins)
              }

              gam_train_df["mean_response"] = as.data.frame(h2o::h2o.predict(object = m, newdata = h2o::as.h2o(gam_train_df)))$predict

          pdp_vals <- append(pdp_vals, list(gam_train_df))
          }
      }

      #Loop over list elements of pdp_vals and assign those to our output dataset
      for (n_p in seq_along(pdp_vals)){
          p <- pdp_vals[[n_p]]
          #Get name of variable
          p_nam <- colnames(p)[1]

          #Assign this to the model's df
          model_df <- model_df %>%
                      add_row(Variable = p_nam,
                                VariableValue = unname(unlist(p[p_nam])),
                                ResponseValue = p$mean_response
                                )
          model_df["Model"] <- m_nam
      }

      #Assign this to the entire out dataframe
      pdp_df <- pdp_df %>%
                  add_row(model_df)
  }


  #Make a plot from this
  pdp_plot_df <- pdp_df %>%
    dplyr::mutate(VariableFactor = factor(Variable, levels = predictors, labels = c(1:length(predictors))))
  (p_out <- ggplot2::ggplot(ggplot2::aes(x = VariableValue, y = ResponseValue), data = pdp_plot_df) +
          ggplot2::geom_line(ggplot2::aes(color = factor(Model,
                        levels = model_name_list, labels = model_name_list))) +
          ggplot2::theme_light() +
          ggplot2::theme(strip.text.x = ggplot2::element_text(color ="black", angle = 0)) +
          ggplot2::labs(x = "Variable value", y = target_name, color = "Model"))

  #Calculate dataframe of training data to potentially add on
  scatter_df <- as.data.frame(h2o_train) %>%
        tidyr::pivot_longer(cols = all_of(predictors),
        names_to = "Variable",
        values_to = "VariableValue") %>%
    dplyr::mutate(VariableFactor = factor(Variable,
                                    levels = predictors,
                                   labels = c(1:length(predictors))))

  #IF we want to add point density heatmap
  if (point_density == T){
      p_out <- p_out -
        ggplot2::geom_hex(ggplot2::aes(x = VariableValue, y = get(target)),
                            data = scatter_df, alpha = alpha, size = 0.5) +
            ggplot2::scale_fill_distiller(palette = "Spectral",limits = c(0, 100),
                                 na.value = "red") +
            ggplot2::labs(fill = "Point\ndensity")
  }

  #If we want to add scatter of points
  if (scatter == T){
      p_out <- p_out - ggplot2::geom_point(ggplot2::aes(x = VariableValue, y = get(target)),
            data = scatter_df, alpha = 0.05, size = 0.5)
  }

  #If we want to add a rug of the observation values
  if (rug == T){
      p_out <- p_out + ggplot2::geom_rug(ggplot2::aes(x = VariableValue, y = get(target)),
            data = scatter_df, alpha = 0.02)
  }

  p_out <- p_out +
    ggplot2::facet_wrap(VariableFactor~.,
               labeller = label_bquote(.(as.expression(eval(str2expression(paste0('predictor_name_list', '$`',VariableFactor, '`')),
                                                              envir = environment())))),
               scale = "free_x")

  #Save the plot or return it or the pdp_df
  out_list <- list()
  if (return_df == TRUE){
    out_list <- append(out_list, list(pdp_df))
  } else {
    #Save the pdp_df
    out_name <- paste(model_id, "pdp_values", Sys.Date(), ".csv", sep = "_")
    save_path <- paste(out_path, out_name, sep = "/")
    write.csv(pdp_df, save_path, row.names = FALSE)
    #Return path to this file
    out_list <- append(out_list, save_path)
  }
  if (return_plot == FALSE){
    out_name <- paste(model_id, "pdp_plot", Sys.Date(), ".png", sep = "_")
    save_path <- paste(out_path, out_name, sep = "/")
    height <- ifelse(length(predictors) > 3, 150, 70)
    ggplot2::ggsave(save_path, p_out, width = 200, height = height, dpi = 300, unit = "mm")
  } else {out_list <- append(out_list, list(p_out))}

  if (length(out_list) != 0){
    return(out_list)
  }

}

