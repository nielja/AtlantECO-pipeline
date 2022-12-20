library(h2o)
library(ggplot2)
library(here)
library(fgpt)


varimp_heatmap <- function(model_list, model_name_list,
                           model_id,
                           predictors,
                           predictor_name_list = list(),
                           h2o_test = data.frame(),
                           permutation= TRUE,
                           out_path = here::here("Plots")){
    #' Create a heatmap of variable importance
    #'
    #' This function creates a heatmap of variable importance across models (similar to h2o::h2o.varimp_heatmap)
    #' but allowing for GAMs with different variable names (due to the smoothers) and GLMs with polynomial terms
    #'
    #' @param model_list              list of the models
    #' @param model_name_list         list of model names for plotting
    #' @param model_id                Integer model id, only relevant in case of multiple models trained at the same time
    #' @param predictors              vector of names of predictors to use
    #' @param predictor_name_list     vector of bquote predictor names
    #' @param h2o_test                H2O testing dataframe
    #' @param permutation             Bool, if True, base the variable importance calculation on a permutation analysis
    #' @param out_path                Path to save the plots to
    #'
    #' @return                        Saves the variable importance plot
    #' @export

    #If no specific name list for the predictors is defined, assign the regular names
    #If no specific name list for the predictors is defined, assign the regular names
    if (length(predictor_name_list) == 0){
        predictor_names_df <- c("Temperature", "Chlorophyll-a", "MLD","EKE", "Oxygen")
        names(predictor_names_df) <-  c("Temperature_MLD", "log_Chl_a", "log_MLD",
                                        "log_EKE", "Oxygen_depth_200")
        if (all(predictors %in% names(predictor_names_df))){
            predictor_name_list <-  list()
            reordered_p <- list()
            for (p in names(predictor_names_df)){
                if (p %in% predictors){
                    predictor_name_list <- append(predictor_name_list,
                                                  predictor_names_df[p])
                    reordered_p <- append(reordered_p, p)
                }
            }
            predictor_name_list <- unlist(predictor_name_list)
            predictors <- unlist(reordered_p)
        } else {
            predictor_name_list = predictors
        }
    }
    names(predictor_name_list) <- c(1:length(predictors))

    #Create dataframe of variable importance values
    varimp_df <- data.frame(Predictor = rep(predictors, times = length(model_list)),
                            PredictorNames = rep(predictor_name_list,
                                                 times = length(model_list)),
                            Model = rep(model_name_list, each = length(predictors)),
                            Importance = NaN)


    #Loop over models to calculate importance
    for (n in seq_along(model_list)){
        m <- model_list[[n]]
        m_nam <- model_name_list[[n]]
        print(m_nam)

        #Baseline RMSE
        base_RMSE <- h2o::h2o.performance(m, newdata = h2o_test)@metrics$RMSE

        #Calculate the importance metric
        if (permutation == TRUE){
            if (m_nam == "GLM" | m_nam == "GAM"){
                #Do calculation manually
                predictor_importance <- vector(mode = "numeric", length = length(predictors))
                names(predictor_importance) <- predictors
                #Loop over predictors, shuffle them and add quadratic columns
                for (p0 in predictors){
                    print(p0)
                    h2o_pred_glm <- h2o_test
                    h2o_pred_glm[p0] <- h2o::as.h2o(fyshuffle(as.vector(h2o_pred_glm[p0])))
                    #Add quadratic columns
                    for (p in predictors){
                        p_nam = paste(p, "^2", sep = "")
                        h2o_pred_glm[p_nam] <- h2o_pred_glm[p]^2
                    }
                    #Calculate prediction RMSE
                    predictor_importance[p0] = abs(h2o::h2o.performance(m, newdata = h2o_pred_glm)@metrics$RMSE - base_RMSE)

                }

                #Make proper shape
                importance_df <- data.frame(Names = names(predictor_importance),
                                            vals_orig = unname(predictor_importance)) %>%
                  dplyr::mutate(Values = vals_orig/sum(vals_orig)) %>%
                  dplyr::select(-vals_orig)
            } else {
                importance_values <- h2o::h2o.permutation_importance(m,
                                                                newdata = h2o_test,
                                                                metric = "RMSE",
                                                                seed = 123)

                importance_df <- data.frame(Names = importance_values$Variable,
                                            Values = importance_values$Percentage)
            }



        } else { #Calculate regular importance values -> TODO: careful, this does not work for GLM with quadratic terms
            importance_values <- h2o::h2o.varimp(m)
            importance_df <- data.frame(Names = importance_values$variable,
                                       Values = importance_values$percentage)
        }

        # If this is a GLM or GAM, add up the values per predictor
        if ((m@algorithm == "gam") | ((m@algorithm == "glm") & (length(m@model$coefficients_table$names) > length(predictors)+1))){
            for (p in predictors){
                    ind_p <- unlist(lapply(importance_df$Names, function(x) grepl(p, x, fixed = TRUE)))
                    varimp_df[varimp_df$Predictor == p & varimp_df$Model == m_nam, "Importance"] <- sum(importance_df$Values[ind_p])
            }
        } else { #For all other models
            for (p in predictors){
                varimp_df[varimp_df$Predictor == p & varimp_df$Model == m_nam, "Importance"] = importance_df[importance_df$Names == p, "Values"]
            }
        }
    }

    fill_label <- ifelse(permutation == TRUE, "Variable\nImportance by\nPermutation",
                         "Variable\nImportance")

    #Now make a plot from this
    p <- ggplot2::ggplot(ggplot2::aes(x = factor(Model, levels = model_name_list,
      labels = model_name_list), y = factor(PredictorNames, levels = rev(predictor_name_list),
                                            labels = rev(predictor_name_list)), fill = Importance),
                 data = varimp_df) +
            ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlBu")),
                           breaks = seq(0, 1, by = 1/20))  +
      ggplot2::geom_text(ggplot2::aes(label = round(Importance, 2)), size = 15*varimp_df$Importance) +
      ggplot2::labs(x = "Model", y = "Environmental driver", fill = fill_label) +
      ggplot2::theme_light() +
      ggplot2::theme(axis.line.x = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            #axis.ticks = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank())

    #Save plot
    name_title <- ifelse(permutation == TRUE, "VariableImportancePermutation",
                         "VariableImportance")[1]
    out_name <- paste(model_id, name_title, Sys.Date(), ".png", sep = "_")
    save_path <- paste(out_path, out_name, sep = "/")
    height <- length(predictors)*20
    ggplot2::ggsave(save_path, p, width = 150, height = height, dpi = 300, unit = "mm")


}
