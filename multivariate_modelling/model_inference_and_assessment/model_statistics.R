library(h2o)
library(scales)
library(hydroGOF)

model_statistics <- function(model_stat,
                             h2o_train,
                             h2o_test,
                             pdp_variation_sd = NA,
                             predictors = character(),
                             target = target,
                             report_importance = TRUE,
                             report_coefs = TRUE,
                             report_permutation = TRUE,
                             standardized = TRUE,
                             scaled_importance = FALSE,
                             percentage_importance = TRUE,
                             scaled_permutation = FALSE,
                             percentage_permutation = TRUE,
                             include.rmse.train = TRUE,
                             include.rmse.test = TRUE,
                             include.rsquared.train = TRUE,
                             include.rsquared.test = TRUE,
                             include.meanerror.train = TRUE,
                             include.meanerror.test = TRUE,
                             include.deviance = TRUE,
                             include.aic = TRUE,
                             include.devexp = TRUE,
                             include.NSE = TRUE,
                             include.nSD = TRUE,
                             include.pdpvariation = TRUE, ...){
  #' Model statistics calculation for one model
  #'
  #' This function computes a range of model statistics for one model. It is called by the wrapper function model_stats_list.
  #'
  #' @param model_stat              H2O Model to run the statistics analysis for
  #' @param h2o_train               H2O training dataframe
  #' @param h2o_test                H2O testing dataframe
  #' @param pdp_variation_sd        DataFrame with the PDP calculations, output of the pdp_plot_multi function
  #' @param predictors              vector of names of predictors to use
  #' @param target                  vector of name of target variable
  #' @param report_importance       Bool, if True, return variable importance
  #' @param report_coefs            Bool, if True, return the predictor coefficients
  #' @param report_permutation      Bool, if True, compute additional variable importance based on permutation
  #' @param standardized            Bool, if True, compute standardized predictor coefficients
  #' @param scaled_importance       Bool, if True, compute scaled importance (either this or percentage or none)
  #' @param percentage_importance   Bool, if True, compute percentage importance (either this or scaled_importance or none)
  #' @param scaled_permutation      Bool, if True, compute scaled permutation importance (either this or percentage_permutation or none)
  #' @param percentage_permutation  Bool, if True, compute percentage permutation importance (either this or scaled_permutation or none)
  #' @param include.rmse.train      Bool, if True, report RMSE on training data set
  #' @param include.rmse.test       Bool, if True, report RMSE on testing data set
  #' @param include.rsquared.train  Bool, if True, report R squared on training data set
  #' @param include.rsquared.test   Bool, if True, report R squared on testing data set
  #' @param include.meanerror.train Bool, if True, report mean absolute error on training data set
  #' @param include.meanerror.test  Bool, if True, report mean absolute error on testing data set
  #' @param include.deviance        Bool, if True, report deviance
  #' @param include.aic             Bool, if True, report AIC
  #' @param include.devexp          Bool, if True, report deviance explained
  #' @param include.NSE             Bool, if True, report Nash-Sutcliffe-efficiency
  #' @param include.nSD             Bool, if True, report normalized standard deviation between predictions
  #' @param include.pdpvariation    Bool, if True, report standard deviation between PDP curves of the different models
  #'
  #' @returns                       A data frame with statistics names, values and model name
  #' @export

  #Get model and model_id from model_stat
  model <- model_stat[[1]]
  model_name <- model_stat[[2]]
  model_id <- model_stat[[3]]

  #Performance on test set
  test_performance <- h2o::h2o.performance(model = model, newdata = h2o_test)

  ### ==== Model coefficients =====
  # If we have a GAM or a GLM with quadratic terms, add up the coefficients of the different smoothers
  if ((model@algorithm == "gam") | ((model@algorithm == "glm") & (length(model@model$coefficients_table$names) > length(predictors)+1))){
    coef_table <- model@model$coefficients_table
    coefs <- numeric()
    coefnames <- predictors
    #Loop over predictors and add up
    if (standardized == TRUE){
      for (p in predictors){
        ind_p <- unlist(lapply(coef_table$names, function(x) grepl(p, x, fixed = TRUE)))
        coefs <- c(coefs, sum(coef_table$standardized_coefficients[ind_p]))
      }
    } else {
      for (p in predictors){
        ind_p <- unlist(lapply(coef_table$names, function(x) grepl(p, x, fixed = TRUE)))
        coefs <- c(coefs, sum(coef_table$coefficients[ind_p]))
      }
    }
  } else {
    # for all other models, just extract coefnames from model
    coefnames <- predictors #model@model$names[1:length(model@model$names)-1]
    coef_table <- model@model$coefficients_table
    coefs <- unlist(lapply(predictors, function(x) coef_table[coef_table$names == x, "coefficients"]))
    if (is.null(coefs)){
      coefs <- rep(NaN, times = length(coefnames))
    }
  }

  ### ===== Variable importance ======
  importance_names <- paste(coefnames, "Importance")
  coefs_importance <- numeric()
  varimp_baseline <- numeric()
  #Again, check per model type
  if ((model@algorithm == "gam") | ((model@algorithm == "glm") & (length(model@model$coefficients_table$names) > length(predictors)+1))){
    varimp = h2o::h2o.varimp(model)
    for (p in predictors){
      #Baseline importance to then scale
      ind_p <- unlist(lapply(varimp$variable, function(x) grepl(p, x, fixed = TRUE)))
      varimp_baseline <- c(varimp_baseline, sum(varimp$relative_importance[ind_p]))}
    if (scaled_importance == TRUE) {
      coefs_importance <- rescale(varimp_baseline, to = 1, from = 0)
    } else if (percentage_importance == TRUE){
      for (p in predictors){
        ind_p <- unlist(lapply(varimp$variable, function(x) grepl(p, x, fixed = TRUE)))
        coefs_importance <- c(coefs_importance, sum(varimp$percentage[ind_p]))}
    } else {
      coefs_importance <- varimp_baseline
    }
  } else { #For all other models
    varimp <- h2o::h2o.varimp(model) %>% dplyr::arrange(match(variable, coefnames))
    if (scaled_importance == TRUE) {
      coefs_importance <- varimp$scaled_importance
    } else if (percentage_importance == TRUE){
      coefs_importance <- varimp$percentage
    } else {
      coefs_importance <- varimp$relative_importance
    }
  }

  ### ===== Permutation importance (Fisher-Yeast-algorithm on test set) =====
  permutation_names <- paste(coefnames, "Permutation Importance")
  coefs_permutation <- numeric()
  perm_baseline <- numeric()
  #Baseline RMSE
  base_RMSE <- h2o::h2o.performance(model, newdata = h2o_test)@metrics$RMSE
  if (model_name == "GLM" | model_name == "GAM"){
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
      predictor_importance[p0] = abs(h2o::h2o.performance(model, newdata = h2o_pred_glm)@metrics$RMSE - base_RMSE)

    }

    #Make proper shape
    importance_df <- data.frame(Names = names(predictor_importance),
                                perm_baseline = unname(predictor_importance)) %>%
      dplyr::mutate(percentage_permutation = perm_baseline/sum(perm_baseline))
    importance_df["scaled_permutation"] = scales::rescale(importance_df$perm_baseline, to = c(0, 1))

    if (scaled_permutation == TRUE) {
      coefs_permutation <- importance_df$scaled_permutation
    } else if (percentage_permutation == TRUE){
      coefs_permutation <- importance_df$percentage_permutation
    } else {
      coefs_permutation <- importance_df$perm_baseline
    }
  } else { #For all other models
    perm <- h2o::h2o.permutation_importance(model, newdata = h2o_test, metric = "RMSE", seed = 123) %>% dplyr::arrange(match(Variable, coefnames))
    if (scaled_permutation == TRUE) {
      coefs_permutation <- perm$`Scaled Importance`
    } else if (percentage_permutation == TRUE){
      coefs_permutation <- perm$Percentage
    } else {
      coefs_permutation <- perm$`Relative Importance`
    }
  }

  #Decide whether we report coefficients or importance value
  coefs_report.names <- character()
  coefs_report <- character()
  if (report_coefs == TRUE){
    coefs_report.names <- c(coefs_report.names, coefnames)
    coefs_report <- c(coefs_report, coefs)
  }
  if (report_importance == TRUE){
    coefs_report.names <- c(coefs_report.names, importance_names)
    coefs_report <- c(coefs_report, coefs_importance)
  }
  if (report_permutation == TRUE){
    coefs_report.names <- c(coefs_report.names, permutation_names)
    coefs_report <- c(coefs_report, coefs_permutation)
  }


  # create empty GOF vectors and subsequently add GOF statistics from model:
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.rmse.train == TRUE) {
    rmse_train <- model@model$training_metrics@metrics$RMSE
    gof <- c(gof, rmse_train)
    gof.names <- c(gof.names, "RMSE Train")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rmse.test == TRUE) {
    rmse_test <- test_performance@metrics$RMSE
    gof <- c(gof, rmse_test)
    gof.names <- c(gof.names, "RMSE Test")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rsquared.train == TRUE) {
    r2_train <- model@model$training_metrics@metrics$r2
    gof <- c(gof, r2_train)
    gof.names <- c(gof.names, "R^2 Train")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rsquared.test == TRUE) {
    r2_test <- test_performance@metrics$r2
    gof <- c(gof, r2_test)
    gof.names <- c(gof.names, "R^2 Test")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.meanerror.train == TRUE) {
    mae_train <- model@model$training_metrics@metrics$mae
    gof <- c(gof, mae_train)
    gof.names <- c(gof.names, "Mean Absolute Error Train")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.meanerror.test == TRUE) {
    mae_test <- test_performance@metrics$mae
    gof <- c(gof, mae_test)
    gof.names <- c(gof.names, "Mean Absolute Error Test")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    nulldev <- model@model$training_metrics@metrics$null_deviance
    if (is.null(nulldev)){nulldev = NaN}
    resdev <- model@model$training_metrics@metrics$residual_deviance
    if (is.null(resdev)){resdev = NaN}
    gof <- c(gof, nulldev, resdev)
    gof.names <- c(gof.names, "Null Deviance", "Residual Deviance")
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
  }
  if (include.devexp == TRUE) {
    devexp <- (nulldev - resdev) / nulldev
    gof <- c(gof, devexp)
    gof.names <- c(gof.names, "Deviance Explained")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.NSE == TRUE){
    #Make predictions on train and test set separately
    #For GLM, add polynomial terms
    if (model_name == "GLM"){
      #Adding polynomial terms
      h2o_pred_glm <- h2o_train
      predictors_glm <- predictors
      for (p in predictors){
        p_nam = paste(p, "^2", sep = "")
        predictors_glm <- append(predictors_glm, p_nam)
        h2o_pred_glm[p_nam] <- h2o_pred_glm[p]^2
      }
      #Prediction for GLM
      pr_model <- h2o::h2o.predict(object = model, newdata = h2o_pred_glm)
    } else {
      #For all other model types, regular prediction for H2O models
      pr_model <- h2o::h2o.predict(object = model, newdata = h2o_train)
    }

    pred_train <- as.data.frame(pr_model)$predict
    obs_train <- as.data.frame(h2o_train)[target][[1]]
    #Calculate NSE
    nse_train <- NSE(pred_train, obs_train)

    #Testing dataset
    #For GLM, add polynomial terms
    if (model_name == "GLM"){
      #Adding polynomial terms
      h2o_pred_glm <- h2o_test
      predictors_glm <- predictors
      for (p in predictors){
        p_nam = paste(p, "^2", sep = "")
        predictors_glm <- append(predictors_glm, p_nam)
        h2o_pred_glm[p_nam] <- h2o_pred_glm[p]^2
      }
      #Prediction for GLM
      pr_model <- h2o::h2o.predict(object = model, newdata = h2o_pred_glm)
    } else {
      #For all other model types, regular prediction for H2O models
      pr_model <- h2o::h2o.predict(object = model, newdata = h2o_test)
    }
    pred_test <- as.data.frame(pr_model)$predict
    obs_test <- as.data.frame(h2o_test)[target][[1]]

    #Calculate NSE
    nse_test <- NSE(pred_test, obs_test)
    gof <- c(gof, nse_train, nse_test)
    gof.names <- c(gof.names, "NSE_Train", "NSE_Test")
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
  }
  if (include.aic == TRUE) {
    aic <- model@model$training_metrics@metrics$AIC
    if (is.null(aic)){aic = NaN}
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  #If we want to include the average standard deviation of all pdp curves for all predictors
  if (include.pdpvariation == TRUE){
    gof <- c(gof, pdp_variation_sd)
    gof.names <- c(gof.names, "PDP_SD")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  # create data.frame
  model_stats_df <- data.frame(Values = t(c(model_id, model_name, as.numeric(coefs_report), gof))) %>%
    `colnames<-` (c("ModelID", "ModelName", coefs_report.names, gof.names))

  return(model_stats_df)
}

#Apply model_statistics for a list of models
model_stats_list <- function(model_list,
                             model_names,
                             model_id,
                             out_path,
                             predictors,
                             target,
                             df_choice,
                             remove_outliers,
                             weight_option,
                             region,
                             h2o_train,
                             h2o_test,
                              ave_nSD,
                              pdp_path = NA){
  #' Calculate model statistics for a list of models
  #'
  #' This is a wrapper function that calls the model_statistics function to assess a list of models.
  #'
  #' @param model_list       list of the models
  #' @param model_names      list of model names for plotting
  #' @param model_id         Integer model id, only relevant in case of multiple models trained at the same time
  #' @param out_path         path to folder where the plots should be saved
  #' @param predictors       vector of names of predictors to use
  #' @param target           vector of name of target variable
  #' @param df_choice        character name of the dataframe
  #' @param remove_outliers  Bool, whether outliers were removed in the modelling
  #' @param weight_option    character name of the column where the weights are stored in
  #' @param region           Character, "full", "SO", "NA", "NP" depending on whether we want a focus on one of the CPR regions
  #' @param h2o_train        H2O training dataframe
  #' @param h2o_test         H2O testing dataframe
  #' @param ave_nSD          data frame of the average standard deviation between model predictions (output of the prediction_map_plot function)
  #' @param pdp_path         path to the pdp dataframe (output of the pdp_plot_multi function)
  #'
  #' @returns                saves the model statistics calculated into a csv
  #'
  #' @export


  #Set names of model_list as model_id_list
  model_stat_list <- vector(mode = "list", length = length(model_list))
  #model_id_name <- paste(as.character(model_id), model_names, sep = "_")
  for (i in seq_along(model_list)){
    m <- model_list[[i]]
    m_nam <- model_names[i]
    m_id <- model_id
    model_stat_list[[i]] <-list(m, m_nam, m_id)

  }

  if (!is.na(pdp_path)){
    pdp <- read.csv(pdp_path)

    #Evaluate standard deviation
    pdp_variation_sd <- pdp %>%
     dplyr::group_by(Variable, VariableValue) %>%
     dplyr::summarise(model_sd = sd(ResponseValue, na.rm = TRUE)) %>%
      #Group by Variable
     dplyr::group_by(Variable) %>%
     dplyr::summarise(variable_sd_abs = mean(model_sd)) %>%
      #Get final mean relative sd
     dplyr::summarise(mean_sd_all_vars = mean(variable_sd_abs))

    pdp_variation_sd <- pdp_variation_sd[[1]]
  } else {pdp_variation_sd = NA}

  #List apply
  out_stats <- lapply(model_stat_list,
                      FUN = model_statistics,
                      h2o_test = h2o_test,
                      h2o_train = h2o_train,
                      predictors  = predictors,
                      target = target,
                      pdp_variation_sd = pdp_variation_sd
                      )

  #List of dfs into one df
  stat_df <- dplyr::bind_rows(out_stats) %>%
    #Add target, predictors, model_id, dataset_name, outliers, weight, weight_column  column
    dplyr::mutate(target = target,
           predictors = rep(as.character(paste(predictors, collapse = ", "), times = n())),
            model_id = model_id,
            df_choice = df_choice,
            remove_outliers = remove_outliers,
            weight_option = weight_option,
            region = region) %>%
    #Include average nSD
    dplyr::mutate(nSD_incl_MESS = ave_nSD[1],
           nSD_excl_MESS = ave_nSD[2])

  #Save this as a csv
  out_name <- paste(model_id, "_statistics_", Sys.Date(), ".csv", sep = "")
  out_name <- paste(out_path, out_name, sep = "/")
  write.csv(stat_df, out_name)
}
