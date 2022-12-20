library(tidyr)
library(h2o)
library(dplyr)
library(here)

#Set source for other code
source(here::here("General_Pipeline_PCA/Code/multivariate_modelling/model_training/rf_h2o_tuning_cleaning.R"))
source(here::here("General_Pipeline_PCA/Code/multivariate_modelling/model_training/gbm_h2o_tuning_cleaning.R"))
source(here::here("General_Pipeline_PCA/Code/multivariate_modelling/model_training/dl_h2o_tuning_cleaning.R"))
source(here::here("General_Pipeline_PCA/Code/multivariate_modelling/model_preparation_and_tools/h2o_clean.R"))

h2o_model_train <- function(predictors,
                            target,
                            out_path,
                            h2o_train,
                            h2o_test,
                            weights_column = "model_weights",
                            df_name = "",
                            model_id = 0,
                            n_models = 10,
                            tune = FALSE,
                            nfolds = 5,
                            operative_system = "linux"){
  #' Trains (and tunes) the five H2O models
  #'
  #' @param predictors          vector of names of predictors to use
  #' @param target              vector of name of target variable
  #' @param out_path            path to folder where the models should be saved
  #' @param h2o_train           H2O training dataframe
  #' @param h2o_test            H2O testing dataframe
  #' @param weights_column      character name of the column where the weights are stored in
  #' @param df_name             character name of the dataframe
  #' @param model_id            ID of the model, relevant only when training multiple models after each other
  #' @param n_models            Int number of models to be tuned at the same time -> this depends on your memory and has to be assessed with trial and error
  #' @param tune                Bool, whether to tune the models
  #' @param nfolds              Int, number of CV folds to use
  #' @param operative_system    Character, if linux, use xgboost, if windows, use gbm
  #'
  #' @returns                   Returns a named vector of the paths where the models are saved

  print(paste0("DF: ", df_name, "; Predictors: ", paste(predictors, collapse = ", "), "; Target: ", target))


  ### ===== Actual model setup ===========
  h2o::h2o.show_progress()

  #Get keys for the training and testing dataset
  train_test_keys <- h2o::h2o.ls()$key

  ### === Create GLM =======
  #Adding polynomial terms
  h2o_train_glm <- h2o_train
  predictors_glm <- predictors
  for (p in predictors){
    p_nam = paste(p, "^2", sep = "")
    predictors_glm <- append(predictors_glm, p_nam)
    h2o_train_glm[p_nam] <- h2o_train_glm[p]^2
  }
  print("GLM Model")
  glm_model <- h2o::h2o.glm(x = predictors_glm,
                       y = target[1],
                       training_frame = h2o_train_glm,
                       nfolds = nfolds,
                       weights_column = weights_column,
                       #fold_column = k_fold_column,
                       model_id = "GLM_Model")

  #Save model
  glm_path <- h2o::h2o.saveModel(object = glm_model,
                path = paste(out_path, paste(model_id, Sys.Date(), sep = "_"), sep = "/"),
                force = TRUE)

  rm_h2o_elements(train_test_keys)

  ### === Create GAM ====
  print("GAM Model")
  gam_model <- h2o::h2o.gam(y = target,
                       x = predictors,
                       training_frame = h2o_train,
                       gam_columns = predictors,
                       family = "gaussian",
                       model_id = "GAM_Model",
                       keep_gam_cols = FALSE,
                       weights_column = weights_column,
                        nfolds = nfolds)

  #Save model
  gam_path <- h2o::h2o.saveModel(object = gam_model,
                path = paste(out_path, paste(model_id, Sys.Date(), sep = "_"), sep = "/"),
                force = TRUE)

  rm_h2o_elements(train_test_keys)

  ### ==== Create tuned RF ======
  print("RF Model")
  rf_model <- rf_train_h2o(h2o_train = h2o_train,
                           target = target,
                           predictors = predictors,
                           train_test_keys = train_test_keys,
                           n_models = n_models,
                           out_path = out_path,
                           tune_model = tune,
                           nfolds = nfolds,
                           weights_column = weights_column,
  model_id_n = model_id)
                            #fold_column = k_fold_column)

  #Save model
  rf_path <- h2o::h2o.saveModel(object = rf_model,
                           path = paste(out_path, paste(model_id, Sys.Date(), sep = "_"), sep = "/"),
                            force = TRUE)

  rm_h2o_elements(train_test_keys)

  ### ==== Create tuned GBM  ======
  print("GBM Model")
  gbm_model <- gbm_train_h2o(h2o_train = h2o_train,
                             target = target,
                             predictors = predictors,
                             n_models = n_models,
                             train_test_keys = train_test_keys,
                             out_path = out_path,
                             tune_model = tune,
                             nfolds = nfolds,
                             weights_column = weights_column,
                              operative_system = operative_system,
  model_id_n = model_id)

  #Save model
  gbm_path <- h2o::h2o.saveModel(object = gbm_model,
                            path = paste(out_path, paste(model_id, Sys.Date(), sep = "_"), sep = "/"),
                force = TRUE)

  rm_h2o_elements(train_test_keys)

  ### ====== Create DL model ============================
  print("DL Model")
  dl_model <- dl_train_h2o(h2o_train,
                           target,
                           predictors,
                           train_test_keys = train_test_keys,
                           n_models = n_models,
                           out_path = out_path,
                           tune_model = tune,
                           nfolds = nfolds,
                           weights_column = weights_column,
                            model_id_n = model_id)
  #Save model
  dl_path <- h2o::h2o.saveModel(object = dl_model,
                           path = paste(out_path, paste(model_id, Sys.Date(), sep = "_"), sep = "/"),
                            force = TRUE)

  rm_h2o_elements(train_test_keys)

  ### ======= Return named list of paths where the models are saved ============
  model_paths <- c(glm_path, gam_path, rf_path, gbm_path, dl_path)
  names(model_paths) <- c("GLM", "GAM", "RF", "GBM", "DL")

  #Return this
  return(model_paths)
}