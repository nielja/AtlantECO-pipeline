#RF tuning function to be called in multivariate model
library(dplyr)
#library(xgboost)
library(h2o)
source(here::here("./General_Pipeline_PCA/Code/multivariate_modelling/model_preparation_and_tools/h2o_clean.R"))

gbm_train_h2o <- function(h2o_train, target, predictors, out_path,model_id_n,
                          train_test_keys,
                          n_models = 10,
                          tune_model = FALSE,
                         nfolds = 0, model_id = "GBM_Model", weights_column = "model_weights",
                           operative_system = "linux"){
    #' Function to train a boosted regression tree with the H2O package
    #'
    #' @param h2o_train           H2O training dataframe
    #' @param target              vector of name of target variable
    #' @param predictors          vector of names of predictors to use
    #' @param out_path            path to folder where the models should be saved
    #' @param model_id_n          if "none", defaults to H2O given one, else will be set
    #' @param train_test_keys     list of h2o keys that just refer to the training and testing datasets, i.e. things that have to be kept
    #' @param n_models            Int number of models to be tuned at the same time -> this depends on your memory and has to be assessed with trial and error
    #' @param tune_model          Bool, whether to tune the model
    #' @param nfolds              Int, number of CV folds to use
    #' @param model_id            Character to give to the H2O ID
    #' @param weights_column      character name of the column where the weights are stored in
    #' @param operative_system    Character, if linux, use xgboost, if windows, use gbm
    #'
    #'
    #' @returns A trained H2O GBM model; if tuning is enabled, also saves a csv of the tuning performance



  #print("running gbm models...")
  # create training & validation sets
  split <- h2o::h2o.splitFrame(h2o_train, ratios = 0.75, seed = 123)
  train <- split[[1]]
  valid <- split[[2]]

  if (operative_system == "linux"){
    model_type = "xgboost"
  } else if (operative_system == "windows"){
    model_type = "gbm"
  }

  #Train default model for xgboost
  set.seed(123)
  #Note that this doesn't run on windows
  if (tune_model == FALSE){
    # training basic GBM model with defaults
    if (model_type == "xgboost"){
      h2o_default <- h2o::h2o.xgboost(
        x = predictors,
        y = target,
        training_frame = train,
        validation_frame = valid,
        nfolds = nfolds,
        #ntrees = 5000,
        score_each_iteration = TRUE,
        stopping_rounds = 10, #if there is no improvement after 10 rounds, stop
        stopping_tolerance = 0,
        seed = 123,
        model_id = model_id,
        weights_column = weights_column
      )
    } else if (model_type == "gbm"){
      h2o_default <- h2o::h2o.gbm(
        x = predictors,
        y = target,
        training_frame = train,
        validation_frame = valid,
        nfolds = nfolds,
        #ntrees = 5000,
        score_each_iteration = TRUE,
        stopping_rounds = 10, #if there is no improvement after 10 rounds, stop
        stopping_tolerance = 0,
        seed = 123,
        model_id = model_id,
        weights_column = weights_column
      )
    }

    #Return the default model
    return(h2o_default)

  } else if (tune_model == T){
    #print("tuning the GBM models...")
    #H2O grid search for parameters

    # create hyperparameter grid
    hyper_grid <- list(
      max_depth = c(1, 3, 5),
      min_rows = c(1, 5, 10),
      learn_rate = c(0.01, 0.05, 0.1),
      #learn_rate_annealing = c(.99, 1),
      sample_rate = c(.5, .75, 1),
      col_sample_rate = c(.8, .9, 1)
    )


    # full grid search strategy
    search_criteria <- list(
      strategy = "Cartesian")

    #Calulate number of options
    expand_grid <- expand.grid(hyper_grid)  %>%
      dplyr::arrange(desc(max_depth), min_rows, learn_rate, desc(sample_rate), desc(col_sample_rate))
    n_options <-nrow(expand.grid(hyper_grid))


    #Get subset of hyper_grid
    print(paste("tuning", n_models, "at the same time to assess", n_options, "hyperparameter options"))
    for (i in c(1:(n_options/n_models))){
      print(i)
      #Selection of rows of hyperparameter grid
      mini_grid <- expand_grid[(((i-1)*n_models+1):(i*n_models)), ]

      #Make this to a list in the form of hyper_grid
      mini_grid_list <- vector(mode = "list", length = ncol(mini_grid))
      names(mini_grid_list) = colnames(mini_grid)
      for (c in colnames(mini_grid)){
        mini_grid_list[c] <- unique(mini_grid[c])
      }

      # perform grid search
      gbm_grid <- h2o::h2o.grid(
        algorithm = model_type,
        grid_id = "gbm_grid",
        x = predictors,
        y = target,
        nfolds = nfolds,
        training_frame = train,
        validation_frame = valid,
        weights_column = weights_column,
        hyper_params = mini_grid_list,
        search_criteria = search_criteria, # add search criteria
        score_each_iteration = TRUE,
        parallelism = 0,
        ntrees = 500,
        seed = 123
      )

      # collect the results and sort by our model performance metric
      # of choice
      grid_perf <- h2o::h2o.getGrid(
        grid_id = "gbm_grid",
        sort_by = "RMSE",
        decreasing = FALSE
      )

      #Keep the hyperparameter evaluation grid and append it to earlier rounds
      if (i == 1){
        hyper_param_grid <- as.data.frame(grid_perf@summary_table)
      } else {
        hyper_param_grid <- rbind(hyper_param_grid, as.data.frame(grid_perf@summary_table))
        hyper_param_grid <- hyper_param_grid %>% dplyr::arrange(rmse)
      }

      #Get model id of the optimal model and save this
      gbm_final <- grid_perf@model_ids[[1]]
      gbm_final <- h2o::h2o.getModel(gbm_final)

      #Save the best model and delete the grid
      h2o::h2o.saveModel(gbm_final, paste0(out_path, "/", model_id_n, "_", Sys.Date(), "_gbm_grid"))

      #Print length of elements in h2o::h2o.ls()
      print(length(h2o::h2o.ls()$key))
      #Delete everything from memory
      rm_h2o_elements(train_test_keys)
      print(length(h2o::h2o.ls()$key))
      #Clearing things
      rm(grid_perf)
      rm(gbm_grid)
      #Clear everything
      gc()
      gc()
      # tell back-end cluster nodes to do three back-to-back JVM full GCs.
      h2o:::.h2o.garbageCollect()
      h2o:::.h2o.garbageCollect()
      h2o:::.h2o.garbageCollect()

    }

    #In the end, evaluate everything, take the best model, reload this and return it
    final_model_id <- hyper_param_grid %>%
      dplyr::filter(rmse == min(rmse)) %>%
      dplyr::select(model_ids)


    if (nrow(final_model_id) > 1){
      #If we have more than one final model, choose the first one but print a warning
      print(paste("more than one optimal model with same rmse: ", final_model_id))
      print("choosing the first one, but remember this")
      final_model_id <- final_model_id[1,]
    }



    #Save the tuning evaluation
    write.csv(hyper_param_grid,
              paste(out_path, paste(model_id_n, Sys.Date(), "GBM_tuning_evaluation.csv", sep = "_"), sep = "/"))

    #Reload the final model
    final_model <- h2o::h2o.loadModel(paste0(out_path, "/", model_id_n, "_", Sys.Date(), "_gbm_grid/", final_model_id[[1]]))

    #Return final model
    return(final_model)

    }
 }