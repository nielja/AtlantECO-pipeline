library(h2o)
source(here::here("./General_Pipeline_PCA/Code/multivariate_modelling/model_preparation_and_tools/h2o_clean.R"))

dl_train_h2o <- function(h2o_train, target, predictors, train_test_keys,
                         out_path, model_id_n,n_models = 10, tune_model = FALSE,
                         nfolds = 0, model_id = "DL_Model", weights_column = "model_weights"){
  #' Function to train a deep neural network with the H2O package
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
  #'
  #' @returns A trained H2O Deep Neural Network model; if tuning is enabled, it also stores a csv of the tuning performance

  
  ### ==== H2o Model
  # create training & validation sets
  split <- h2o::h2o.splitFrame(h2o_train, ratios = 0.75, seed = 123)
  train <- split[[1]]
  valid <- split[[2]]
  
  #set the response column to target variable
  response <- target[1]
  
  # set the predictor names
  n_features <- length(predictors)
  
  #Default baseline DL Model with 5 hidden layers
  if (tune_model == FALSE){
  h2o_default <-  h2o::h2o.deeplearning(x = predictors,
                                 y = response,
                                 training_frame = train,
                                 validation_frame = valid,
                                 standardize = T,
                                 model_id = model_id,
                                 activation = "Rectifier",
                                 epochs = 100,
                                 seed = 123,
                                 nfolds = nfolds,
                                 hidden = 5,
                                 weights_column = weights_column,
                                 variable_importances = T)

    #Return the default model
    return(h2o_default)
  }
  if (tune_model == T){

    ### === random grid search ======
    # hyperparameter grid
    hyper_grid <- list(activation = c("Rectifier",
                                      "RectifierWithDropout",
                                      "Tanh",
                                      "Maxout",
                                      "MaxoutWithDropout"),
                          hidden = list(c(5, 5), c(10,10),c(15,15), c(20, 20),c(50,50,50)),
                          l1 = c(0,1e-3,1e-5),
                          l2 = c(0,1e-3,1e-5))


    # full grid search strategy
    search_criteria <- list(
      strategy = "Cartesian")

    #Calulate number of options
    expand_grid <- expand.grid(hyper_grid)
    n_options <-nrow(expand.grid(hyper_grid))


    #Get subset of hyper_grid
    print(paste("tuning", n_models, "at the same time to assess",
                n_options, "hyperparameter options"))
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

      #Perform grid search on this subset
      # perform grid search
      grid_search_params <- h2o::h2o.grid(
        algorithm = "deeplearning",
        grid_id = "dl_grid",
        x = predictors,
        y = response,
        training_frame = train,
        hyper_params = mini_grid_list,
        validation_frame = valid,
        weights_column = weights_column,
        nfolds = nfolds,
        seed = 123,
        parallelism = 0,
        search_criteria = search_criteria,
        epochs = 100
      )

      # collect the results and sort by our model performance metric
      # of choice
      grid_perf <- h2o::h2o.getGrid(
        grid_id = "dl_grid",
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
      dl_final <- grid_perf@model_ids[[1]]
      dl_final <- h2o::h2o.getModel(dl_final)

      #Save the best model and delete the grid
      h2o::h2o.saveModel(dl_final, paste0(out_path, "/", model_id_n, "_", Sys.Date(), "_dl_grid"))

      #Print length of elements in h2o::h2o.ls()
      print(length(h2o::h2o.ls()$key))
      #Delete everything from memory
      rm_h2o_elements(train_test_keys)
      print(length(h2o::h2o.ls()$key))
      #Clearing things
      rm(grid_perf)
      rm(grid_search_params)
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
              paste(out_path, paste(model_id_n, Sys.Date(), "DL_tuning_evaluation.csv", sep = "_"), sep = "/"))

    #Reload the final model
    final_model <- h2o::h2o.loadModel(paste0(out_path, "/", model_id_n, "_", Sys.Date(), "_dl_grid/", final_model_id[[1]]))


    #Return the final model
    return(final_model)
  }
}