library(here)
library(parallel)
library(ggplot2)
library(h2o)
source(here::here("General_Pipeline_PCA/Code/multivariate_modelling/model_preparation_and_tools/preprocess_dataset.R"))
source(here::here("General_Pipeline_PCA/Code/multivariate_modelling/model_training/models_h2o.R"))
source(here::here("General_Pipeline_PCA/Code/multivariate_modelling/model_inference_and_assessment/assess_model_performance_h2o.R"))
source(here::here("General_Pipeline_PCA/Code/multivariate_modelling/model_preparation_and_tools/dendrogram_correlation.R"))
source(here::here("General_Pipeline_PCA/Code/multivariate_modelling/model_inference_and_assessment/uncertainty_analysis.R"))
source(here::here("General_Pipeline_PCA/Code/plotting_tools/baseline_env_plot_region.R"))
source(here::here("General_Pipeline_PCA/Code/multivariate_modelling/model_preparation_and_tools/predictor_name_func.R"))


model_pipeline <- function(dataset_choices,
                           predictor_choices,
                           target,
                           tune_model = FALSE,
                           evaluation_mode = FALSE,
                           model_paths = c(),
                           model_name_list = c("GLM", "GAM", "RF", "GBM", "DL"),
                           unit_small = bquote('#'),
                           unit_large = bquote('#'),
                           unit_conversion = 1,
                           monthly = FALSE,
                           seasonal = TRUE,
                           n_models = 15,
                           pca = FALSE,
                           remove_highly_correlated = TRUE,
                           max_color = "",
                           #Variable importance heatmap
                           varimp_heatmap_func = TRUE,
                           #PDP plot
                           pdp_plot_func = TRUE,
                           #Model statistics as csv
                           model_statistics_func = TRUE,
                           #Prediction map
                           prediction_map_func = TRUE,
                           external_field_path = "",
                            #Whether to use the layouted predictor names
                           use_pred_names = TRUE,
                           #Which operating system to use
                           operative_system = "linux"){

  #' Pipeline function to preprocess the data, train multivariate models and assess them
  #'
  #' This function conducts all steps of the multivariate modelling process. The dataset is preprocessed, i.e.
  #' split into training and testing and five models (GLM, GAM, RF, GBM, DL) are trained. The models are then
  #' evaluated using a variable importance analysis, PDP curves and calculating performance statistics and model
  #' residuals. Additionally, they are used to extrapolate the target variable to the entire region of interest.
  #'
  #' @param dataset_chocies         Dataframe with columns data_path, out_path, df_name, region, weights_column, weights_option
  #' @param predictor_chocies       List with vectors of different predictor combinations
  #' @param tune_model              Bool, whether to tune the RF, GBM and DL, otherwise they are just fitted
  #' @param evaluation_mode         Bool, if FALSE (default), actually trains the models, otherwise loads existing models and evaluates them; note this only works for one setup at a time
  #' @param model_paths             list of paths to where the models are stored, only needs to be specified if evaluation_mode = TRUE
  #' @param model_name_list         list of the model names, don't modify this
  #' @param unit_small              bquote version of the unit the original target variable is in
  #' @param unit_large              bquote version of an aggregate version of the small unit to mention for the annual aggregate
  #' @param unit_conversion         numerical value giving the conversion factor between unit_small and unit_large
  #' @param monthly                 Bool, if True, create monthly plots of the extrapolated prediction
  #' @param seasonal                Bool, if True, create seasonally averaged plots of the extrapolated prediction
  #' @param n_models                Int number of models to be tuned at the same time -> this depends on your memory and has to be assessed with trial and error
  #' @param pca                     Bool, if True, compute PCA-transformation of the environmental predictors and train the models on these
  #' @param remove_highly_correlated Bool, if True, remove predictor sets with maximum covariance > 0.7
  #' @param max_color               numeric value for the maximum color of the color bar in the prediction plot, if "" it will be calculated automatically
  #' @param varimp_heatmap_func     Bool, if True, the variable importance heatmap is created and saved
  #' @param pdp_plot_func           Bool, if True, partial dependence plots are created and saved
  #' @param model_statistics_func   Bool, if True, model performance statistics and annual aggregations of the target variable are computed and saved
  #' @param prediction_map_func     Bool, if True, extrapolated maps of the target variable are computed and saved
  #' @param external_field_path     Character string of path to external field data set, if "", use WOA fields
  #' @param use_pred_names          Bool, if True, we use fancy layouting for the variable names in the PDP, not implemented for all variables yet
  #'
  #' @returns Returns a trained models and a ton of plots and csv files depending on what options are enabled
  #'
  #' @export


  #Number of dataset choices
  dataset_options <- nrow(dataset_choices)
  target_options <- length(target)
  #Get overview
  print(paste("Number of models to train:",
              dataset_options*length(predictor_choices)*length(target)*5))

  #Add predictor choices and target variable to the dataset
  model_setup <- dataset_choices %>%
    #Add target
    dplyr::slice(rep(1:n(), each = length(target))) %>%
    dplyr::mutate(target = rep(target, times = dataset_options)) %>%
    #Add predictors
    dplyr::slice(rep(1:n(), each = length(predictor_choices))) %>%
    dplyr::mutate(predictors = rep(predictor_choices, times = dataset_options * target_options),
                  model_id = seq(1:n()))

  #Check correlation between predictors and if wanted, save a dendrogram
  ###TODO: add PCA option
  if (pca == FALSE){
    model_setup_check <- dendro_correlation(model_setup,
                                            full_predictor_set = FALSE,
                                            dendrogram = FALSE)

    #We don't want to run the model for wrong choices of predictors, i.e. give a warning
    # and potentially drop those rows with correlations that are too high
    if (any(model_setup_check$correlation_below_threshold == FALSE)){
      print(paste("Correlations too high for predictor set", model_setup_check[model_setup_check$correlation_below_threshold == FALSE, "predictors"]))
      print(paste("Maximum correlation value ", model_setup_check[model_setup_check$correlation_below_threshold == FALSE, "max_correlation"]))

    if (remove_highly_correlated == TRUE){
    model_setup <- model_setup_check %>%
      dplyr::filter(correlation_below_threshold == TRUE)}
    } else {model_setup <- model_setup_check}
  }
  if (evaluation_mode == FALSE){
    print(paste("After corrections, training", nrow(model_setup)*5, "models"))
  } else {print("Evaluating models, no training")}

  #Start H2O session
  h2o::h2o.show_progress()
  h2o::h2o.init(enable_assertions = FALSE, ignore_config = TRUE)
  h2o::h2o.removeAll()

  ### ====== Model training and evaluation loop ========
  #Loop through rows of the model_setup dataset, train models for each and evaluate them
  system.time(
  for (r in row.names(model_setup)){ #needs approximately 1.5 mins per row
    print(paste(r, "/", nrow(model_setup)))
    #Select that row from the dataframe
    model_choices <- model_setup[r,]

   #Load and preprocess the dataset
    df_preprocessed <- preprocessing_modelling(predictors = model_choices$predictors[[1]],
                        target = c(model_choices$target[[1]]),
                        data_path = model_choices$data_path[[1]],
                        pca = pca,
                        remove_outliers = TRUE,
                        region = model_choices$region[[1]],
                        weights_column = model_choices$weights_column[[1]],
                        weight_option = model_choices$weight_options[[1]])
    h2o_train <- df_preprocessed[[1]]
    h2o_test <- df_preprocessed[[2]]
    #Also get pca element and eventually pca-transformed variables
    res.pca <- df_preprocessed[[3]]
    predictors <- df_preprocessed[[4]]

    #Reorder predictors and get proper names
    predictor_preprocessing <- pred_name_func(predictors,
                                              use_pred_names = use_pred_names)
    predictors <- predictor_preprocessing[[1]]
    predictor_name_list <- predictor_preprocessing[[2]]
    assign("predictor_name_list", predictor_name_list, envir = .GlobalEnv)

    #Train the models
    if (evaluation_mode == FALSE){
      model_paths <- h2o_model_train(predictors = predictors,
                                    target = c(model_choices$target[[1]]),
                                    out_path = model_choices$out_path[[1]],
                                    h2o_train = h2o_train,
                                    h2o_test = h2o_test,
                                    weights_column = "model_weights",
                                    df_name = model_choices$df_name[[1]],
                                    model_id = model_choices$model_id[[1]],
                                    tune = tune_model,
                                    n_models = n_models,
                                    nfolds = 5,
                                    operative_system = operative_system)

      model_name_list <- names(model_paths)
    }

    #Assess the models
    assess_model_performance_h2o(model_paths = model_paths,
                                 model_name_list = model_name_list,
                                 model_id = model_choices$model_id[[1]],
                                 predictors = predictors,
                                 model_predictors = model_choices$predictors[[1]],
                                 target = model_choices$target[[1]],
                                 h2o_train = h2o_train,
                                 h2o_test = h2o_test,
                                 data_path = model_choices$data_path[[1]],
                                 out_path = model_choices$out_path[[1]],
                                 unit_small = unit_small,
                                 unit_large = unit_large,
                                 unit_conversion = unit_conversion,
                                 df_choice = model_choices$df_name[[1]],
                                 remove_outliers = TRUE,
                                 weight_option = model_choices$weight_options[[1]],
                                 region = model_choices$region[[1]],
                                 predictor_name_list = predictor_name_list,
                                 #Variable importance heatmap
                                 varimp_heatmap_func = varimp_heatmap_func,
                                 varimp_permutation = TRUE,
                                 #PDP plot
                                 pca = pca,
                                 res.pca = res.pca,
                                 pdp_plot_func = pdp_plot_func,
                                 pdp_point_density = FALSE,
                                 pdp_rug = TRUE,
                                 pdp_scatter = TRUE,
                                 #Model statistics as csv
                                 model_statistics_func = model_statistics_func,
                                 #Prediction map
                                 prediction_map_func = prediction_map_func,
                                 env_path = here::here("./Data/21_10_18_environmental_data/all_env_vars_all_months/"),
                                 external_field_path = external_field_path,
                                 prediction_monthly = monthly,
                                 prediction_seasonally = seasonal,
                                 prediction_MESS = TRUE,
                                 prediction_mess_plot = TRUE,
                                 prediction_std_plot = TRUE,
                                 prediction_save_plot = TRUE,
                                 prediction_region = model_choices$region[[1]],
                                 prediction_min_col = 0,
                                 prediction_max_col = max_color,

                                 prediction_add_points = TRUE,
                                 residual_analysis_func = TRUE)

    #Remove all elements from H2O server to keep it at full capacity
    h2o::h2o.removeAll()

  }
  )

  if (prediction_map_func == TRUE){
    #Check if we have a different subfolder for the external environmental predictor field
    if ((external_field_path != "") & (dir.exists(paste0(dataset_choices$out_path, "/",
                                                         stringr::str_extract(basename(external_field_path), "[^.]+"))))){
      out_path <- paste0(dataset_choices$out_path, "/", stringr::str_extract(basename(external_field_path), "[^.]+"))
    } else {out_path <- unique(dataset_choices$out_path)}
    #Gather flux and statistics dataframes and create joint ones
    flux_base <- paste("_total_stocks_", Sys.Date(), ".csv", sep = "")
    flux_df <- vector(mode = "list", length = nrow(model_setup))
    #Open all files and merge into one
    for (r in model_setup$model_id){
      #Proper df name
      flux_name <- paste(out_path, paste(r, flux_base, sep = ""), sep = "/")
      #Load csv and append to flux_df
      flux_df[[r]] <- as.data.frame(read.csv(flux_name))
    }

    flux_df <- dplyr::bind_rows(flux_df) %>%
      dplyr::select(-title) %>%
      dplyr::mutate(model_id = rep(model_setup$model_id, each = 5))

    #Save this
    write.csv(flux_df, paste(out_path,
                    paste("total_stocks_comparison_", Sys.Date(), "_.csv", sep = ""), sep ="/"))
  }

  #Gather all statistical model evaluations
  stat_base <- paste("_statistics_", Sys.Date(), ".csv", sep = "")
  stat_df <- vector(mode = "list", length = nrow(model_setup))
  #Open all files and merge into one
  for (r in model_setup$model_id){
   #Proper df name
   stat_name <- paste(out_path, paste(r, stat_base, sep = ""), sep = "/")
   #Load csv and append to stat_df
   stat_df[[r]] <- as.data.frame(read.csv(stat_name))
  }

  #Bind list of dfs into df
  stat_df <- dplyr::bind_rows(stat_df)

  #Save it
  write.csv(stat_df, paste(out_path,
                           paste("statistics_comparison_", Sys.Date(), "_.csv", sep = ""), sep ="/"))


  #Remove all elements in H2o server and shut pipeline down
  h2o::h2o.removeAll()
  #h2o::h2o.shutdown(prompt = FALSE)

}