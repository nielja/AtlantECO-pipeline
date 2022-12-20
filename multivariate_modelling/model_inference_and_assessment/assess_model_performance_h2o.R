#Script to load the trained H2O models and assess their performance and their predictions
library(h2o)
library(here)
#Set source for other code
source(here::here("./General_Pipeline_PCA/Code/multivariate_modelling/model_inference_and_assessment/prediction_map_plot.R"))
source(here::here("./General_Pipeline_PCA/Code/multivariate_modelling/model_inference_and_assessment/varimp_heatmap.R"))
source(here::here("./General_Pipeline_PCA/Code/multivariate_modelling/model_inference_and_assessment/pdp_plot_h2o.R"))
source(here::here("./General_Pipeline_PCA/Code/multivariate_modelling/model_inference_and_assessment/model_statistics.R"))
source(here::here("./General_Pipeline_PCA/Code/multivariate_modelling/model_inference_and_assessment/residual_analysis.R"))

assess_model_performance_h2o <- function(model_paths,
                                         model_name_list,
                                         model_id,
                                         predictors,
                                         target,
                                         h2o_train,
                                         h2o_test,
                                         out_path,
                                         df_choice,
                                         remove_outliers,
                                         weight_option,
                                         region,
                                         predictor_name_list,
                                         data_path,
                                         #Variable importance heatmap
                                         varimp_heatmap_func = TRUE,
                                         varimp_permutation = TRUE,
                                         #PDP plot
                                         pdp_plot_func = TRUE,
                                         pdp_point_density = FALSE,
                                         pdp_rug = TRUE,
                                         pdp_scatter = FALSE,
                                         pdp_return_df = FALSE,
                                         #Model statistics as csv
                                         model_statistics_func = TRUE,
                                         #Prediction map
                                         prediction_map_func = TRUE,
                                         env_path = here::here("./Data/21_10_18_environmental_data/all_env_vars_all_months/"),
                                         external_field_path = "",
                                         pca = FALSE,
                                         res.pca = NA,
                                         model_predictors = c(),
                                         prediction_monthly = FALSE,
                                         prediction_seasonally = TRUE,
                                         prediction_MESS = TRUE,
                                         prediction_mess_plot = TRUE,
                                         prediction_std_plot = TRUE,
                                         prediction_save_plot = TRUE,
                                         prediction_region = "full",
                                         prediction_min_col = "",
                                         prediction_max_col = "",
                                         prediction_add_points = FALSE,
                                         unit_small = bquote('#'),
                                         unit_large = bquote('#'),
                                         unit_conversion = 1,
                                          #Residual analysis
                                         residual_analysis_func = TRUE){
    #' Wrapper to run all model assessment and inference functions
    #'
    #' This function runs the variable importance analysis, PDP curve calculation, model statistics calculation,
    #' spatial extrapolation of results and residual analysis. The single analyses can be switched on and off as wanted.
    #'
    #' @param model_paths             list of paths to where the models are stored
    #' @param model_name_list         list of model names
    #' @param model_id                Integer model id, only relevant in case of multiple models trained at the same time
    #' @param predictors              vector of names of predictors to use
    #' @param target                  vector of name of target variable
    #' @param h2o_train               H2O training dataframe
    #' @param h2o_test                H2O testing dataframe
    #' @param out_path                path to folder where the plots should be saved
    #' @param df_choice               character name of the dataframe
    #' @param remove_outliers         Bool, whether outliers were removed in the modelling
    #' @param weight_option           character name of the column where the weights are stored in
    #' @param region                  Character, "full", "SO", "NA", "NP" depending on whether we want a focus on one of the CPR regions
    #' @param predictor_name_list     vector of bquote predictor names for the PDP plots
    #' @param data_path               path to the plankton observation dataset, from here() directory down
    #' @param varimp_heatmap_func     Bool, if True, run the variable importance analysis
    #' @param varimp_permutation      Bool, if True, base the variable importance calculation on a permutation analysis
    #' @param pdp_plot_func           Bool, if True, run the partial depence plot analysis
    #' @param pdp_point_density       Bool, if True, show a point density plot for the PDP
    #' @param pdp_rug                 Bool, if True, show a rug of observation values in the PDP
    #' @param pdp_scatter             Bool, if True, show scatter points of the observation data in the PDP
    #' @param pdp_return_df           Bool, if True, return the pdp dataframe as output of the function -> only for testing purposes
    #' @param model_statistics_func   Bool, if True, run the model statistics analysis
    #' @param prediction_map_func     Bool, if True, run the extrapolation of the target variable to total region
    #' @param env_path                path to the folder with the environmental climatologies
    #' @param external_field_path     Character string of path to external field data set, if "", use WOA fields
    #' @param pca                     Bool, if True, compute PCA-transformation of the environmental predictors and train the models on these
    #' @param res.pca                 PCA results element to compute transformations
    #' @param model_predictors        Vector of original predictor names, only applicable if we have pca == TRUE
      #' @param prediction_monthly      Bool, if True, produce monthly plots of the prediction value
      #' @param prediction_seasonally   Bool, if True, produce seasonal plots of the prediction value
      #' @param prediction_MESS         Bool, if True, run a MESS analysis and include stippling in the prediction plots
      #' @param prediction_mess_plot    Bool, if True, produce a map showing the number of months with MESS < 0
      #' @param prediction_std_plot     Bool, if True, produce a map plot for the standard deviation between the model predictions
      #' @param prediction_save_plot    Bool, if True, save the prediction plots
      #' @param prediction_region       Character, "full", "SO", "NA", "NP" , region to conduct the extrapolation for
      #' @param prediction_min_col      Numeric, minimum value for colorscale of prediction plot, if "", will be set automatically
      #' @param prediction_max_col      Numeric, maximum value for colorscale of prediction plot, if "", will be set automatically
      #' @param prediction_add_points   Bool, if True, add observation points to the prediction maps
      #' @param unit_small              bquote version of the unit the original target variable is in
      #' @param unit_large              bquote version of an aggregate version of the small unit to mention for the annual aggregate
      #' @param unit_conversion         numerical value giving the conversion factor between unit_small and unit_large
      #' @param residual_analysis_func  Bool, if True, run the analysis of the residual patterns
      #'
      #' @returns                       Saves a number of plots depending on the chosen subfunctions
      #'
      #' @export


  h2o::h2o.no_progress()

  #Create empty list to load models into
  model_list <- vector(mode = "list", length = length(model_name_list))
  names(model_list) <- model_name_list
  #Load the models
  for (i in seq_along(model_name_list)){
    m <- model_name_list[i]
    model_list[[m]] <- h2o::h2o.loadModel(model_paths[[i]])
  }

  #Check if we need to adapt the output folder
  if (external_field_path != ""){
    env_df <- read.csv(external_field_path)
    #Check if we can use this
    if (all(c(predictors, "Longitude", "Latitude", "Month") %in% colnames(env_df))){
      out_path <- paste0(out_path, "/", stringr::str_extract(basename(external_field_path), "[^.]+"))
      print(paste("using out_path", out_path))
      if (dir.exists(out_path) == FALSE){
        dir.create(out_path)
      }
  }}

  ### ======== Plot variable importance ==================
  if (varimp_heatmap_func){
    print("varimp plot")
    varimp_heatmap(model_list = model_list,
                   model_name_list = model_name_list,
                   model_id = model_id,
                   predictors = predictors,
                   h2o_test = h2o_test,
                   permutation = varimp_permutation,
                    out_path = out_path)
  }


  #### ====== MAKE PARTIAL DEPENDENCE PLOTSSSSSS ===========
  if (pdp_plot_func){
    print("pdp plot")
    envir = environment()
    pdp_path <- pdp_plot_multi(model_list = model_list,
                   model_name_list = model_name_list,
                   model_id = model_id,
                   predictors = predictors,
                   target = target,
                   h2o_train = h2o_train,
                   point_density = pdp_point_density,
                   rug = pdp_rug,
                   scatter = pdp_scatter,
                   out_path = out_path,
                  return_df = pdp_return_df,
                  predictor_name_list = predictor_name_list)
    pdp_path = pdp_path[[1]]
  } else {pdp_path = NA}


  ### ==== Prediction plots ==========
  if (prediction_map_func){
    print("prediction map function")
    ave_nSD <- prediction_map_plot(predictors = predictors,
                        target = target,
                        model_list = model_list,
                        model_name_list = model_name_list,
                        model_id = model_id,
                        h2o_train = h2o_train,
                        data_path = data_path,
                        pca = pca,
                        res.pca = res.pca,
                        model_predictors = model_predictors,
                        external_field_path = external_field_path,
                        unit_small = unit_small,
                        unit_large = unit_large,
                        unit_conversion =  unit_conversion,
                        df_choice = df_choice,
                        remove_outliers = remove_outliers,
                        weight_option = weight_option,
                        env_path = here::here("./Data/21_10_18_environmental_data/all_env_vars_all_months/"),
                        monthly = prediction_monthly,
                        seasonal = prediction_seasonally,
                        MESS = prediction_MESS,
                         region = prediction_region,
                        mess_plot = prediction_mess_plot,
                        std_plot = prediction_std_plot,
                        save_pred_plot = prediction_save_plot,
                        min_col = prediction_min_col,
                        max_col = prediction_max_col,
                        out = out_path,
    add_points = prediction_add_points)
  } else {ave_nSD = NaN}

  #Performance output as data.frame
  if (model_statistics_func){
    print("model stats")
    model_stats_list(model_list,
                     h2o_test = h2o_test,
                     h2o_train = h2o_train,
                     ave_nSD = ave_nSD,
                     pdp_path = pdp_path,
                     model_names = model_name_list,
                     model_id = model_id,
                     out_path = out_path,
                     predictors = predictors,
                     target = target,
                     df_choice = df_choice,
                     remove_outliers = remove_outliers,
                     weight_option = weight_option,
                     region = region)
  }

  ### ==== Residual analysis ==============
  if (residual_analysis_func){
    residual_analysis(h2o_train = h2o_train,
                      h2o_test = h2o_test,
                      predictors = predictors,
                      target = target,
                      model_list = model_list,
                      model_id = model_id,
                      model_name_list = model_name_list,
                      region = prediction_region,
                      out = out_path,
                      unit_small = unit_small)
  }

}