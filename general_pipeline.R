library(here)
#library(devtools)
source(here::here("./General_Pipeline_PCA/Code/preprocessing.R"))
#source(here::here("./General_Pipeline_PCA/Code/plotting_tools/map_plot_basic.R"))
source(here::here("./General_Pipeline_PCA/Code/univariate_modelling/varimp_plot_with_dendrogram.R"))
source(here::here("./General_Pipeline_PCA/Code/univariate_modelling/univariate_analyses_h2o.R"))
source(here::here("./General_Pipeline_PCA/Code/multivariate_modelling/h2o_model_pipeline.R"))

#Install the package and take the functions from there
#setwd("/net/meso/work/nknecht/Masterarbeit")
#install("atlantecoSDM")
### === Pipeline to analyse and extrapolate a data set in AtlantECO format ====
#'Add the name of the dataset deposited in "Data/plankton_data". This has to be either .RData or .csv
df_name <- "AtlantECO-BASEv1_dataset_Jellyfish_abundances_22_04_22.RData"
#prefix to add to the preprocessed saved dataset
save_name <- "jellyfish"
#output path -> make sure this exists and create the folders UnivariateModels and Models within
out_path <- here::here("./General_Pipeline_PCA/Model_Output/Jellyfish_PCA")
#set the unit that your measurements are in -> for a nicer layout, use bquote() but you can also just put a characters string
unit_measurement <- bquote("ind."~m^{-3})

### ==== Pre-processing of the dataset ====
#' Pre-processing: drop all NaN MeasurementValues, eventDate, Depth, decimalLatitude and decimalLongitude
#' Re-gridding to WOA18, summing and averaging over the specified depth depth_ave
#' We also collocate all environmental variables that we have collected and compute a log10(x + 1)
#' transformation of the target_col, this will be stored as log10_MV and is the default
#' target for the following analyses. You can change it under the definition of target in the next
#' section.

df <- preprocessing(df_name, #Dataframe path
                    df_save_name = save_name, #Prefix to add to the saved data set
                    depth_ave = 200, #Depth in m over which we aggregate
                    target_col = "MeasurementValue") #Target column that should be preserved


#### ==== Run univariate analysis ====
#' Univariate analysis to find out which environmental predictors are most
#' important for the given target variable. The results of the analysis will
#' be saved in the folder UnivariateModels. Choose a subset of predictors to
#' analyse in the line below. Also set the target variable of interest.

#Choose predictor variables we want to work with and do univariate analysis
predictors <- c("Temperature_MLD", "log_Chl_a", "log_Nitrate_MLD", "log_Phosphate_MLD",
                "Salinity_MLD","Oxygen_depth_200", "log_MLD", "log_EKE")
target <- "log10_MV"

#Different combination of spatial aggregation levels (10°, 5°, 1°) and latitudinal or not
aggregation_level = c(10, 5, 1)
latitudinal = c(TRUE, FALSE)
option_matrix <- expand.grid(aggregation_level, latitudinal)
colnames(option_matrix) <- c("aggregation_level", "latitudinal")

#Apply the univariate analysis function function to the grid of options defined above
mapply(univariate_analysis_env_space,
       df = list(df),
       aggregation_level = option_matrix$aggregation_level,
       predictors = list(predictors),
       target = target,
       #Output path where the results of the univariate analysis will be stored
       out_path = paste0(out_path, "/UnivariateModels/"),
       pdp_calculation = TRUE,
       env_background = FALSE,
       latitudinal = option_matrix$latitudinal,
       remove_outliers = TRUE)


#Plot the actual deviance explained along with a dendrogram
#' Note: the red dashed line indicates a correlation of |r| = 0.7. Predictors in the final
#' model should not be correlated more highly than this.
varimp_plot_univariate_dendro(df = df,
                              df_path = paste0(out_path, "/UnivariateModels"),
                              predictors = predictors,
                              d_nam = Sys.Date())

### === Multivariate modelling ========
#' Create an ensemble of five multivariate models to predict the target variable as
#' a function of environmental predictors. Based on the dendrogram and the percentage
#' of deviance each environmental predictor explains, choose the final predictor dataset.
#' The results of the multivariate modelling will be saved in the /Models/no_tuning folder,
#' but you can also adapt this. Default is no model tuning as this is faster and the results
#' generally don't differ significantly.

#Select the final predictors you want to use and let's run the multivariate models
predictor_choices <- list(c("Temperature_MLD", "log_Chl_a", "Oxygen_depth_200", "log_MLD"))

#If wanted, specify a path to an environmental field that includes the same predictors as used in the model
#external_field_path <- here::here("./General_Pipeline_PCA/Data/external_env_climatologies/fields_for_ml_corrected.csv")
external_field_path <- ""

#Modify this if wanted
pca_choice = TRUE #whether to run a PCA on the environmental predictors before model training

### ==== Don't modify =======================
### Preparation for the multivariate model setup, don't modify this
target <- c("log10_MV")
target_options <- length(target)
base_path <- "General_Pipeline_PCA/Data/plankton_data/"

### Setup of the model basics
dataset_choices <- data.frame(#Path to the surface aggregated plankton data frame with collocated environmental data
                              data_path = paste0(save_name, "_data_surface_aggregated_w_env_data.csv"),
                              #Path to store the results
                              out_path = paste0(out_path, "/Models"),
                              #Name given to the dataframe, not relevant for analysis, just for saving
                              df_name = c("all_w_zeros"),
                              #Region of analysis: full, SO (SouthernOcean), NAP (North Atlantic and Pacific)
                              region = c("full"),
                              #Name of a column that we want to use to weight the observations -> unless there is
                              #a very good reason to do this, don't use it, i.e. set weights_column = FALSE and
                              #weights_options = "none"
                              weights_column = FALSE,
                              weight_options = "none") %>%
  mutate(data_path = paste(base_path, data_path, sep = ""))

### === Run the model training and model assessment ====
model_pipeline(dataset_choices,
               predictor_choices,
               target,
               tune_model = FALSE, #Change this to TRUE if you want to conduct model tuning -> but it takes very long!
               seasonal = TRUE, #Seasonally averaged output plots as addition to the yearly averages
               monthly = FALSE, #Monthly output plots -> only useful to make a gif
               external_field_path = external_field_path, #path to external nvironmental predictor field
               pca = pca_choice, #whether to run a PCA on the environmental predictors before modelling
               remove_highly_correlated = TRUE, #don't set this to FALSE unless you have a good reason
               unit_small = unit_measurement, #unit that the biomass values are in
               operative_system = "linux")

### === Alternative: run the model in evaluation mode only once we have already trained it ====
#' NOTE: the model_paths has to be adapted to where the models are saved, this only works after you have run
#' the model training once
model_paths <- "./General_Pipeline_PCA/Salps_wo_CPR/Models/1_2022-06-02/"

#Don't modify this part
model_names <- c("GLM_Model", "GAM_Model", "RF_Model", "GBM_Model", "DL_Model")
model_paths <- paste0(model_paths, model_names)
model_pipeline(dataset_choices,
               predictor_choices,
               target,
               model_paths = model_paths,
               evaluation_mode = TRUE, #do not train the models, just evaluate
               tune_model = FALSE, #Change this to TRUE if you want to conduct model tuning -> but it takes very long!
               seasonal = TRUE, #Seasonally averaged output plots as addition to the yearly averages
               monthly = FALSE, #Monthly output plots -> only useful to make a gif
               external_field_path = external_field_path, #path to external nvironmental predictor field
               unit_small = unit_measurement)


