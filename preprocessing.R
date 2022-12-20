library(plyr)
library(dplyr)
library(lubridate)
library(tidyverse)
library(here)
source(here::here("./General_Pipeline_PCA/Code/plotting_tools/load_env_data.R"))

preprocessing <- function(df_name, df_save_name, depth_ave = 200,
                          target_col = "MeasurementValue",
                          env_path = here::here("./General_Pipeline_PCA/Data/environmental_climatologies/"),
                          env_vars = c("BBP443", "Chl_a", "Dist2Coast", "EKE", "Kd_490",
                                       "MLD_MLD", "PAR", "PIC", "Wind", "Silicate_MLD", "Silicate_MLD_gradient",
                                       "Silicate_surface","Nitrate_MLD", "Nitrate_MLD_gradient", "Nitrate_surface",
                                       "Oxygen_depth_200", "Oxygen_MLD", "Oxygen_surface",
                                       "Phosphate_MLD", "Phosphate_MLD_gradient", "Phosphate_surface",
                                       "Salinity_MLD", "Salinity_surface", "Temperature_MLD",
                                       "Temperature_MLD_gradient", "Temperature_surface", "z_eu",
                                       "log_Nitrate_surface", "log_Nitrate_MLD", "log_Phosphate_surface",
                                       "log_Phosphate_MLD", "log_Silicate_surface", "log_Silicate_MLD",
                                       "log_MLD", "log_Chl_a", "log_PIC", "log_BBP443", "log_EKE",
                                       "log_Dist2Coast", "log_Kd_490","tAlk", "DIC", "pCO2", "RevelleFactor",
                                       "Omega_Ca", "Omega_Ar", "log_tAlk","P_star_surface",
                                       "P_star_MLD", "Si_star_surface", "Si_star_MLD", "log_Si_star_surface",
                                       "log_Si_star_MLD", "log_P_star_surface", "log_P_star_MLD")){

  #' Function to conduct preprocessing on a plankton data set.
  #'
  #' Preprocessing involves loading the data frame, checking whether all units are the same,
  #' dropping all values with a NaN in  MeasurementValues, eventDate, Depth, decimalLatitude and decimalLongitude
  #' Also we conduct re-gridding to WOA18, summing and averaging over the specified depth depth_ave
  #' Last, we connect environmental observations to the gridded dataset and drop coastal values with salinity_surface < 30.
  #' The dataset is returned and also saved to Data/plankton_data/plankton_data_surface_aggregated_w_env_data.csv

  #' @param df_name      path to plankton data frame in AtlantECO format
  #' @param df_save_name character name to add to the save name
  #' @param depth_ave    numeric value, depth of the surface in meter over which we want to aggregate, defaults to 200m,
  #'                     if this is NA, depth averaging has already been conducted
  #' @param target_col   column name of the column we want to use as target variable in the modelling
  #' @param env_path     character string of path to environmental datasets
  #' @param env_vars     vector of characters to specify which environmental variables to load, defaults to all
  #'
  #' @returns            returns a preprocessed dataset with the columns the original and log10(x + 1) transformed target variable

  #Load the data set
  if (grepl(".csv", df_name, fixed = TRUE)){
    df <- read.csv(paste(here::here("General_Pipeline_PCA/Data/plankton_data"), df_name, sep = "/"))
  } else if (grepl(".RData", df_name, fixed = TRUE)){
    loadRData <- function(fileName){
      #' Loads an RData file, and returns it
      #'
      #' @param fileName    path to the file
      #' @returns           data frame under that path
      #'

      load(fileName)
      get(ls()[ls() != "fileName"])
    }
    df <- loadRData(paste(here::here("General_Pipeline_PCA/Data/plankton_data"), df_name, sep = "/"))
  } else if (grepl(".txt", df_name, fixed = TRUE)){
    df <- read.table(paste(here::here("General_Pipeline_PCA/Data/plankton_data"), df_name, sep = "/"))
  }

  #Check whether all units are correct
  if (!is.na(depth_ave) & (length(unique(df$MeasurementUnit)) > 1)){
      print("Potential problem, there is more than one unique unit.")
      print(unique(df$MeasurementUnit))
  } else {

    #Get original number of observaions
    orig_length <- nrow(df)

    if (!(is.na(depth_ave))){
      #Fill in any missing Depth values where we do have max and min depth
      df[is.na(df$Depth),"Depth"] = rowMeans(df[is.na(df$Depth),c("MinDepth", "MaxDepth")])
      #Check columns
      c_cols <- c("eventDate", "Depth", "decimalLatitude", "decimalLongitude", target_col)
    } else {c_cols <- c("decimalLatitude", "decimalLongitude", target_col)}

    #Drop all observations with NaNs in the relevant columns
    df <- df %>%
      dplyr::mutate(eventDate = lubridate::make_date(Year, Month, Day)) %>%
      tidyr::drop_na(any_of(c_cols))
    print(paste("Dropped", orig_length- nrow(df), "data points because of quality issues, i.e. NaNs where they shouldn't be."))

    #Conduct regridding to WOA18
    df <- df %>%
      dplyr::mutate(lat_gridded = (plyr::round_any(decimalLatitude * 4 + 2, 4) - 2)/4,
             lon_gridded = (plyr::round_any(decimalLongitude * 4 + 2, 4) - 2)/4)

    #If we don't have a depth averaging yet, conduct it
    if (!(is.na(depth_ave))){
      #Conduct surface value aggregation for the entire data set
      df_agg <- df %>%
        #Drop all values from outside the defined depth range
        dplyr::filter(Depth <= depth_ave) %>%
        #Add distinct samples from the same tow
        dplyr::group_by(decimalLatitude, decimalLongitude, eventDate,
                 orig_occurrenceID, Depth, Year, Month, Day, lon_gridded, lat_gridded) %>%
        dplyr::summarise(MeasurementValue = sum(get(target_col))) %>%
        #Average over distinct sampling events at each gridded point and month
        dplyr::ungroup() %>%
        dplyr::group_by(lat_gridded, lon_gridded, Month) %>%
        dplyr::summarise(MeasurementValue = mean(MeasurementValue))

      print(paste("Surface aggregated dataset now has", nrow(df_agg), "data points."))
    } else {
      df_agg <- df %>%
            dplyr::mutate(target_intermediate = get(target_col)) %>%
            dplyr::select(lat_gridded, lon_gridded, Month, target_intermediate) %>%
            dplyr::mutate(MeasurementValue = target_intermediate) %>%
            dplyr::select(-target_intermediate)
    }

    #Load environmental climatologies for the desired variables
    df_env <- load_env_data(predictor_list = env_vars,
                            env_path = env_path)
    df_env <- df_env %>%
      dplyr::mutate(lon_gridded = Longitude,
             lat_gridded = Latitude)

    #Collocate environmental data to plankton observations
    df_col <- df_agg %>%
      dplyr::left_join(., df_env, by = c("lon_gridded", "lat_gridded", "Month"))

    #Drop coastal values with surface salinity < 30
    df_col <- df_col %>%
      dplyr::filter(Salinity_surface >= 30)

    #Add log10 abundance
    df_col <- df_col %>%
      dplyr::mutate(log10_MV = log10(MeasurementValue + 1))
    colnames(df_col)[colnames(df_col) == "MeasurementValue"] = target_col

    #Save this collocated dataframe and also return it
    save_path <- paste(here::here("./General_Pipeline_PCA/Data/plankton_data/"),
                       paste0(df_save_name, "_data_surface_aggregated_w_env_data.csv"), sep = "/")
    print(paste0("saving the preprocessed dataset to ", save_path))
    write.csv(df_col, save_path)

    return(df_col)
  }
}