library(ncdf4)
library(dplyr)
library(tidyr)
library(here)

load_env_data <- function(predictor_list,
                          decimalCoords = FALSE,
                          env_path = here::here("./Data/21_10_18_environmental_data/all_env_vars_all_months/")){
    #' Load environmental data files
    #'
    #' This function loads environmental data for all predictors in predictor_list and return a dataframe with them
    #'
    #' @param predictor_list    vector with names of predictors to be loaded
    #' @param decimalCoords     Bool, if True, columns for Longitude and Latitude will be named decimalLongitudet etc.
    #' @param env_path          path to where the environmental data files are
    #'
    #' @return                  dtaframe with all environmental variables in the wide form

  #Lon and Lat names
  if (decimalCoords == TRUE){
    lon_name = "decimalLongitude"
    lat_name = "decimalLatitude"
  } else {
    lon_name = "Longitude"
    lat_name = "Latitude"
  }

  #Load environmental variables
  setwd(env_path)
  #Map variable names (i.e. column names in the dataframe to env var dataframe names)
  map_var_paths = data.frame(VarName = c("BBP443", "Chl_a", "Dist2Coast", "EKE", "Kd_490",
                                         "MLD_MLD", "PAR", "pco2_vars", "PIC", "Wind", "Silicate_MLD", "Silicate_MLD_gradient",
                                         "Silicate_surface","Nitrate_MLD", "Nitrate_MLD_gradient", "Nitrate_surface",
                                         "Oxygen_depth_200", "Oxygen_MLD", "Oxygen_surface",
                                         "Phosphate_MLD", "Phosphate_MLD_gradient", "Phosphate_surface",
                                         "Salinity_MLD", "Salinity_surface", "Temperature_MLD",
                                         "Temperature_MLD_gradient", "Temperature_surface", "z_eu"),
                             Path = c('bbp_443_gsm_SeaWIFS_allmonths.nc',
                                      'chlor_a_SeaWIFS_allmonths.nc',
                                      'dist2coast_allmonths.nc',
                                      'EKE_allmonths.nc',
                                      'Kd_490_SeaWIFS_allmonths.nc',
                                      'MLD_SODA.nc',
                                      'par_SeaWIFS_allmonths.nc',
                                      'pco2_related_vars.nc',
                                      'pic_SeaWIFS_allmonths.nc',
                                      'wind_allmonths.nc',
                                      'woa18_all_i_allmonths_mld_average_SODA3.4.2.nc',
                                      'woa18_all_i_allmonths_mld_gradient.nc',
                                      'woa18_all_i_allmonths_surface.nc',
                                      'woa18_all_n_allmonths_mld_average_SODA3.4.2.nc',
                                      'woa18_all_n_allmonths_mld_gradient.nc',
                                      'woa18_all_n_allmonths_surface.nc',
                                      'woa18_all_o_allmonths_depth_200.nc',
                                      'woa18_all_o_allmonths_mld_average_SODA3.4.2.nc',
                                      'woa18_all_o_allmonths_surface.nc',
                                      'woa18_all_p_allmonths_mld_average_SODA3.4.2.nc',
                                      'woa18_all_p_allmonths_mld_gradient.nc',
                                      'woa18_all_p_allmonths_surface.nc',
                                      'woa18_decav_s_allmonths_mld_average_SODA3.4.2.nc',
                                      'woa18_decav_s_allmonths_surfaceSSS.nc',
                                      'woa18_decav_t_allmonths_mld_average_SODA3.4.2.nc',
                                      'woa18_decav_t_allmonths_mld_gradient.nc',
                                      'woa18_decav_t_allmonths_surfaceSST.nc',
                                      'Zeu_lee_SeaWIFS_allmonths.nc'))

  map_log_paths = data.frame(VarName = c("log_Nitrate_surface",
                                         "log_Nitrate_MLD",
                                         "log_Phosphate_surface",
                                         "log_Phosphate_MLD",
                                         "log_Silicate_surface",
                                         "log_Silicate_MLD",
                                         "log_MLD",
                                         "log_Chl_a",
                                         #"log_tAlk",
                                         "log_PIC",
                                          "log_BBP443",
                                         "log_EKE",
                                          "log_Dist2Coast",
                                         "log_Kd_490"))
  map_log_paths["VarNameAdapted"] = substr(map_log_paths$VarName, 5, 25)
  map_log_paths[map_log_paths$VarName == "log_MLD", "VarNameAdapted"] = "MLD_MLD"
  #map_log_paths[map_log_paths$VarName == "log_tAlk", "VarNameAdapted"] = "talk"

  #Star paths
  map_star_paths = data.frame(VarName = c("P_star_surface",
                                          "P_star_MLD",
                                          "Si_star_surface",
                                          "Si_star_MLD",
                                          "log_Si_star_surface",
                                          "log_Si_star_MLD",
                                          "log_P_star_surface",
                                          "log_P_star_MLD"))


  #Make dataframe of the pCO2 variables and their equivalent name in the nc-files
  pco2_var_df <- data.frame(orig = c("tAlk", "DIC", "pCO2", "RevelleFactor",
                                     "Omega_Ca", "Omega_Ar", "log_tAlk"),
                            nc_name = c("talk", "dic", "spco2", "revelle_factor",
                                        "omega_ca", "omega_ar", "talk"))

  #Loop over all desired predictors
  star_vars <- c()
  for (nv in seq_along(predictor_list)){
    var <- predictor_list[nv]
    log_val <- FALSE
    #If we want to plot pco2-related vars
    if (var %in% pco2_var_df$orig){
      var_choice <- "pco2_vars"
      if (var == "log_tAlk"){
        log_val <- TRUE
      }
      #Map the nc_name to var
      var <- pco2_var_df[pco2_var_df$orig == var, "nc_name"]


    } else {var_choice <- var}

    #If we have a log transformed variable
    if (var %in% map_log_paths$VarName){
      var_choice <- map_log_paths[map_log_paths$VarName == var, "VarNameAdapted"]
      log_val <- TRUE
    }

    #If we have a star variable
    if (var %in% map_star_paths$VarName){
      #Skip this variable and calculate it later
      star_vars <- append(star_vars, var)
      var_choice <- NaN
    }

    if (!is.na(var_choice)){
      #Load necessary netcdf files
      var_nc <- ncdf4::nc_open(map_var_paths[map_var_paths$VarName == var_choice, "Path"])
      #Get varname and load the var values
      var_nam <- names(var_nc$var)

      if (length(var_nam)>1){
        var_nam <- var
      }
      var_val <- ncvar_get(var_nc, var_nam)
      if (log_val == T){
        var_val <- log10(var_val)
      }

      #Check dimensions so that they are all in the same order (month, lat, lon)
      if (sum((dim(var_val) == c(360, 180, 12))) == 3){
        var_val = var_val
        env_sub_df = data.frame(c(var_val))
        colnames(env_sub_df)[1] <- var
        #Add lon and lat columns
        env_sub_df[lon_name] = rep(-179.5:179.5)
        env_sub_df[lat_name] = rep(-89.5:89.5, each = 360)
        #Add a month column to env_df
        env_sub_df["Month"] = rep(1:12, each = 360*180)
        #Sort it
        env_sub_df <- env_sub_df[
          with(env_sub_df, order(get(lon_name), get(lat_name), Month)),
        ]

      } else if (sum((dim(var_val) == c(12, 360, 180))) == 3){
        var_val = aperm(var_val, c(2, 3, 1))
        env_sub_df = data.frame(c(var_val))
        colnames(env_sub_df)[1] <- var
        #Add lon and lat columns
        env_sub_df[lon_name] = c(var_nc$dim$lon$vals) #rep(-179.5:179.5)
        env_sub_df[lat_name] = rep(var_nc$dim$lat$vals, each = length(var_nc$dim$lon$vals)) #rep(89.5:-89.5, each = 360)
        #Add a month column to env_df
        env_sub_df["Month"] = rep(1:12, each = 360*180)
        #Sort it
        env_sub_df <- env_sub_df[
          with(env_sub_df, order(get(lon_name), get(lat_name), Month)),
        ]
      }
      if (nv == 1){
        env_df <- env_sub_df
      } else {
        if ((var == "talk") & (log_val == TRUE)){
          var_name = "log_tAlk"
        } else (var_name = var)
        env_df[var_name] <- env_sub_df[var]}
    }
  }

  #Rename the columns that have something to do with pco2 variables
  pco2_var_df <- pco2_var_df[pco2_var_df$orig != "log_tAlk",]
  for (p_nam in pco2_var_df$nc_name){
    colnames(env_df)[colnames(env_df) == p_nam] <- pco2_var_df[pco2_var_df$nc_name == p_nam, "orig"]
  }

  #If we have star-variables
  if (length(star_vars) > 0){
    if ("P_star_surface" %in% star_vars){
      env_df["P_star_surface"] <- env_df["Phosphate_surface"] - env_df["Nitrate_surface"]/16
    }
    if ("P_star_MLD" %in% star_vars){
      env_df["P_star_MLD"] <- env_df["Phosphate_MLD"] - env_df["Nitrate_MLD"]/16
    }
    if ("Si_star_surface" %in% star_vars){
      env_df["Si_star_surface"] <- env_df["Silicate_surface"]/env_df["Nitrate_surface"]
    }
    if ("Si_star_MLD" %in% star_vars){
      env_df["Si_star_MLD"] <- env_df["Silicate_MLD"]/env_df["Nitrate_MLD"]
    }
    if ("log_Si_star_surface" %in% star_vars){
      env_df["log_Si_star_surface"] <- log10(env_df["Silicate_surface"]/env_df["Nitrate_surface"])
    }
    if ("log_Si_star_MLD" %in% star_vars){
      env_df["log_Si_star_MLD"] <- log10(env_df["Silicate_MLD"]/env_df["Nitrate_MLD"])
    }
    if ("log_P_star_surface" %in% star_vars){
      env_df["log_P_star_surface"] <- log10(env_df["Phosphate_surface"] - env_df["Nitrate_surface"]/16)
    }
    if ("log_P_star_MLD" %in% star_vars){
      env_df["log_P_star_MLD"] <- log10(env_df["Phosphate_MLD"] - env_df["Nitrate_MLD"]/16)
    }
  }

  #Replace inf with nan
  env_df[mapply(is.infinite, env_df)] <- NA

  #Round all coordinates properly to .5
  env_df <- env_df %>%
    dplyr::mutate(lon_name_val = round(get(lon_name), 1),
           lat_name_val = round(get(lat_name), 1))
  env_df[lon_name] = env_df$lon_name_val
  env_df[lat_name] = env_df$lat_name_val
  env_df <- env_df %>%
    dplyr::select(-lon_name_val, -lat_name_val)

  #Set working directory back to normal
  setwd(here::here())

  #Return the full dataset
  return(env_df)
}