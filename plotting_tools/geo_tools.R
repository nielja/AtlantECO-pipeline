#Tools to compute distance between grid lines and calculate grid sizes in m2
library(plyr)
library(dplyr)


deg2rad <- function(deg){
  #' Degree to radians
  #'
  #' @param deg     numeric value in degrees
  #'
  #' @returns       numeric value in radians

   (deg * pi) / (180)
}


haversine_lat <- function(lat){
  #' Haversine distance between longitudinal lines
  #'
  #' This function calculates the haversine distance between two lines of longitude at a given latitude.
  #'
  #' @param lat     Numeric value of degree latitudes
  #'
  #' @returns       Numeric value of the distance
  #' @export

  # convert decimal degrees to radians
  lat1 <- deg2rad(lat)
  lat2 <- lat1

  # haversine formula
  dlon = deg2rad(1)
  dlat = lat2 - lat1
  a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
  c = 2 * asin(sqrt(a))
  r = 6371 * 1000 # Radius of earth in m
  return(c * r)}

area_lat <- function(lat){
  #' Calculate the area in m2 of a given grid cell
  #'
  #' This function calculates the area of a 1 x 1 degree grid cell at a given latitude as a basic
  #' trapeze.
  #'
  #' @param lat     Numeric value of degree latitudes
  #'
  #' @returns       Numeric value of the area in m2
  #' @export

  #Basic distance between two latitudes in m (this is the same
  #as the distnace between two longitudes at the equator, i.e. ~111km)
  dlat_base <- haversine_lat(0)

  return(dlat_base * haversine_lat(lat))
}


aggregate_function <- function(df, aggregation_level,
                               mean_cols, latitudinal  = FALSE){
  #' Spatial aggregation function
  #'
  #' This function calculates means of specified columns over a certain specified grid cell area.
  #'
  #' @param df                  data frame on which to conduct calculations
  #' @param aggregation_level   Int number of degrees that the grid cell should have at each side
  #' @param mean_cols           Vector of column names to calculate the mean of
  #' @param latitudinal         Bool, if True, calculate latitudinal averages, otherwise cell-wise
  #'
  #' @returns                   data frame with spatially averaged data in the mean_cols
  #' @export

  #Check for proper coordinate names
  if ("decimalLongitude" %in% colnames(df)){
    lon_name <- "decimalLongitude"
    lat_name <- "decimalLatitude"
  } else if ("Longitude" %in% colnames(df)){
    lon_name <- "Longitude"
    lat_name <- "Latitude"
  } else {
    lon_name <- "lon_gridded"
    lat_name <- "lat_gridded"
  }

  df_out <- df %>%
    dplyr::ungroup() %>%
    dplyr::mutate(lon_gridded = plyr::round_any(get(lon_name), aggregation_level),
           lat_gridded = plyr::round_any(get(lat_name), aggregation_level)) %>%
    dplyr::group_by(lon_gridded, lat_gridded, Month) %>%
    dplyr::summarise_at(mean_cols, mean, na.rm = TRUE) %>%
    dplyr::ungroup()

  #If we want latitudinal means, group just by lat_gridded and Month
  if (latitudinal == TRUE){
    df_out <- df %>%
      dplyr::mutate(lon_gridded = plyr::round_any(get(lon_name), aggregation_level),
             lat_gridded = plyr::round_any(get(lat_name), aggregation_level)) %>%
      dplyr::group_by(lat_gridded, Month) %>%
      dplyr::summarise_at(mean_cols, mean, na.rm = TRUE) %>%
      dplyr::ungroup()
  }

  return(df_out)}