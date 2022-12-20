library(ggplot2)
library(dismo)
library(ncdf4)
library(viridis)
#library(modEvA)
library(raster)
library(tidyr)
#library(DescTools)
library(dplyr)
library(maps)
library(here)
library(h2o)
library(scales)
library(matrixStats)
#library(parallel)
#library(lubridate)
#library(units)
#library(pracma)
#library(gridExtra)
#library(rgdal)
library(stringr)
source(here::here("./General_Pipeline_PCA/Code/plotting_tools/shift_legend_ggplot.R"))
source(here::here("./General_Pipeline_PCA/Code/plotting_tools/geo_tools.R"))
source(here::here("./General_Pipeline_PCA/Code/plotting_tools/load_env_data.R"))
source(here::here("./General_Pipeline_PCA/Code/plotting_tools/ggplot_tools.R"))


prediction_map_plot <- function(predictors,target,
                                model_list,
                                model_name_list,
                                model_id,
                                h2o_train,
                                data_path,
                                pca = FALSE,
                                res.pca = NA,
                                model_predictors = c(),
                                external_field_path = "",
                                unit_small = bquote('#'),
                                unit_large = bquote('#'),
                                unit_conversion = 1,
                                region = "full",
                                df_choice = "",
                                remove_outliers = "",
                                weight_option = "",
                                env_path = here::here("./General_Pipeline_PCA/Data/environmental_climatologies/"),
                                basin_path = here::here("./General_Pipeline_PCA/Data/basin_data/basin_data_regional.RData"),
                                monthly = FALSE,
                                seasonal = TRUE,
                                MESS = TRUE,
                                mess_plot = FALSE,
                                std_plot = FALSE,
                                save_pred_plot = TRUE,
                                min_col = "",
                                max_col = "",
                                out_path = here::here("./Plots/"),
                                add_points = FALSE){
  #' Create a map of predicted target values
  #'
  #' This function predicts the target variable in the entire specified region based on the environmental conditions.
  #'
  #' @param predictors              vector of names of predictors to use
  #' @param target                  vector of name of target variable
  #' @param model_list              list of the models
  #' @param model_name_list         list of model names for plotting
  #' @param h2o_train               H2O training dataframe
  #' @param data_path               path to the plankton observation dataset, from here() directory down
  #' @param pca                     Bool, if True, compute PCA-transformation of the environmental predictors and train the models on these
  #' @param res.pca                 PCA results element to compute transformations
  #' @param model_predictors        Vector of original predictor names, only applicable if we have pca == TRUE
  #' @param external_field_path     Character string of path to external field data set, if "", use WOA fields
  #' @param unit_small              bquote version of the unit the original target variable is in
  #' @param unit_large              bquote version of an aggregate version of the small unit to mention for the annual aggregate
  #' @param unit_conversion         numerical value giving the conversion factor between unit_small and unit_large
  #' @param region                  Character, "full", "SO", "NA", "NP" , region to conduct the extrapolation for
  #' @param df_choice               character name of the dataframe
  #' @param remove_outliers         Bool, whether outliers were removed in the modelling
  #' @param weight_option           character name of the column where the weights are stored in
  #' @param env_path                path to the folder with the environmental climatologies
  #' @param basin_path              path to the basin definition files
  #' @param monthly                 Bool, if True, produce monthly plots of the prediction value
  #' @param seasonal                Bool, if True, produce seasonal plots of the prediction value
  #' @param MESS                    Bool, if True, run a MESS analysis and include stippling in the prediction plots
  #' @param mess_plot               Bool, if True, produce a map showing the number of months with MESS < 0
  #' @param std_plot                Bool, if True, produce a map plot for the standard deviation between the model predictions
  #' @param save_pred_plot          Bool, if True, save the prediction plots
  #' @param min_col                 Numeric, minimum value for colorscale of prediction plot, if "", will be set automatically
  #' @param max_col                 Numeric, maximum value for colorscale of prediction plot, if "", will be set automatically
  #' @param out_path                Path to save the plots to
  #' @param add_points              Bool, if True, add observation points to the prediction maps
  #'
  #' @returns                       Saves the chosen plots and prediction data frames and returns a data frame of the average standard deviation between model predictions
  #' @export




  #Load environmental predictors
  print("load env predictors")
  if (pca == FALSE){
    if (external_field_path == ""){
      env_df <- load_env_data(predictors, env_path = env_path, decimalCoords = FALSE)
      env_df <- env_df %>%
        dplyr::arrange(Latitude, Longitude)
    } else {#Load external predictor fields
      env_df <- read.csv(external_field_path)
      #Also append out_path then
      if (!all(c(predictors, "Longitude", "Latitude", "Month") %in% colnames(env_df))){
        print("not all columns present in the environmental data field")
        print(colnames(env_df))
        print("using training fields instead")
        env_df <- load_env_data(predictors, env_path = env_path, decimalCoords = FALSE)
      } else {
        env_df <- env_df %>%
          dplyr::select(all_of(c(predictors, "Longitude", "Latitude", "Month"))) %>%
          dplyr::arrange(Latitude, Longitude)
        }
    }
  } else if (pca == TRUE){
    print("computing PCA transformations")
    if (external_field_path == ""){
      env_df <- load_env_data(model_predictors, env_path = env_path, decimalCoords = FALSE)
      env_df <- env_df %>%
        dplyr::arrange(Latitude, Longitude)
    } else {#Load external predictor fields
      env_df <- read.csv(external_field_path)
      #Also append out_path then
      if (!all(c(model_predictors, "Longitude", "Latitude", "Month") %in% colnames(env_df))){
        print("not all columns present in the environmental data field")
        print(colnames(env_df))
        print("using training fields instead")
        env_df <- load_env_data(model_predictors, env_path = env_path, decimalCoords = FALSE)
      } else {
        env_df <- env_df %>%
          dplyr::select(all_of(c(model_predictors, "Longitude", "Latitude", "Month"))) %>%
          dplyr::arrange(Latitude, Longitude)
      }
    }

      env_df_dims <- as.data.frame(predict(res.pca, newdata = env_df))
      #Rename columns to "Dim.1" etc.
      colnames(env_df_dims) <- stringr::str_replace(colnames(env_df_dims), "PC","Dim.")
      #Select only those that we want and append to env_df
      env_df <- cbind(env_df, env_df_dims[predictors]) %>%
        dplyr::select("Longitude", "Latitude", "Month", all_of(predictors))

  }

  #MESS analysis if wanted
  if (MESS == T){
    print("calculating MESS")
    #Create new dataset from training data with only relevant environmental columns
    mess_df_train <- as.data.frame(h2o_train) %>%
      dplyr::select(all_of(predictors))

    #Create raster file from environmental data
    env_df_reloc <- env_df %>%
      dplyr::relocate(Longitude, Latitude)
    env_df_reloc["MESS"] = NA

    #Make a new data set following the shape of env_df_reloc
    env_df_full <- env_df_reloc[0,]


    #Loop over months and append
    for (m in unique(env_df_reloc$Month)){
      raster_env_df <- env_df_reloc %>%
        dplyr::filter(Month == m) %>%
        dplyr::select(Longitude, Latitude, all_of(predictors)) %>%
        raster::rasterFromXYZ()

      #MESS analysis
      #time_start <- Sys.time()
      mess_out <- suppressWarnings(dismo::mess(raster_env_df, mess_df_train, full = FALSE))
      #print(Sys.time()-time_start)

      #Raster dataset to points
      mess_out <- as.data.frame(raster::rasterToPoints(mess_out))
      #Replace inf with nan
      mess_out[mapply(is.infinite, mess_out)] <- NA
      #Append month variable
      mess_out["Month"] = m
      #Reorder
      mess_out <- mess_out %>%
        dplyr::arrange(y, x)

      #Append to full dataset
      env_df_reloc[env_df_reloc$Month == m, "MESS"] = mess_out$mess
    }
    env_df <- env_df_reloc
    env_df["InsideTraining"] = env_df_reloc$MESS >= 0

  }

  #Depending on the projection, drop values outside of the desired latitudinal range
  load(basin_path)
  env_df_na <- basin_df %>%
    dplyr::filter(Depth == 200) %>%
    dplyr::select(-Basin) %>%
    dplyr::right_join(., env_df, by = c("Latitude", "Longitude"))

  #Merge basin names, CPR_regions and gridcellsize at 200 m to env_df
  if (region == "full"){
    #Mask land points
    env_df_na <- env_df_na %>%
      dplyr::filter(!is.na(basin_name) & !(basin_name %in% c("SuluSea", "RedSea", "MediterraneanSea",
                                                             "BalticSea", "HudsonBay", "KaraSea")))

  } else if (region == "SO"){
    #Cut off at 40°S -> only values south of this
    env_df_na <- env_df_na %>%
      dplyr::filter(CPR_regions == "SO_CPR")
  } else if (region == "NA"){
    #North Atlantic model
    env_df_na <- env_df_na %>%
      dplyr::filter(CPR_regions == "NA_CPR")
  } else if (region == "NP"){
    #North Pacific model
    env_df_na <- env_df_na %>%
      dplyr::filter(CPR_regions == "NP_CPR") %>%
      #Change the longitude column for better plotting
      dplyr::mutate(Longitude = ifelse(Longitude >= 0, Longitude - 360, Longitude))
  }

  #Drop the basin name columns again
  env_df_na <- env_df_na %>%
    dplyr::select(-basin_name, -CPR_regions, -Depth) %>%
    #Drop all points where one of the predictors is NaN
    tidyr::drop_na(any_of(predictors))


  #Predict target variable for entire var space
  for (i in seq_along(model_list)){
    m <- model_list[i][[1]]
    m_nam <- model_name_list[i]
    env_df_h2o <- env_df_na %>%
      dplyr::select(all_of(predictors)) %>%
      h2o::as.h2o(.)

    #For GLM, add polynomial terms
    if (m_nam == "GLM"){
      #Adding polynomial terms
      h2o_pred_glm <- env_df_h2o
      predictors_glm <- predictors
      for (p in predictors){
        p_nam = paste(p, "^2", sep = "")
        predictors_glm <- append(predictors_glm, p_nam)
        h2o_pred_glm[p_nam] <- h2o_pred_glm[p]^2
      }
      #Prediction for GLM
      pr_model <- h2o::h2o.predict(object = m, newdata = h2o_pred_glm)
    } else {
      #For all other model types, regular prediction for H2O models
      pr_model <- h2o::h2o.predict(object = m, newdata = env_df_h2o)
    }
    #Append prediction column to env_df
    env_df_na[m_nam] = as.data.frame(pr_model)$predict
  }


  #Pivot model predictions to long form and calculate standard deviation between models
  final_pred <- env_df_na %>%
    #Calculate standard deviation between model predictions
    dplyr::mutate(target_std_abs = matrixStats::rowSds(as.matrix(env_df_na[model_name_list]), na.rm = TRUE),
    #Add the size of the grid cell in m2 to later calculate absolute standing stocks
           gridcellsize_m2 = area_lat(Latitude)) %>%
    #Pivot to long form
    tidyr::pivot_longer(cols = all_of(model_name_list),
                 names_to = "Model",
                 values_to = target)

  ### ==== Baseline plot ==========
  #Plotting preparation
  world <- ggplot2::map_data("world")

  #Make baseline plot with map etc. -> make these before std and MESS plots
  p_base <- function(df = final_pred){
    if (region == "full"){
      p_base <- ggplot2::ggplot(df) +
        ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, group = group),
                     data = world[world$long <= 180,], fill = "grey85", colour = "grey50", size = 0.3) +
        ggplot2::scale_x_continuous(name = "", breaks = c(-180,-120,-60,0,60,120,180),
                           labels = c("180°","120°W","60°W","GM","60°E","120°E","180°"),
                           expand = c(0,0)) +
        ggplot2::scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
                           labels = c("-90°N","-60°N","-30°N","Eq","30°N","60°N","90°N"), expand = c(0,0)) +
        ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),legend.key = ggplot2::element_rect(fill = "grey50"),
              panel.grid.major = ggplot2::element_line(colour = "white",linetype = "dashed")) +
        ggplot2::coord_quickmap() +
        ggplot2::theme(legend.position = "bottom",
              axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
    } else if (region == "SO"){
      #For polar projections
      # Defines the x axes required
      x_lines <- seq(-120,180, by = 60)
      #Create extra data.frames for lines
      lines_df <- data.frame(x_1 = 180,
                             y_1 = seq(-45, -85, by = -10),
                             label_1 = paste0(seq(-45, -85, by = -10), "°N"))

      xlines_df <- data.frame(x_lines = x_lines,
                              y = -32,
                              label =  c("120°W", "60°W", "0°", "60°E", "120°E", "180°W") )

      p_base <- ggplot2::ggplot(df) +
        ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, group = group),
                     data = world[(world$long <= 180) & (world$lat <= -40),], fill = "grey85", colour = "grey50", size = 0.3) +
        # Convert to polar coordinates
        ggplot2::coord_map("ortho", orientation = c(-90, 0, 0)) +
        ggplot2::scale_y_continuous(breaks = seq(-40, -90, by = -5), labels = NULL) +
        # Removes Axes and labels
        ggplot2::scale_x_continuous(breaks = NULL) +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        # Adds labels
        ggplot2::geom_text(ggplot2::aes(x = x_1, y = y_1, hjust = -0.2, label = label_1), data = lines_df, size = 2.5) +
        ggplot2::geom_text(ggplot2::aes(x = x_lines, y = y, label = label), data = xlines_df, size = 2.5) +
        # Adds axes
        ggplot2::geom_hline(ggplot2::aes(yintercept = -40), size = 1)  +
        ggplot2::geom_segment(ggplot2::aes(y = -40, yend = -90, x = x_lines, xend = x_lines), data = xlines_df, linetype = "dashed") +
        # Change theme to remove axes and ticks
        ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
              panel.grid.major = ggplot2::element_line(size = 0.25, linetype = 'dashed',
                                              colour = "black"),
              axis.ticks = ggplot2::element_blank(),
              legend.key = ggplot2::element_rect(fill = "grey50"),
              legend.position = "bottom")

    } else if (region == "NA"){
      #"NA_CPR", ifelse((((basin_name == "PacificOcean") | (basin_name == "ArcticOcean")) & ((Latitude > 40) & (Latitude < 65)) & ((Longitude < - 100) | (Longitude > 150))),
      #Make north Atlantic map
      (p_base <- ggplot2::ggplot(df) +
          ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, group = group),
                data = world,
                       fill = "grey85", colour = "grey50", size = 0.3) +
          ggplot2::coord_quickmap(xlim = c(-100, 30), ylim = c(40, 70)) +
          ggplot2::scale_x_continuous(name = "", breaks = c(-90,-60, -30, 0, 30),
                labels = c("90°","60°W","30°W","GM","30°E")) +
          ggplot2::scale_y_continuous(name = "", breaks = seq(40, 70, by = 5),
                labels = paste0(seq(40, 70, by = 5), "°N")) +
          ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),legend.key = ggplot2::element_rect(fill = "grey50"),
              panel.grid.major = ggplot2::element_line(colour = "white",linetype = "dashed")) +
          ggplot2::theme(legend.position = "bottom",
              axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)))

    } else if (region == "NP"){
      world_pacific <- world[(world$long <= -120) | (world$long >= 120), ]
      world_pacific$long <- ifelse(world_pacific$long > 0, world_pacific$long - 360, world_pacific$long)
      #Make north Pacific map
      (p_base <- ggplot2::ggplot(df) +
        ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, group = group),
                     data = world_pacific,
                     fill = "grey85", colour = "grey50", size = 0.3) +
        ggplot2::coord_quickmap(xlim = c(-240, -120), ylim = c(38, 67)) +
        ggplot2::scale_x_continuous(name = "", breaks = c(-240, -210, -180, -150, -120),
                           labels = c("120°E", "150°E", "180°E/W","150°W","120°W")) +
        ggplot2::scale_y_continuous(name = "", breaks = seq(40, 65, by = 5),
                           labels = paste0(seq(40, 65, by = 5), "°N")) +
        ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", size = 5),legend.key = ggplot2::element_rect(fill = "grey50"),
              panel.grid.major = ggplot2::element_line(colour = "white",linetype = "dashed")) +
        ggplot2::theme(legend.position = "bottom",
              axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)))
    }
    return(p_base)
  }

  #Create MESS stipling to add to plot if wanted
  if(MESS == T){
    #Set the resolution of stippling points depending on the number of facets
    if (length(model_name_list) > 1){
      res_stip = 6
    } else {res_stip = 4}

    #Differentiate between monthly resolution and yearly average
    #Keep every month, as this might change over the course of the year
      env_df_stiple_monthly <- env_df_na %>%
       dplyr::filter(InsideTraining == FALSE) %>%
        #Keep only some lat and lons to not make it too full
       dplyr::filter((((Latitude-0.5) %% res_stip == 0) & ((Longitude - 0.5) %% res_stip == 0)) | (((Latitude - 0.5)  %% res_stip == res_stip/2) & ((Longitude - 0.5) %% res_stip == res_stip/2)))

      #If we want yearly average plot, sum the values up over months and if
      #there is a six month or more value outside training, plot it
      env_df_stiple_yearly <- env_df_na %>%
        dplyr::group_by(Latitude, Longitude) %>%
        dplyr::summarise(InsideTraining = sum(as.integer(InsideTraining))) %>%
       dplyr::filter(InsideTraining <= 6) %>%
        #Keep only some lat and lons to not make it too full
       dplyr::filter((((Latitude-0.5) %% res_stip == 0) & ((Longitude - 0.5) %% res_stip == 0)) | (((Latitude - 0.5)  %% res_stip == res_stip/2) & ((Longitude - 0.5) %% res_stip == res_stip/2))) %>%
        dplyr::ungroup()

      #Make seasonally averaged stiple
    #If we want yearly average plot, sum the values up over months and if
    #there is a six month or more value outside training, plot it
    env_df_stiple_seasonal <- env_df_na %>%
      dplyr::mutate(Season = factor(Month, levels = c(1:12), labels = c(rep("DJF", times = 2), rep("MAM", times = 3),
                                                                 rep("JJA", times = 3), rep("SON", times = 3), "DJF"))) %>%
      dplyr::group_by(Latitude, Longitude, Season) %>%
      dplyr::summarise(InsideTraining = sum(as.integer(InsideTraining))) %>%
     dplyr::filter(InsideTraining < 3) %>%
      #Keep only some lat and lons to not make it too full
     dplyr::filter((((Latitude-0.5) %% res_stip == 0) & ((Longitude - 0.5) %% res_stip == 0)) | (((Latitude - 0.5)  %% res_stip == res_stip/2) & ((Longitude - 0.5) %% res_stip == res_stip/2))) %>%
      dplyr::ungroup()
    }


  ### ==== MESS, std plot and varimp plot =======
  #If wanted, make a mess plot
  if (mess_plot == T){
    print("mess plot")
    #MESS plot as binary: within or not within training df
    mess_plot_df <- env_df_na %>%
      tidyr::drop_na(InsideTraining) %>%
      #Sum up over months
      dplyr::group_by(Latitude, Longitude) %>%
      dplyr::summarise(InsideTraining = sum(as.numeric(InsideTraining)))

    #Create plot
    m_plot <- p_base() -
      ggplot2::geom_tile(ggplot2::aes(x = Longitude, y = Latitude, fill = InsideTraining),
                data = mess_plot_df) +
      viridis::scale_fill_viridis() +
      ggplot2::labs(fill = "Months with \nenvironmental space \nin training set")


    #Save the plots
    nam_mess <- paste(model_id, "_MESS_", Sys.Date(), ".png", sep = "")
    nam_mess_out <- paste(out_path, nam_mess, sep = "/")
    ggplot2::ggsave(nam_mess_out, m_plot,  width = 200, unit = "mm", height = 150, dpi = 300)

  }

  #Average nSD is Nan unless we calculate it
  ave_nSD <- c(NaN, NaN)
  #If desired, make map of standard deviation between models
  if (std_plot == T){
    print("std plot")
    #Calculate yearly average
    p_std <- final_pred %>%
      dplyr::group_by(Latitude, Longitude) %>%
      dplyr::summarise(target = mean(get(target), na.rm = TRUE),
                AbsoluteStd = mean(target_std_abs, na.rm = TRUE)) %>%
      dplyr::mutate(RelativeStd = AbsoluteStd/abs(target)*100)

    #Add yearly MESS stiple to plot if wanted
    if (MESS == TRUE){
      p_std_base <- p_base() +
        ggplot2::geom_point(ggplot2::aes(x = Longitude, y = Latitude),
                   data = env_df_stiple_yearly,
                   color = "black",
                   shape = 20,
                   alpha = 0.7,
                   size = 0.005)
    } else {p_std_base <- p_base()}

    #Plot of absolute SD based on p_base()
    p_std_abs <- p_std_base -
      ggplot2::geom_tile(ggplot2::aes(x = Longitude, y = Latitude, fill = AbsoluteStd),
                alpha = 1, data = p_std) +
      ggplot2::scale_fill_gradient2(low = "white",
                          high = "darkred", mid = "orange",
                          midpoint = mean(c(quantile(p_std$AbsoluteStd, 0.025, na.rm = TRUE),
                                    quantile(p_std$AbsoluteStd, 0.975, na.rm = TRUE))),
                          limits = c(0,
                                      quantile(p_std$AbsoluteStd, 0.975, na.rm = TRUE)),
                          oob = scales::squish) +
      ggplot2::theme(legend.position = "right",
            legend.title.align=0.5) +
      ggplot2::labs(fill = bquote(SD[abs]~(.(unit_small))))

    #Plot of relative SD based on p_base()
    p_std_rel <-  p_std_base -
      ggplot2::geom_tile(ggplot2::aes(x=Longitude, y=Latitude, fill=RelativeStd), alpha=1, data = p_std) +
      ggplot2::scale_fill_gradient2(low = "white",
                          high = "darkred", mid = "orange",
                          midpoint = mean(c(quantile(p_std$RelativeStd, 0.025, na.rm = TRUE),
                          quantile(p_std$RelativeStd, 0.975, na.rm = TRUE))),
                          limits = c(0,
                                     quantile(p_std$RelativeStd, 0.975, na.rm = TRUE)),
                          oob = scales::squish) +
      ggplot2::theme(legend.position = "right",
            legend.title.align=0.5) +
      ggplot2::labs(fill = bquote(SD[rel]~("%")))

    #Calculate average normalized nSD including MESS areas
    ave_nSD <- c(mean(p_std$RelativeStd)/100)
    #excluding MESS areas
    ave_nSD_excl_MESS <- final_pred %>%
      dplyr::filter(InsideTraining == TRUE) %>%
      dplyr::group_by(Longitude, Latitude) %>%
      dplyr::summarise(target = mean(get(target), na.rm = TRUE),
                       AbsoluteStd = mean(target_std_abs, na.rm = TRUE)) %>%
      dplyr::mutate(RelativeStd = AbsoluteStd/abs(target))
    ave_nSD <- c(ave_nSD, mean(ave_nSD_excl_MESS$RelativeStd))

    #Save plots
    nam_std <- paste0(model_id, "_abs_and_rel_sd_", Sys.Date(), "_.png", sep = "")
    out_nam_std <- paste(out_path, nam_std, sep = "/")
    g <- gridExtra::arrangeGrob(rbind(ggplotGrob(p_std_abs), ggplotGrob(p_std_rel), size="first"))
    ggplot2::ggsave(out_nam_std, g, width = 200, unit = "mm", height = 150)
  }


  ### ======= Calculate annual biomass stocks ================
  print("calculate annual biomass stocks")
  #Create copy of final_pred so we can leave out problematic values for flux calculation but leave them in for plotting
  final_pred_fluxes <- final_pred
  #Calculate total gridcellsize covered for weighted global average calculation
  total_monthly_gridcellsize_covered_m2 <- final_pred %>%
    dplyr::group_by(Month, Model) %>%
    dplyr::summarise(gridcellarea_monthly_m2 = sum(gridcellsize_m2))
    #Different definitions of total ocean area depending on plankton type and choice of area
  if (region == "full"){
      #total_ocean_area_m2_min <- 322 * 10^6 * 10^6 #Bednarsek (for pteropods, but we also use it for forams)
      total_ocean_area_m2 <- 362 * 10^6 * 10^6
    } else {
      #For all other cases, we just sum up the areas (this will include uncertainty, but still do it)
      total_ocean_area_m2 <- total_monthly_gridcellsize_covered_m2 %>%
        dplyr::group_by(Month) %>%
        dplyr::summarise(total_area = max(gridcellarea_monthly_m2)) %>%
        sum(.)
    }

    #If MESS is calculated, set all standing stocks to NA where we are outside of training dataset
    if (MESS == T){
      final_pred_fluxes[((final_pred_fluxes$InsideTraining == FALSE) & (!is.na(final_pred_fluxes$InsideTraining))), target] = NA
    }


    #Calculate max area covered by models
    area_covered <- final_pred_fluxes %>% dplyr::group_by(Model) %>%
      dplyr::summarise(total_area = sum(gridcellsize_m2)) %>%
      dplyr::select(total_area) %>%
      base::unique(.)
    print("total standing stock calculation")
    total_stocks_yearly_cellwise_means <- final_pred_fluxes %>%
      #Get gridcellwise mean daily flux across year, weighted by percentage of total ocean area
      dplyr::group_by(Latitude, Longitude, Model) %>%
      dplyr::summarise(standing_stock_unitsmall_TC_yearly_mean_m2 = mean(200 * (10^(get(target)) - 1), na.rm = TRUE),
                      percentage_area = gridcellsize_m2/area_covered[[1]]) %>%
      #Calculate global yearly mean daily flux per model
      dplyr::group_by(Model) %>%
      dplyr::summarise(standing_stock_unitsmall_TC_yearly_mean_global_sum = sum(standing_stock_unitsmall_TC_yearly_mean_m2 * percentage_area, na.rm = TRUE)) %>%
      #Calculate total magnitude of flux by multiplying global yearly mean daily flux with total ocean area and 365 days
      dplyr::mutate(standing_stock_unitlarge_TC_yearly_mean = standing_stock_unitsmall_TC_yearly_mean_global_sum * total_ocean_area_m2 * unit_conversion) %>%
      #Drop irrelevant columns
      dplyr::select(-standing_stock_unitsmall_TC_yearly_mean_global_sum) %>%
      #Create title columns
      dplyr::mutate(title = paste(Model, " - Mean stock: ", round(standing_stock_unitlarge_TC_yearly_mean,2), " ",
                                           unit_large, sep = ""))


  ### ===== Plot of the number of months that we have biomass estimates for =====
  print("final pred month cover")
  final_pred_month_cover <- final_pred %>%
    dplyr::left_join(env_df, by = all_of(c("Longitude", "Latitude", "Month"))) %>%
    dplyr::filter(!is.na(get(target))) %>%
    dplyr::group_by(Longitude, Latitude) %>%
    dplyr::summarise(n_months = n()/length(model_name_list))

  p_no_months <- p_base(final_pred_month_cover) +
    ggplot2::geom_tile(ggplot2::aes(x = Longitude, y = Latitude, fill = n_months)) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(fill = "# Months with predictions")

  #Save this
  nam <- paste0(model_id, "_no_months_predicted_", Sys.Date(), ".png")
  nam_out <- paste(out_path, nam, sep = "/")
  ggplot2::ggsave(nam_out, p_no_months, dpi = 300, width = 200,height = ifelse(region == "SO", 285, 150), units = "mm")

  #Save final_pred dataframe
  save(final_pred, file = paste(out_path, paste0(model_id, "_final_pred_df", Sys.Date(), ".RData"), sep = "/"))

  #### ====== Biomass plots ===========
  print("starting biomass plots")
  #Adapt color max and mins for prediction plot
  if (min_col ==""){
    min_col = stats::quantile(final_pred[,target][[1]], probs = 0.05, na.rm = TRUE)
  }
  if (max_col == ""){
    max_col = stats::quantile(final_pred[,target][[1]], probs = 0.99, na.rm = TRUE)
  }

  print("prediction plotting")
  #Facet wrap per model
  p_base_plot <- p_base() +
    ggplot2::facet_wrap(factor(Model, levels = model_name_list, labels = model_name_list)~.,
               ncol = 2) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 6),
          strip.placement = "outside",
          strip.background = ggplot2::element_rect(color = "white", size = 3))

  #For regular, non-monthly plot, make extra data frame with values averaged over months
  final_pred_year <- final_pred %>%
    dplyr::group_by(Latitude, Longitude, Model) %>%
    dplyr::summarise(target = mean(get(target), na.rm = TRUE))

  #Add yearly stiple if wanted
  if (MESS == TRUE){
    p_base_plot <- p_base_plot +
      ggplot2::geom_point(ggplot2::aes(x = Longitude, y = Latitude),
                 data = env_df_stiple_yearly,
                 color = "black",
                 shape = 20,
                 alpha = 0.7,
                 size = 0.005) +
      ggplot2::labs(fill = bquote(log[10]~Total~(.(unit_small))))
  }

  #Make plot of target variable
  p_pred <- p_base_plot -
      ggplot2::geom_tile(ggplot2::aes(x=Longitude, y=Latitude, fill=target), data = final_pred_year,
                alpha=1) +
      ggplot2::scale_fill_gradient2(limits = c(min_col, max_col),
        low = "cornflowerblue",
        mid = "yellow",
        high = "red",
        na.value = "black",
        midpoint = mean(c(min_col, max_col)),
        oob = scales::squish,
        guide = guide_colorbar(legend.title = target,
        title.position = "top", barheight = 1,
                               barwidth = 10))

  ### ==== Also make model-average plot ================
  #For regular, non-monthly plot, make extra data frame with values averaged over months
  final_pred_year_ave <- final_pred %>%
    dplyr::group_by(Latitude, Longitude) %>%
    dplyr::summarise(target = mean(get(target), na.rm = TRUE))

  p_average <- p_base(df = final_pred_year_ave) +
    ggplot2::geom_tile(ggplot2::aes(x=Longitude, y=Latitude, fill=target), data = final_pred_year,
              alpha=1) +
    ggplot2::scale_fill_gradient2(limits = c(min_col, max_col),
                         low = "cornflowerblue",
                         mid = "yellow",
                         high = "red",
                         na.value = "black",
                         midpoint = mean(c(min_col, max_col)),
                         oob = scales::squish,
                         guide = guide_colorbar(legend.title = target,
                                                title.position = "top", barheight = 1,
                                                barwidth = 10))

  #Add MESS if wanted
  if (MESS == TRUE){
    p_average <- p_average +
      ggplot2::geom_point(ggplot2::aes(x = Longitude, y = Latitude),
                 data = env_df_stiple_yearly,
                 color = "black",
                 shape = 20,
                 alpha = 0.7,
                 size = 0.005) +
      ggplot2::labs(fill = bquote(log[10]~Total~(.(unit_small)))  )
  }


  #Save plots if wanted
  if (save_pred_plot == T){
    p_shifted <- shift_legend(p_pred)
    nam <- paste(model_id, "_predictions_", Sys.Date(), ".png", sep = "")
    nam_out <- paste(out_path, nam, sep = "/")
    ggplot2::ggsave(nam_out, p_shifted, dpi = 300, width = 200,height = ifelse(region == "SO", 285, 150), units = "mm")

    #Model mean
    nam <- paste(model_id, "_predictions_model_average_", Sys.Date(), ".png", sep = "")
    nam_out <- paste(out_path, nam, sep = "/")
    ggplot2::ggsave(nam_out, p_average, dpi = 300, width = 200,height = ifelse(region == "SO", 285, 150), units = "mm")
  }

  ### ===== Make seasonally averaged plots if wanted ============
  if (seasonal == TRUE){
    #Add observation points on same colorscale if wanted
    if (add_points == TRUE){
      #Load dataset
      obs_val <- read.csv(paste(here::here(),data_path, sep = "/"))
      #Make it seasonal
      obs_val_season <- obs_val %>%
        dplyr::mutate(Season = factor(Month, levels = c(1:12), labels = c(rep("DJF", times = 2), rep("MAM", times = 3),
                                                                   rep("JJA", times = 3), rep("SON", times = 3), "DJF")))

    }

    final_pred_season <- final_pred %>%
      dplyr::mutate(Season = factor(Month, levels = c(1:12), labels = c(rep("DJF", times = 2), rep("MAM", times = 3),
                                                                 rep("JJA", times = 3), rep("SON", times = 3), "DJF"))) %>%
      dplyr::group_by(Longitude, Latitude, Season) %>%
      dplyr::summarise(target = mean(get(target), na.rm = TRUE))

      #rm(p_season)
    #windows()
    #Make seasonal plots with facet_wrap
    p_season <- p_base(df = final_pred_season) -
        ggplot2::geom_tile(ggplot2::aes(x=Longitude, y=Latitude, fill=target),
                  data = final_pred_season,
                  alpha=1) +
        ggplot2::facet_wrap(factor(Season, levels = c("DJF", "MAM", "JJA", "SON"),
                          labels = c("DJF", "MAM", "JJA", "SON"))~.) +
        ggplot2::scale_fill_gradient2(limits = c(min_col, max_col),
                             low = "cornflowerblue",
                             mid = "yellow",
                             high = "red",
                             #na.value = "black",
                             midpoint = mean(c(min_col, max_col)),
                             oob = scales::squish,
                             guide = guide_colorbar(legend.title = target,
                                                    title.position = "top", barheight = 1,
                                                    barwidth = 10))     +
        #ggtitle(s) +
        ggplot2::labs(fill = bquote(log[10]~Total~(.(unit_small)))    )

      if (MESS == TRUE){
        p_season <- p_season +
          ggplot2::geom_point(ggplot2::aes(x = Longitude, y = Latitude),
                     data = env_df_stiple_seasonal,
                     color = "black",
                     shape = 20,
                     alpha = 0.7,
                     size = 0.005)
      }
      nam <- paste0(model_id, "_predictions_model_average_seasonal_", Sys.Date(), ".png")
      out_nam <- paste(out_path, nam, sep = "/")
      ggplot2::ggsave(out_nam, p_season, dpi = 300, width = 200, height = ifelse(region == "SO", 285, 150), units = "mm")

      #Add points and save that version if we want
      if (add_points == TRUE){
        p_season <- p_season +
          ggplot2::geom_point(ggplot2::aes(x = lon_gridded, y = lat_gridded, fill = get(target)),
                     stroke = 0.5, shape = 21, size = 1,
                     data = obs_val_season)

        #Save this one too
        nam <- paste(model_id, "_predictions_model_average_seasonal_plus_obs_", Sys.Date(), ".png", sep = "")
        out_nam <- paste(out_path, nam, sep = "/")
        ggplot2::ggsave(out_nam, p_season, dpi = 300, width = 200, height = ifelse(region == "SO", 285, 150), units = "mm")


    }
  }


  ### ===== Make monthly model average plot if wanted ============
  if (monthly == TRUE){
    #Add observation points on same colorscale if wanted
    if (add_points == TRUE){
      #Load dataset
      obs_val <- read.csv(paste(here::here(),data_path, sep = "/"))
    }

    final_pred_monthly <- final_pred %>%
      dplyr::group_by(Longitude, Latitude, Month) %>%
      dplyr::summarise(target = mean(get(target), na.rm = TRUE))

    for (m in seq(1, 12)){
      p_month <- p_base(df = final_pred_monthly) -
        ggplot2::geom_tile(ggplot2::aes(x=Longitude, y=Latitude, fill=target),
                  data = final_pred_monthly[final_pred_monthly$Month == m,],
                  alpha=1) +
        ggplot2::scale_fill_gradient2(limits = c(min_col, max_col),
                             low = "cornflowerblue",
                             mid = "yellow",
                             high = "red",
                             na.value = "black",
                             midpoint = mean(c(min_col, max_col)),
                             oob = scales::squish,
                             guide = guide_colorbar(legend.title = target,
                                                    title.position = "top", barheight = 1,
                                                    barwidth = 10))     +
        ggtitle(month.name[m]) +
        ggplot2::labs(fill = bquote(log[10]~Total~(.(unit_small)))    )

      if (MESS == TRUE){
        p_month <- p_month +
          ggplot2::geom_point(ggplot2::aes(x = Longitude, y = Latitude),
                     data = env_df_stiple_monthly[env_df_stiple_monthly$Month == m, ],
                     color = "black",
                     shape = 20,
                     alpha = 0.7,
                     size = 0.005)
      }
      nam <- paste(model_id, "_predictions_model_average_monthly_", ifelse(m<10, paste(0, m, sep = ""), m), "_", Sys.Date(), ".png", sep = "")
      out_nam <- paste(out_path, nam, sep = "/")
      ggplot2::ggsave(out_nam, p_month, dpi = 300, width = 200, height = ifelse(region == "SO", 285, 150), units = "mm")

      #Add points and save that version if we want
      if (add_points == TRUE){
        p_month <- p_month +
          ggplot2::geom_point(ggplot2::aes(x = lon_gridded, y = lat_gridded, fill = get(target)),
                     stroke = 0.5, shape = 21, size = 1,
                     data = obs_val[obs_val$Month == m,])

        #Save this one too
        nam <- paste0(model_id, "_predictions_model_average_monthly_plus_obs_", ifelse(m<10, paste(0, m, sep = ""), m), "_", Sys.Date(), ".png")
        out_nam <- paste(out_path, nam, sep = "/")
        ggplot2::ggsave(out_nam, p_month, dpi = 300, width = 200, height = ifelse(region == "SO", 285, 150), units = "mm")

      }
    }

  }

  #Save the predictions
  if (MESS == T){
    #If MESS is calculated, exclude the values outside of training range from mean calculation
      final_pred[((final_pred$InsideTraining == FALSE) & (!is.na(final_pred$InsideTraining))), target] = NA
    }

  #Create dataframe to save
  final_pred_save <- final_pred %>%
    dplyr::group_by(Model, Month) %>%
    dplyr::summarise(target_val = mean(get(target), na.rm = TRUE)) %>%
    #Add target, predictors, model_id, dataset_name, outliers, weight, weight_column  column
    dplyr::mutate(target_name = target,
           predictors = rep(as.character(paste(predictors, collapse = ", "), times = n())),
           model_id = model_id,
           df_choice = df_choice,
           remove_outliers = remove_outliers,
           weight_option = weight_option,
           region = region)

  out_csv <- paste(out_path, paste(model_id, "_prediction_values_", Sys.Date(), ".csv", sep = ""), sep = "/")
  write.csv(final_pred_save, file = out_csv)

  #Save the yearly stock calculations
  total_stocks_yearly <- total_stocks_yearly_cellwise_means %>%
      #Add target, predictors, model_id, dataset_name, outliers, weight, weight_column  column
      dplyr::mutate(target = target,
             predictors = rep(as.character(paste(predictors, collapse = ", "), times = n())),
             model_id = model_id,
             df_choice = df_choice,
             remove_outliers = remove_outliers,
             weight_option = weight_option,
             region = region)

    out_csv <- paste(out_path, paste0(model_id, "_total_stocks_", Sys.Date(), ".csv"), sep = "/")
    write.csv(total_stocks_yearly, out_csv, row.names = FALSE)


  #Return average nSD
  return(ave_nSD)
}


