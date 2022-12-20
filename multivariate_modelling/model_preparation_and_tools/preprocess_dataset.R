library(dplyr)
library(here)
library(h2o)
library(tidyr)
library(factoextra)
library(caret)

preprocessing_modelling <- function(predictors = c(),
                          target = c(),
                          data_path = "",
                          pca = FALSE,
                          region = "full",
                          basin_df_path = here::here("./General_Pipeline_PCA/Data/basin_data/basin_data_regional.RData"),
                          remove_outliers = TRUE,
                          weights_column = FALSE,
                          weight_option = "log"){

  #' Loads and preprocesses the dataset for multivariate modelling
  #'
  #' @param predictors      vector of names of predictors to use
  #' @param target          vector of name of target variable
  #' @param data_path       path to the training dataset, from here() directory down
  #' @param pca             Bool, if True, compute PCA-transformation of the predictors
  #' @param region          Character, "full", "SO", "NA", "NP" depending on whether we want a focus on one of the CPR regions
  #' @param basin_df_path   path to the basin file that contains SO, NA, NP designations
  #' @param remove_outliers Bool, if True, remove outliers of target variable with z-score criterion (z>3)
  #' @param weights_column  Bool, if True, use a column as weights -> not recommended
  #' @param weight_option   Character out of various options how to calculate the weights, needs to be adapted if you want to use it
  #'
  #' @return                Returns a list of h2o_train, h2o_test, and eventually the PCA-element and PCA-transformed predictor name vector

  ### ======= Load data =============
  df = read.csv(paste(here::here(),data_path, sep = "/"))

  #Load basin file
  load(here::here(basin_df_path))

  ### ======= Preprocess data ===========
  #Replace -Inf or Inf with NaN
  df_model <- df %>%
    dplyr::select(c(all_of(predictors), all_of(target),
                    "lat_gridded", "lon_gridded", "Month")) %>%
    dplyr::mutate_if(is.numeric, list(~na_if(., Inf))) %>%
    dplyr::mutate_if(is.numeric, list(~na_if(., -Inf))) %>%
    na.omit()

  #Merge the basin file
  df_model_basin <- basin_df %>%
    dplyr::filter(Depth == 200) %>%
    dplyr::select(-Basin) %>%
    dplyr::mutate(lat_gridded = Latitude,
           lon_gridded = Longitude) %>%
    dplyr::select(-Latitude, -Longitude) %>%
    right_join(., df_model, by = c("lat_gridded", "lon_gridded"))

  #Merge basin names, CPR_regions and gridcellsize at 200 m to env_df
  if (region == "full"){
    #Mask land points
    df_model_basin <- df_model_basin %>%
      dplyr::filter(!is.na(basin_name) & !(basin_name %in% c("SuluSea", "BayofBengal", "RedSea", "MediterraneanSea",
                                                             "BalticSea", "HudsonBay", "KaraSea")))
  } else if (region == "SO"){
    #Cut off at 40°S -> only values south of this
    df_model_basin <- df_model_basin %>%
      dplyr::filter(CPR_regions == "SO_CPR")
  } else if (region == "NA"){
    #North Atlantic model
    df_model_basin <- df_model_basin %>%
      dplyr::filter(CPR_regions == "NA_CPR")
  } else if (region == "NP"){
    #North Pacific model
    df_model_basin <- df_model_basin %>%
      dplyr::filter(CPR_regions == "NP_CPR") %>%
      #Change the longitude column for better plotting
      dplyr::mutate(Longitude = ifelse(lon_gridded >= 0, lon_gridded - 360, lon_gridded))
  }

  #Drop the basin columns
  df_model <- df_model_basin %>%
    dplyr::select(-basin_name, -CPR_regions, -Depth)

  #If we want to add weights as inverse to the sampling density, do it
  if (weights_column == TRUE){
    if (weight_option == "invlogcounts"){
      df_model <- df_model %>%
        dplyr::mutate(model_weights = 1/log10(Counts + 1))
    } else if (weight_option == "invorigcounts"){
      df_model <- df_model %>%
        dplyr::mutate(model_weights = 1/Counts)
    } else if (weight_option == "invsqrtcounts"){
      df_model <- df_model %>%
        dplyr::mutate(model_weights = 1/sqrt(Counts))
    } else if (weight_option == "origcounts"){
      df_model <- df_model %>%
        dplyr::mutate(model_weights = Counts)
    } else if (weight_option == "logcounts"){
      df_model <- df_model %>%
        dplyr::mutate(model_weights = log10(Counts + 1))
    } else if (weight_option == "sqrtcounts"){
      df_model <- df_model %>%
        dplyr::mutate(model_weights = sqrt(Counts))
    } else if (weight_option == "invorigsd"){
      df_model <- df_model %>%
        dplyr::mutate(model_weights = 1/(TC_sd+1))
    } else if (weight_option == "invlogsd"){
      df_model <- df_model %>%
        dplyr::mutate(model_weights = 1/log10(TC_sd + 2))
    } else if (weight_option == "invsqrtsd"){
      df_model <- df_model %>%
        dplyr::mutate(model_weights = 1/sqrt(TC_sd + 1))
    }
  } else {
    df_model <- df_model %>%
      dplyr::mutate(model_weights = 1)
  }
  #Remove outliers based on z-score if desired
  if (remove_outliers == TRUE){
    #nrow_orig <- nrow(df_model)
    df_model <- df_model %>%
      dplyr::mutate(zscore = (get(target) - mean(get(target), na.rm = TRUE))/
      sd(get(target), na.rm = TRUE),
      outlier_bool = zscore > 3) %>%
      dplyr::filter(outlier_bool == FALSE) %>%
      dplyr::select(-zscore, -outlier_bool)
    #print(paste(nrow_orig - nrow(df_model), "datapoints removed as outliers"))
  }

  if (pca == TRUE){#Compute PCA values for the points with plankton observations
    #Replace infinite values with nan
    df_model <- do.call(data.frame,                      # Replace Inf in data by NA
                  lapply(df_model,
                         function(x) replace(x, is.infinite(x), NA)))

    #Compute PCA on the chosen environmental predictors
    res.pca <- stats::prcomp(df_model[predictors], scale = TRUE)

    #Compute Scree plot
    scree_plot <- factoextra::fviz_eig(res.pca)
    print(scree_plot)
    #Also save this
    if (ifelse(!dir.exists(file.path(here::here(), str_remove(out_path, here::here()))),
               dir.create(file.path(here::here(), str_remove(out_path, here::here()))), FALSE) == TRUE &
      ifelse(!dir.exists(file.path(out_path, "Models")),
             dir.create(file.path(out_path, "Models")), FALSE) == TRUE){
      print("output folders created")
    } else {"problem with the output folders"}

    out_name <- paste0(out_path, "/Models/scree_plot_", Sys.Date(), ".png", sep ="")
    ggsave(out_name, scree_plot, width = 200, height = 150, dpi = 300, unit = "mm")

    #Compute contributions of different predictors to dimensions
    cont_plot <- factoextra::fviz_pca_var(res.pca,
                              col.var = "contrib", # Color by contributions to the PC
                              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                              repel = TRUE     # Avoid text overlapping
    )
    print(cont_plot)
    #Also save this
    out_name <- paste0(out_path, "/Models/contribution_pca_plot_", Sys.Date(), ".png", sep ="")
    ggsave(out_name, cont_plot, width = 200, height = 150, dpi = 300, unit = "mm")


    #Use those dimensions that explain >=80%
    eig.val <- factoextra::get_eigenvalue(res.pca) #eigenvalues
    eig.val <- eig.val %>%
      dplyr::mutate(PCA_Dim = rownames(.))
    rownames(eig.val) <- 1:nrow(eig.val)
    #Get first row where cumulative variance explained is >80
    max_dim <- min(as.integer(rownames(eig.val[eig.val$cumulative.variance.percent >=80,])))
    #Get those dimensions
    print(paste("using", max_dim, "dimensions, as they explain", eig.val[max_dim, "cumulative.variance.percent"], "of the variance"))
    pca_dims <- eig.val[1:max_dim,]
    pca_predictors <- pca_dims$PCA_Dim

    #Return PCA coordinates for those dimensions for all environmental observations
    res.ind <- factoextra::get_pca_ind(res.pca)

    # Results for Variables
    res.var <- factoextra::get_pca_var(res.pca)
    #Plot heatmap of coordinates
    res.heat <- as.data.frame(res.var$coord) %>%
      mutate(Variable = rownames(.)) %>%
      pivot_longer(cols = all_of(colnames(.)[!(colnames(.) %in% c("ranking","Variable"))]),
                   names_to = "pca_dimension",
                   values_to = "dimension_value") %>%
      mutate(ranking = as.integer(substr(pca_dimension, 5, 7)))

    #Order of ranking
    res.heat.order <- res.heat %>%
      dplyr::select(pca_dimension, ranking) %>%
      dplyr::arrange(ranking) %>%
      dplyr::distinct()

    res.heat.plot <- ggplot(aes(x = factor(pca_dimension, levels = res.heat.order$pca_dimension,
                                           labels = res.heat.order$pca_dimension),
                                y = Variable), data = res.heat) +
      geom_tile(aes(fill = dimension_value)) +
      scale_fill_gradient2(low = "red", high = "blue", mid = "white") +
      labs(x = "PCA Dimension", y = "", fill = "Coordinate") +
      ggplot2::geom_text(ggplot2::aes(label = round(dimension_value, 2)))
    #Also save this
    out_name <- paste0(out_path, "/Models/pca_variable_heatmap.png", Sys.Date(), ".png", sep ="")
    ggsave(out_name, res.heat.plot, width = 200, height = 150, dpi = 300, unit = "mm")

    #Append to env_df
    for (p in pca_predictors){
      df_model[p] <- as.data.frame(res.ind$coord)[p]
    }

  } else{#Keep empty/original elements
    res.pca <- NA
    pca_predictors <- predictors
  }

  # Train-test split
  set.seed(123)
  df_model$id <- 1:nrow(df_model)
  df_train <- df_model %>% dplyr::slice_sample(n = .75*nrow(df_model))
  df_test  <- dplyr::anti_join(df_model, df_train, by = 'id')
  #Convert variables to H2O objects
  h2o_train <- h2o::as.h2o(df_train)
  h2o_test <- h2o::as.h2o(df_test)

  #Return h2o_train and h2o_test and the name of the dataframe to set as dataset_before
  return(list(h2o_train, h2o_test, res.pca, pca_predictors))
}
