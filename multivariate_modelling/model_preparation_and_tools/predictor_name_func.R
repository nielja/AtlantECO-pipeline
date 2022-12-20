library(envnames)

pred_name_func <- function(predictors, use_pred_names){
  #' Function to assign predictor names as bquote expressions with proper layout.
  #'
  #' This function assigns layouted bquote expressions to the environmental predictors. Currently it only works for
  #' temperature, chlorophyll, MLD, EKE and oxygen.
  #'
  #' @param predictors      vector of predictor names for which to return the expressions
  #' @param use_pred_names  Bool, if True, we try to assign these expressions, otherwise just return the predictors
  #'
  #' @returns               A list, first element is the original predictor names, second is the bquote expressions

  if (use_pred_names == TRUE){
    print(paste0(envnames::environment_name(environment()), ":", envnames::address(environment())))
    #If no specific name list for the predictors is defined, assign the regular names
    predictor_names_df <- c((Temperature~(degree*C)),
                            (log[10]~~Chlorophyll*"-"*a~(mg~m^{-3})),
                            (log[10]~~MLD~(m)),
                            (log[10]~~EKE~(m^{2}~s^{-2})),
                            (Oxygen~(µmol~kg^{-1})))
    names(predictor_names_df) <-  c("Temperature_MLD", "log_Chl_a", "log_MLD",
                                    "log_EKE", "Oxygen_depth_200")
    if (all(predictors %in% names(predictor_names_df))){
      predictor_name_list <-  list()
      reordered_p <- list()
      for (p in names(predictor_names_df)){
        if (p %in% predictors){
          predictor_name_list <- append(predictor_name_list,
                                        predictor_names_df[p])
          reordered_p <- append(reordered_p, p)
        }
      }
      predictor_name_list <- unlist(predictor_name_list)
      predictors <- unlist(reordered_p)
    } else {
      predictor_name_list = as.list(predictors)
  }} else {predictor_name_list = as.list(predictors)}

  names(predictor_name_list) <- as.character(c(1:length(predictors)))

  #Return both the reordered predictors and the predictor_name_list
  return(list(predictors, predictor_name_list))
}