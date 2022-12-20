#Load the prediction RData files for a number of model runs and calculate the contribution to final predictions of each
#setup parameter for one region
###TODO: set this up as function call for model_pipeline
#Parameters: predictor choice, model algorithm choice, w/wo zeros, w/wo outliers,
library(here)
library(h2o)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

#Testing
#Choose date
#d_nam <- Sys.Date() - 1
##Choose plankton type
#plankton_type <- "foraminifera"
#
##Define ANOVA factors
#anova_fac <- c("df_choice", "Model", "predictors", "TIC_TC_choice", "growth_rate_choice")
#anova_names <- c("Dataframe choice", "Model choice", "Predictor choice",  "TIC-TC factor", "Growth rate choice")
#

uncertainty_contributions <- function(plankton_type,
                                      in_path,
                                      anova_fac = c("df_choice", "Model", "predictor_options", "TIC_TC_choice", "growth_rate_choice"),
                                      anova_names = c("Dataframe choice", "Model choice", "Predictor choice",  "TIC-TC factor", "Growth rate choice"),
                                      d_nam = Sys.Date()){

  #' Function to run a multivariate ANOVA and assess contribution of different factors to inter-model differences
  #'
  #' @param plankton_type       Character, either "foraminifera", or "pteropoda"
  #' @param anova_fac           Vector of names from the prediction_values csv to be used for this analysis
  #' @param anova_names         Vector of characters to be used as names for the anova_fac in the plotting
  #' @param d_nam               Date input to specify which prediction_value files should be used for analysis, defaults to today's date

  #'
  #' @return will save two plots: independent boxplots for the difference of the C_fluxes between different levels of the
  #' anova_fac; and effect strength of the different factors and their interactions as calculated by MANOVA
  #' @export
  #'
  #' @examples

  #Load relevant files
  #in_path <- here::here("./Models")

  #Get list of files and load the relevant RData files
  lf <- list.files(paste(in_path))
  #Load all those relevant ones
  pred_files <- sort(lf[grepl(paste("_prediction_values_", d_nam, sep = ""), lf, fixed = TRUE)])
  p_df <- lapply(paste(in_path, pred_files, sep = "/"), read.csv) %>%
    dplyr::bind_rows() %>%
    dplyr::select(-X)

  #Rename the columns
  for (c in anova_fac){
    colnames(p_df)[colnames(p_df) == c] = anova_names[anova_fac == c]
  }


  #Do the next plots for all the regions analysed in the analyses
  regions_unique <- unique(p_df$region)
  for (re in regions_unique){
    print(re)

    #Check if we have different datasets and predictor sets, otherwise, leave this out
    if (length(unique(p_df$`Dataframe choice`)) == 1){
      anova_fac = anova_fac[anova_fac != "df_choice"]
      anova_names = anova_names[anova_names != "Dataframe choice"]
    }
    if (length(unique(p_df$predictors)) == 1){
      anova_fac = anova_fac[anova_fac != "predictor_options"]
      anova_names = anova_names[anova_names != "Predictor choice"]
    }

    #First, for full dataset create a boxplot for each of the influencing factors
    p_plot <- p_df %>%
      dplyr::filter(region == re) %>%
      dplyr::mutate("Predictor choice" = factor(predictors, levels = unique(p_df$predictors), labels = c(1:length(unique(p_df$predictors))))) %>%
      #For the predictor column, add a linebreak
      #dplyr::mutate(`Predictor choice` = str_replace_all(`Predictor choice`, ",\\s*", ",\n")) %>%
      tidyr::pivot_longer(cols = all_of(anova_names),
                   names_to = "anova_facs",
                    values_to = "anova_vals")

    #Make plot
    s_boxplot <- p_plot %>%
      ggplot2::ggplot(ggplot2::aes(x = factor(anova_vals, levels = c(unique(p_plot[p_plot$anova_facs == "Predictor choice",]$anova_vals), unique(p_df$`Dataframe choice`),
                                                                                      "min", "mean", "max",
      "GLM", "GAM", "RF", "GBM", "DL"), labels = c(unique(p_plot[p_plot$anova_facs == "Predictor choice",]$anova_vals), unique(p_df$`Dataframe choice`),
                                                                                      "min", "mean", "max",
                                                                                      "GLM", "GAM", "RF", "GBM", "DL")),
                 y = C_flux_daily)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::facet_wrap(factor(anova_facs, levels = anova_names, labels = anova_names)~., scales = "free_x") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
      ggplot2::labs(x = "") +
      ggplot2::theme_bw() +
      ylim(0, quantile(p_plot$C_flux_daily, 0.95, na.rm = TRUE))

    #multiway ANOVA: predictor choice, model
    p_anova <- p_df %>%
      dplyr::mutate(`Predictor choice` = factor(predictors)) %>%
      dplyr::select(all_of(anova_names), "C_flux_daily")

    manova_model <- aov(as.formula(paste0("C_flux_daily ~ `",
                                          paste(anova_names, collapse = "` * `"), "`")),
                        data = p_anova)

    s_ma <- as.data.frame(summary(manova_model)[[1]])

    #Extract coefficients and show as percentage of total SS
    s_ma <- s_ma %>%
      dplyr::mutate(sum_sq = `Sum Sq`,
             percentage_explained = sum_sq/sum(sum_sq)) %>%
      dplyr::arrange(-percentage_explained) %>%
      dplyr::filter(`Pr(>F)`<0.05) %>%
      #Remove trailing whitespaces and "`" signs
      dplyr::mutate(effect_factor = str_trim(str_replace_all(rownames(.), "`", "")))

    #Plot effect strength
    s_effect <- s_ma %>%
      ggplot2::ggplot(ggplot2::aes(x = reorder(effect_factor, percentage_explained), y = 100*percentage_explained)) +
      ggplot2::geom_segment(ggplot2::aes(xend = reorder(effect_factor, percentage_explained), yend=0)) +
      ggplot2::geom_point( size = 2, color="black") +
      ggplot2::theme_bw() +
      ggplot2::coord_flip() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 1)) +
      ggplot2::labs(y = "Explained variability (%)", x = "" ) +
      ylim(0, 50)

    #Save the plots
    out_name_boxplot <- paste(in_path, paste0(re, "_confounder_effects_boxplot_", Sys.Date(), ".png"), sep = "/")
    out_name_effect <- paste(in_path, paste0(re, "_confounder_effects_manova_", Sys.Date(), ".png"), sep = "/")

    ggplot2::ggsave(out_name_boxplot, s_boxplot, width = 200, dpi = 300, unit = "mm")
    ggplot2::ggsave(out_name_effect, s_effect, width = 100,height = 50,  dpi = 300, unit = "mm")

  }

}

#uncertainty_contributions(plankton_type = "foraminifera", d_nam = Sys.Date() - 1)