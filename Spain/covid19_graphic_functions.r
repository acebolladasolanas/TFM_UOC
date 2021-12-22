#**************************************************************************************************
#
#  This script contains some common functions and variable defintions used to create
#  graphics related to the analysis of the development of the COVID-19 in different countries.
#
#  The idea is to separate the steps to create the graphics, for which  specific functions are
#  defined here, from the statistical and predictive analysis of the data.
#
#  The execution of the functions defined here requires some objects that must be available
#  in the script, in which the funciotns are called.
#
#**************************************************************************************************

suppressWarnings(suppressMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(rstan)
  library(data.table)
  library(lubridate)
  library(gdata)
  library(EnvStats)
  library(matrixStats)
  library(scales)
  library(gridExtra)
  library(ggpubr)
  library(bayesplot)
  library(cowplot)
}))


#**************************************************************************************************
#  Functions to create graphics
#**************************************************************************************************

#--------------------------------------------------------------------------------------------------
#
#  Function geom_stepribbon
#  ------------------------
#
#  Function defined in a separate file
#
#--------------------------------------------------------------------------------------------------
source("utils/geom-stepribbon.r")


#--------------------------------------------------------------------------------------------------
#
#  Function f_graph_prediction_vs_data
#  -----------------------------------
#
#  This function is used to create a graphic to visualize real data (as histogram) compared
#  to estimated data (as curve).
#
#
#   Arguments
#   ---------
#     pandemic_dates     list of dates to which the data are related (will be shonwn in the x-axis)
#     model_prediction   predicted data
#     reported_data      real data
#     data_str           string to identify the type of data being displayed (will  be shown in y-axis)
#
#
#   Return value
#   ------------
#     Graphic generated from the data in the function arguments
#
#--------------------------------------------------------------------------------------------------
f_graph_prediction_vs_data <-function(pandemic_dates, model_prediction, reported_data, data_str) {

  # --- Extracta data from the arguments -------------------------------------------------------- #
  N_DATES                    <- ncol(model_prediction)
  predicted_values           <- colQuantiles(model_prediction  [,1:N_DATES], probs=0.5)        # 21Mar21
  predicted_values_95_li     <- colQuantiles(model_prediction  [,1:N_DATES], probs=0.025)
  predicted_values_95_ui     <- colQuantiles(model_prediction  [,1:N_DATES], probs=0.975)
  predicted_values_50_li     <- colQuantiles(model_prediction  [,1:N_DATES], probs=0.25 )      # 21Mar21
  predicted_values_50_ui     <- colQuantiles(model_prediction  [,1:N_DATES], probs=0.75 )      # 21Mar21


  predicted_values_95        <- data.frame(pandemic_dates, predicted_values_95_li, predicted_values_95_ui)
  names(predicted_values_95) <- c("date", "cases_min", "cases_max")
  predicted_values_95$key    <- rep("95%", length(predicted_values_95$date))

  predicted_values_50        <- data.frame(pandemic_dates, predicted_values_50_li, predicted_values_50_ui)     # 21Mar21
  names(predicted_values_50) <- c("date", "cases_min", "cases_max")
  predicted_values_50$key    <- rep("50%", length(predicted_values_50$date)) 

  predicted_values_50_95     <- rbind(predicted_values_95, predicted_values_50)
  levels(predicted_values_50_95$key) <- c("95%", "50%")

  
  # --- Create and return graphic --------------------------------------------------------------- #
  ret_plot <- ggplot(predicted_values_95) +
           
              geom_line( aes(x=pandemic_dates, y=predicted_values), col="black", alpha=0.5 ) +
              
              geom_ribbon( data= predicted_values_50_95, 
                           aes(x=date, ymin=cases_min, ymax=cases_max, fill=key)) +
              
              geom_bar( aes(x=pandemic_dates, y=reported_data), 
                        fill="coral4", stat='identity', alpha=0.5) +
              
              labs( title="", x="", y=data_str ) +    # 21Mar21
              
              scale_x_date( date_breaks="months", labels=date_format("%e %b") ) +               # 6Mar21 Juan
              coord_cartesian(ylim = c(0,1.2*max(max(reported_data),max(predicted_values)))) +    # 21Mar21 to set max value in plot
              scale_fill_manual(name="", labels=c("50%", "95%"),values=c(alpha("deepskyblue4", 0.55),alpha("deepskyblue4", 0.45)) ) + 
              
              theme_pubr() +
              theme( axis.text.x = element_text(angle=45, hjust=1), legend.position="None") + 
              guides(fill=guide_legend(ncol=1))
              
  return(ret_plot)
}


#--------------------------------------------------------------------------------------------------
#
#  Function f_graph_create_plots
#  -----------------------------
#
#  This function is used to create a graphic to visualize real data (as histogram) compared
#  to estimated data (as curve).
#
#
#   Arguments
#   ---------
#     output_graphics_file  file in which the graphics will be saved
#     
#     The function uses additionally data stored in the environment loaded from a 
#     'Rdata' file generated during the processing of the script 'covid19_base.r'.
#
#
#   Return value
#   ------------
#     None
#     The grahics generated will be saved in the filed defined in the
#     function argument 'output_graphics_file'
#
#--------------------------------------------------------------------------------------------------
f_graph_create_plots <- function(output_graphics_file) {
 
  f_log( paste0(" ", Sys.time(), " Start creating graphics") )
  f_log( paste0(" ", Sys.time(), " Output directory  '", dirname (output_graphics_file),  "'") )
  f_log( paste0(" ", Sys.time(), " Output file       '", basename(output_graphics_file),  "'") )
  f_log(" ")

files.pdf <-   
list.files(path = "../albertocebolladasolanas/Documents/Master/TFM/COVID19model_Jun21_noResults_noRDS/results/",
         pattern="pdf")

files.rdata <- 
  c("Italy_covid19_base_mod_2021-12-05_19-32-03.Rdata",
    "Italy_covid19_base_mod_2021-12-05_20-47-35.Rdata",
    "Italy_covid19_base_2021-12-06_00-10-17.Rdata",
    "Italy_covid19_base_2021-12-04_00-40-13.Rdata")


  # --- Extract information from the data loaded from the data file provided as script agrument - #
  # pandemic_dates           <- covid19_country_data$date    # already defined in main script   4Mar21 Juan
  N_DATES                  <- length(pandemic_dates)
  
  model_fit_out            <- rstan::extract(covid19_model_fit)
  model_prediction_cases   <- model_fit_out$prediction_cases
  model_prediction_deaths  <- model_fit_out$prediction_deaths
  model_rt                 <- model_fit_out$Rt
  model_cum_prediction_cases   <- model_fit_out$cum_prediction_cases     # 21Mar21
  model_cum_prediction_deaths  <- model_fit_out$cum_prediction_deaths    # 21Mar21
  model_rt_eff                 <- model_fit_out$Rt_eff          # 21Mar21
  model_prediction_newhospitalized        <- model_fit_out$prediction_newhospitalized   # 27Mar21
  model_prediction_hospitalized           <-  model_fit_out$prediction_hospitalized   # 27Mar21   PENDING: if no such data in stan data, initizalize with zeroes
  model_prediction_IA14        <-  model_fit_out$prediction_IA14       # 1Apr21
  model_prediction_detected    <-  model_fit_out$prediction_detected    # 16Apr21
  model_detRatio               <-  model_fit_out$detRatio      # 17Apr21
                           
  # reported_cases           <- covid19_country_data$Cases   # already defined in main script   4Mar21 Juan
  # reported_deaths          <- covid19_country_data$Deaths  # already defined in main script   4Mar21 Juan
                        


  # TEMP     # 17Apr21
  f_log("Average detection ratio")
  f_log( mean(model_detRatio))
  f_log(" ")




  # --- Create plots ---------------------------------------------------------------------------- # 
  
  #     Cases (infections) - Estimated data vs. reported data
  plot_cases  <- f_graph_prediction_vs_data( pandemic_dates, model_prediction_cases,  
                                             reported_cases, "Daily number of infections" )
                          
  #     Deaths - Estimated data vs. reported data                         
  plot_deaths <- f_graph_prediction_vs_data( pandemic_dates, model_prediction_deaths, 
                                             reported_deaths, "Daily number of deaths" )


 #     Cumulative cases (infections) - Estimated data vs. reported data             # 21Mar21
  plot_cum_cases  <- f_graph_prediction_vs_data( pandemic_dates, model_cum_prediction_cases,
                                             cum_reported_cases, "Cumulative number of infections" )

  #    Cumulative deaths - Estimated data vs. reported data                         # 21Mar21
  plot_cum_deaths <- f_graph_prediction_vs_data( pandemic_dates, model_cum_prediction_deaths,
                                             cum_reported_deaths, "Cumulative number of deaths" )


#     Daily new hosplitalized - Estimated data vs. reported data
  plot_newhospitalized  <- f_graph_prediction_vs_data( pandemic_dates, model_prediction_newhospitalized,
                                             reported_newhospitalized, "New hospitalized" )


#     Total hospitalized in a day - Estimated data vs. reported data
  plot_hospitalized  <- f_graph_prediction_vs_data( pandemic_dates, model_prediction_hospitalized,
                                             reported_hospitalized, "Total hospitalized in a given day" )


#     14-day cumulated incidence per 100K inhabitants - Estimated data vs. reported data
  plot_IA14  <- f_graph_prediction_vs_data( pandemic_dates, model_prediction_IA14,
                                             reported_IA14, "14-day cumulative incidence / 100K inhabitants" )
 
#     Detected cases (detected infections) - Estimated detected data vs. reported data
  plot_detected  <- f_graph_prediction_vs_data( pandemic_dates, model_prediction_detected,
                                             reported_cases, "Estimated detected cases" )






  #     Estimated values for 'Rt'
  
  # Plotting interventions
  predicted_rt       <- colQuantiles(model_rt[,1:N_DATES],probs=0.5)    # 21Mar21
  predicted_rt_95_li <- colQuantiles(model_rt[,1:N_DATES],probs=0.025)
  predicted_rt_95_ui <- colQuantiles(model_rt[,1:N_DATES],probs=0.975)
  predicted_rt_50_li <- colQuantiles(model_rt[,1:N_DATES],probs=0.25 )
  predicted_rt_50_ui <- colQuantiles(model_rt[,1:N_DATES],probs=0.75 )
  
  predicted_rt_95 <- data.frame(pandemic_dates, predicted_rt_95_li, predicted_rt_95_ui)
  names(predicted_rt_95) <- c("date", "rt_min", "rt_max")
  predicted_rt_95$key <- rep("95%", N_DATES)
  
  predicted_rt_50 <- data.frame(pandemic_dates, predicted_rt_50_li, predicted_rt_50_ui)
  names(predicted_rt_50) <- c("date", "rt_min", "rt_max")
  predicted_rt_50$key <- rep("50%", N_DATES)
  
  predicted_values_rt <- rbind(predicted_rt_95, predicted_rt_50)
  predicted_values_rt$key <- factor(predicted_values_rt$key)
  
  plot_rt <- ggplot(predicted_rt_95) +  
             geom_stepribbon( data = predicted_values_rt, aes( x=date, ymin=rt_min, ymax=rt_max, 
                                                               group=key, fill=key )) +
             geom_line( aes(x=pandemic_dates, y=predicted_rt), col="black", alpha= 0.5) +
             geom_hline(yintercept=1, color='black', size= 0.1) + 

             labs( title="", x="", y=expression("Estimated value of" ~ {R[t]}) ) +
                                
             scale_fill_manual( name = "", labels=levels(predicted_values_rt$key),
                                values=c(alpha("seagreen", 0.75), alpha("seagreen", 0.5))) + 
             scale_x_date(date_breaks="months", labels=date_format("%e %b")) +                   # 6Mar21 Juan
   
             theme_pubr() + 
             theme(axis.text.x=element_text(angle=45, hjust=1)) +
             theme(legend.position="right")
  
# plot.deaths.italy.vacunacion.hasta.abril <- 
#   plot_deaths
# plot.deaths.italy.vacunacion.hasta.julio <- 
#   plot_deaths
# plot.deaths.italy.sin.vacunacion.hasta.abril <- 
#   plot_deaths
# plot.deaths.italy.sin.vacunacion.hasta.julio <- 
#   plot_deaths

# plot.rt.italy.vacunacion.hasta.abril <- 
#   plot_rt
# plot.rt.italy.vacunacion.hasta.julio <- 
#   plot_rt
# plot.rt.italy.sin.vacunacion.hasta.abril <- 
#   plot_rt
# plot.rt.italy.sin.vacunacion.hasta.julio <- 
#   plot_rt

# save(list=c("plot.deaths.italy.vacunacion.hasta.abril",
#             "plot.deaths.italy.vacunacion.hasta.julio",
#             "plot.deaths.italy.sin.vacunacion.hasta.julio",
#             "plot.rt.italy.vacunacion.hasta.abril",
#             "plot.rt.italy.vacunacion.hasta.julio",
#             "plot.rt.italy.sin.vacunacion.hasta.julio"),file = "../albertocebolladasolanas/Documents/Master/TFM/PEC3/plots.italia.parte1.Rdata")

 #     Estimated values for 'Effective Rt'

  # Plotting interventions
  predicted_rt_eff       <- colQuantiles(model_rt_eff[,1:N_DATES],probs=0.5)    # 21Mar21
  predicted_rt_eff_95_li <- colQuantiles(model_rt_eff[,1:N_DATES],probs=0.025)
  predicted_rt_eff_95_ui <- colQuantiles(model_rt_eff[,1:N_DATES],probs=0.975)
  predicted_rt_eff_50_li <- colQuantiles(model_rt_eff[,1:N_DATES],probs=0.25 )
  predicted_rt_eff_50_ui <- colQuantiles(model_rt_eff[,1:N_DATES],probs=0.75 )

  predicted_rt_eff_95 <- data.frame(pandemic_dates, predicted_rt_eff_95_li, predicted_rt_eff_95_ui)
  names(predicted_rt_eff_95) <- c("date", "rt_eff_min", "rt_eff_max")
  predicted_rt_eff_95$key <- rep("95%", N_DATES)

  predicted_rt_eff_50 <- data.frame(pandemic_dates, predicted_rt_eff_50_li, predicted_rt_eff_50_ui)
  names(predicted_rt_eff_50) <- c("date", "rt_eff_min", "rt_eff_max")
  predicted_rt_eff_50$key <- rep("50%", N_DATES)

  predicted_values_rt_eff <- rbind(predicted_rt_eff_95, predicted_rt_eff_50)
  predicted_values_rt_eff$key <- factor(predicted_values_rt_eff$key)

  plot_rt_eff <- ggplot(predicted_rt_eff_95) +
             geom_stepribbon( data = predicted_values_rt_eff, aes( x=date, ymin=rt_eff_min, ymax=rt_eff_max,
                                                               group=key, fill=key )) +
             geom_line( aes(x=pandemic_dates, y=predicted_rt_eff), col="black", alpha= 0.5) +
             geom_hline(yintercept=1, color='black', size= 0.1) +

             labs( title="", x="", y=expression("Estimated value of effective" ~ {R[t]}) ) +

             scale_fill_manual( name = "", labels=levels(predicted_values_rt_eff$key),
                                values=c(alpha("seagreen", 0.75), alpha("seagreen", 0.5))) +
             scale_x_date(date_breaks="months", labels=date_format("%e %b")) +                   # 6Mar21 Juan

             theme_pubr() +
             theme(axis.text.x=element_text(angle=45, hjust=1)) +
             theme(legend.position="right")





  # --- Save plots ------------------------------------------------------------------------------ # 
  plot_line1 <- plot_grid(plot_cases, plot_deaths, plot_rt, ncol = 3, rel_widths = c(1, 1, 2))
  plot_line2 <- plot_grid(plot_cum_cases, plot_cum_deaths, plot_rt_eff, ncol = 3, rel_widths = c(1, 1, 2))  
  plot_line3 <- plot_grid(plot_detected, plot_IA14, plot_newhospitalized, plot_hospitalized, ncol = 4, rel_widths = c(1, 1, 1, 1))
  plot_grid <- plot_grid(plot_line1, plot_line2, plot_line3, ncol = 1)
  save_plot(filename=output_graphics_file, plot_grid, base_width=14, base_height = 10.5)    # 21Mar21

  f_log( paste0(" ", Sys.time(), " Graphics created") )
  f_log(" ")
  f_log(" ")
}

