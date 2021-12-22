#**************************************************************************************************
#
#  This script contains some common functions and variable defintions used to obtain data
#  required for the analysis of the development of the COVID-19 in different countries.
#
#  The idea is to separate the steps to read, format, enhance the information, for which
#  specific functions are defined here, from the statistical and predictive analysis of
#  the data.
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
#  Some common definitions
#**************************************************************************************************

#--------------------------------------------------------------------------------------------------
#  Definitions of variables used in the different functions
#--------------------------------------------------------------------------------------------------
SEP_LINE_1 <- paste(replicate(100, "-"), collapse = "")
SEP_LINE_2 <- paste(replicate(100, "+"), collapse = "")
SEP_LINE_3 <- paste(replicate(100, "*"), collapse = "")


#--------------------------------------------------------------------------------------------------
#  Following variables are asumed to be be defined in the environment in which the functions
#  defined here are used. These can be defined for example in the main script using thes
#  functions before calling them.
#--------------------------------------------------------------------------------------------------
#  INPUT_DATA_DIR --> Directory containing the input files
#  LOG_FILE       --> Name of the log file 


#--------------------------------------------------------------------------------------------------
#  Functions to write a text to a log file and additionally to the console
#--------------------------------------------------------------------------------------------------
f_log <- function(text) {

  cat(text, "\n", file=LOG_FILE, append=TRUE)    #  Print to log file
  cat(text, "\n")                                #  Print to console
}

f_log_heading <- function(heading_text) {

   f_log(" ")
   f_log(" ")
   f_log(SEP_LINE_3)
   f_log( paste0(" ", heading_text, " (", Sys.time(), ")") )
   f_log(SEP_LINE_3)
}

#**************************************************************************************************
#  Definitions required to get data from different files
#**************************************************************************************************

#--------------------------------------------------------------------------------------------------
#  Input files
#--------------------------------------------------------------------------------------------------

# COVID-19 daily data for Italy
# Comment lines depend on the analysis that you do want to do
# Analysis until march
input_covid19_data              <- paste0(INPUT_DATA_DIR, '/COVID-19-up-to-date-27-10-21_Italy_run_hasta_abril.csv')
# Analysis until july
# input_covid19_data              <- paste0(INPUT_DATA_DIR, '/COVID-19-up-to-date-27-10-21_Italy_run_hasta_julio.csv')
# Analysis until november
# input_covid19_data              <- paste0(INPUT_DATA_DIR, '/COVID-19-up-to-date-27-10-21_Italy_run_hastaNoviembre.csv')

#  Non pharmacological measures (interventios)
#  ("https://www.ecdc.europa.eu/en/publications-data/download-data-response-measures-covid-19")
# input_covid19_response_measures <- paste0(INPUT_DATA_DIR, '/response_graphs_data_Spain_Clean_run4.csv')
# input_covid19_response_measures <- paste0(INPUT_DATA_DIR, '/response_graphs_data_Spain_Clean_run_hasta_enero.csv')
# input_covid19_response_measures <- paste0(INPUT_DATA_DIR, '/response_graphs_data_Italy_Clean_run2.csv')

# Comment lines depend on the analysis that you do want to do
# Analysis until march
input_covid19_response_measures <- paste0(INPUT_DATA_DIR, '/response_graphs_data_Italy_Clean_run_hasta_marzo.csv')
# Analysis until july
#input_covid19_response_measures <- paste0(INPUT_DATA_DIR, '/response_graphs_data_Italy_Clean_run_hasta_marzo.csv')
# Analysis until november
#input_covid19_response_measures <- paste0(INPUT_DATA_DIR, '/response_graphs_data_Italy_Clean_run_hasta_noviembre.csv')

#  Case fatality rate
input_case_fatality_rate        <- paste0(INPUT_DATA_DIR, '/weighted_fatality.csv')

#  Serial interval
input_serial_interval           <- paste0(INPUT_DATA_DIR, '/serial_interval.csv')

#  Infection to death
input_infection_to_death        <- paste0(INPUT_DATA_DIR, '/infection_to_death.csv')

#  Hospitalization fatality rate
input_case_hospitalization_rate <- paste0(INPUT_DATA_DIR, '/weighted_hospitalization.csv')

#  Vaccination data
# input_vaccination <- paste0(INPUT_DATA_DIR, '/datos.vacunacion.spain.work.csv')
input_vaccination <- paste0(INPUT_DATA_DIR, '/datos.vacunacion.italia.work.csv')

#  Vaccination to inmunization first dose
input_vtimm_first_dose    <- paste0(INPUT_DATA_DIR, '/vtimm.first.dose.csv')

#  Vaccination to inmunization first dose
input_vtimm_second_dose    <- paste0(INPUT_DATA_DIR, '/vtimm.second.dose.csv')

#  IFR
input_ifr  <-  paste0(INPUT_DATA_DIR, '/IFR_by_age.csv')

#--------------------------------------------------------------------------------------------------
#  Function f_search_country_in_list
#  ---------------------------------
#
#  Function to search the value of a ocontry in a list
#
#
#   Arguments
#   ---------
#     select_country     Country to be analyzed
#     country_list       List of countries
#
#
#   Return value
#   ------------
#     List of indices of the list 'country_list' in which the 'select_country' was found.
#     The search if done via this function to allow following criteria:
#      - Ignore case           'Spain' is considered the same as 'SPAIN'
#      - Space and underscore  'United Kingdom' is considered the same as 'United_Kingdom'
#      - Trim                  Remove trailing and leading spaces
#
#     As a summary, for the purpose of this search, following values are considered identical:
#      - "United Kingdom"
#      - "United_Kindgom  "
#      - "   United_Kindgom  "
#      - "___United_Kindgom_"
#      - "UNITED KINGDOM"
#      - ...
#--------------------------------------------------------------------------------------------------
f_get_trimmed_country_str <- function(select_country) {

  select_country <- gsub("_", " ", select_country)
  select_country <- trim(tolower(select_country))
  select_country <- gsub(" ", "_", select_country)
  
  return(select_country)
}


f_search_country_in_list <- function( select_country, country_list) {
    
    select_country <- f_get_trimmed_country_str(select_country)
    country_list   <- lapply(country_list, f_get_trimmed_country_str)
    
    which(country_list==select_country)
}


#--------------------------------------------------------------------------------------------------
#
#  Function f_get_covid_data_for_country
#  -------------------------------------
#
#  Read data with information about the pandemic, select data for the country being
#  processed and determine the start date of the pandemic.
#  Only the data from this date on will be considered for the analysis.
#
#
#   Arguments
#   ---------
#     select_country     Country to be analyzed
#
#     Parameters to determine the start date of the pandemic:
#      - n_deaths_0      Number of deaths used as limit to consider that a pandemic started before
#      - n_days_before   Days before the date when 'n_deaths_0' was reached
#                        used to define the start of the pandemic
#
#      Asuming that the default values are used (see function definition), the start of the
#      pandemic is defined as the date 30 days before the number total of deaths reached the
#      value 10.
#
#                                         # PENDING: if n_death_0=0 start on the first day of the data series   18Apr21
#   Return value
#   ------------
#     Dataframe 'covid19_country_data' containing the relevant data for the selected country
#
#--------------------------------------------------------------------------------------------------
f_get_covid_data_for_country <- function( select_country, n_deaths_0=10, n_days_before=30) {

  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" COVID-19 data")
  f_log(SEP_LINE_2)

  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start loading COVID-19 data") )
  f_log(SEP_LINE_1)
  f_log( paste0(" Read data from file '", basename(input_covid19_data),
                "' for '", select_country, "'") )


  # --- Read data and select and extract relevant information for the country being processed --- #

  #     Read file containing daily data for all countries
  covid19_all_data          <- read.csv(input_covid19_data)


  #     Select data for country to be processed
  country_indices           <- f_search_country_in_list( select_country,
                                                         covid19_all_data$Countries.and.territories )
  covid19_country_data      <- covid19_all_data[country_indices, ]
 
  if (nrow(covid19_country_data) == 0)  {  return(NULL)  }
  

  #     Sort data for selected country by datae
  
  covid19_country_data$date <- as.Date(covid19_country_data$DateRep,format='%d/%m/%Y')
  covid19_country_data$t    <- decimal_date(covid19_country_data$date)
  # covid19_country_data$date <- format(  covid19_country_data$date,"%d/%m/%Y") ### Add Alberto
  covid19_country_data      <- covid19_country_data[order(covid19_country_data$t),]


  #     Find index for date with first deaths
  index_first_death = which(covid19_country_data$Deaths>0)[1]


  #     Find index for date in which the number of deaths defined in 'n_deaths_0' was reached
  index_first_n_deaths = which(cumsum(covid19_country_data$Deaths)>=n_deaths_0)[1]


  #     Pandemic start is defined as 30 days before the number of deaths 'n_deaths_0' is reached
  index_epidemic_start = index_first_n_deaths - n_days_before
  if (index_epidemic_start <=0) {  index_epidemic_start <- 1  }


  #     Select data starting in the index that defines the start of the pandemic
  covid19_country_data_orig <- covid19_country_data
  covid19_country_data <- covid19_country_data[index_epidemic_start:nrow(covid19_country_data),]


  #     Number of days with data since the start of the pandemic
  N = nrow(covid19_country_data)


  # --- Summary of data read -------------------------------------------------------------------- #
  f_get_date_for_index <- function(data, index) { return( data$date[index] ) }

  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Loaded data from file '", basename(input_covid19_data),
                "' for '", select_country, "'") )
  f_log(SEP_LINE_1)

  f_log( paste0(" First death is on day ", index_first_death,
                " (", f_get_date_for_index(covid19_country_data_orig, index_first_death), ")") )

  f_log( paste0(" First ", n_deaths_0, " deaths reached on day ", index_first_n_deaths,
                " (", f_get_date_for_index(covid19_country_data_orig, index_first_n_deaths), ")") )

  f_log( paste0(" Pandermic start asumed 30 days before first ", n_deaths_0, " deaths,",
                " on day ", index_epidemic_start,
                " (", f_get_date_for_index(covid19_country_data_orig, index_epidemic_start), ")") )

  f_log( paste0(" Number of days with data since the start of the pandemic: ", N,
                " (from ", f_get_date_for_index(covid19_country_data, 1),
                " to ", f_get_date_for_index(covid19_country_data, N), ")") )

  f_log(" ")
  f_log(" ")


  # --- Return data ----------------------------------------------------------------------------- #
  return(covid19_country_data)
}


#--------------------------------------------------------------------------------------------------
#
#  Function f_get_measures_for_country
#  -----------------------------------
#
#  Read information about the response measures applied to a country and return the information.
#  The information is enhanced with a column for each date being analyzing, indicating which
#  measures were active in each specific date.
#
#
#   Arguments
#   ---------
#     select_country     Country to be analyzed
#     dates_for_country  List of dates with data available to be analyze for the selected country
#
#
#   Return value
#   ------------
#     Dataframe 'response_measures_country' containing the relevant information about the measures
#
#--------------------------------------------------------------------------------------------------
f_get_measures_for_country <- function( select_country, dates_for_country ) {

  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Response measures (interventions)")
  f_log(SEP_LINE_2)

  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start loading measures data") )
  f_log(SEP_LINE_1)
  f_log( paste0(" Read data from file '", basename(input_covid19_response_measures),
                "' for '", select_country, "'") )


  # --- Read data and select and extract relevant information for the country being processed --- #

  #     Read file containing information about the measures
  response_measures_in <- read.csv(input_covid19_response_measures, stringsAsFactors=T)


  #     Select data for country to be processed
  country_indices           <- f_search_country_in_list( select_country,
                                                         response_measures_in$Country )
  response_measures_country <- response_measures_in[country_indices, ]
  
  
  #     Create dummy measure if none defined - this is required to simplify the implementation, 
  #     to avoid having to code a speficific handling for countries woithout measures.
  #     Even if this is not realistic, a country without measures may result is incosistent data
  #     files are used. In this case, the dummy measure is created with a date in the future to
  #     ensure that the return dataframe *response_measures_country*, but also that the dummy
  #     measure has no impact on the analysis.
  dummy_measure <- NULL
  if ( nrow(response_measures_country) == 0) {
    dummy_measure             <- response_measures_in[1,]
    dummy_measure$Country     <- select_country  
    dummy_measure$date_start  <- "9999-12-31"
    dummy_measure$date_end    <- "9999-12-31"
    # dummy_measure$date_start  <- "31/12/9999"
    # dummy_measure$date_end    <- "31/12/9999"
    response_measures_country <- dummy_measure 
  }
  
  
  #     Add a column for each date to be analyzed and include for each measure if this is
  #     active in the corresponding date
  for(i in 1:nrow(response_measures_country)) {

    measure_i_name <- response_measures_country$Response_measure[i]
    measure_i_from <- response_measures_country$date_start[i]
    measure_i_to   <- response_measures_country$date_end[i]

    if (is.na(measure_i_to)){                                 # 24Apr21
      measure_i_to <- "9999-12-31"                   # 24Apr21
      # measure_i_to <- "31/12/9999"                   # 24Apr21 # add Alberto
      levels(response_measures_country$date_end) <- c(levels(response_measures_country$date_end), measure_i_to )  # 24Apr21 this is needed because this variable is a factor
      response_measures_country$date_end[i] <- measure_i_to    # 24Apr21
    }




    for (j in 1:length(dates_for_country)) {

      response_measures_country[ i, as.character(dates_for_country[j]) ] <- 0

      if( dates_for_country[j] >= as.Date(measure_i_from) &&
          dates_for_country[j] <= as.Date(measure_i_to) ) {
          
        response_measures_country[ i, as.character(dates_for_country[j]) ] <- 1
      }
    }
  }

 
  # PENDING: add here instead of 1, a given factor to use a given fixed value?

  # merge active dates for same measures    20Mar21

    response_measures_country_unique = response_measures_country[1,]   # to initialize response_measures_country_unique list
    for(i in 1:length(unique(response_measures_country$Response_measure))){
        s_measure = unique(response_measures_country$Response_measure)[i]
        k = which(response_measures_country$Response_measure==s_measure)[1]
        response_measures_country_unique[i,] = response_measures_country[k,]
        if (length(which(response_measures_country$Response_measure==s_measure))>1){
		for(j in 2:length(which(response_measures_country$Response_measure==s_measure))){
			k=which(response_measures_country$Response_measure==s_measure)[j]
			response_measures_country_unique[i,as.character(dates_for_country)] = response_measures_country_unique[i,as.character(dates_for_country)] +
			                      	response_measures_country[k,as.character(dates_for_country)] 
			    # PENDING: now if periods with same measures are overlapping, they count only once in covid19_base.stan... do it better?
			if(as.Date(response_measures_country$date_start[k]) < as.Date(response_measures_country_unique$date_start[i])){
			       levels(response_measures_country_unique$date_start) <- c(levels(response_measures_country_unique$date_start), response_measures_country$date_start[k])	
                               response_measures_country_unique$date_start[i] <- response_measures_country$date_start[k]
		        }
			if(as.Date(response_measures_country$date_end[k]) > as.Date(response_measures_country_unique$date_end[i])){
			       levels(response_measures_country$date_end) <- c(levels(response_measures_country$date_end), response_measures_country$date_end[k] )	
                               response_measures_country_unique$date_end[i] <- response_measures_country$date_end[k]
                        }
	        }
	}
    }
    response_measures_country = response_measures_country_unique   # 20Mar21



  #     Sort data by start date of application of the measures
  response_measures_country <- response_measures_country[order(response_measures_country$date_start),]


  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Loaded data from file '", basename(input_covid19_response_measures),
                "' for '", select_country, "'") )
  f_log(SEP_LINE_1)

  n_measure_periods <- nrow(response_measures_country)
  n_unique_measures <- length(unique(response_measures_country$Response_measure))

  if( is.null(dummy_measure) ) { 
    f_log( paste0(" Number of periods with defined measures: ", n_measure_periods) )
    f_log( paste0(" Number of different measures:            ", n_unique_measures) )
  } else {
    f_log( paste0(" No measures defined for '", select_country, "'") )
  }

  f_log(" ")
  f_log(" ")


  # --- Return data ----------------------------------------------------------------------------- #
  return(response_measures_country)
}



#--------------------------------------------------------------------------------------------------
#
#  Function f_get_numeric_info_from_measures
#  -----------------------------------------
#
#  Get numeric information from an object obteined with the function 'f_get_measures_for_country'
#  (see definition above). Objects generated with that function contain information about
#  non pharmacological measures (measure name, start and end dates, and the information for each
#  date about the active measures). This function extracts this last information, returning a
#  dataframe containing one row for each measure period and a colume for each date object of
#  analysis.
#
#
#   Arguments
#   ---------
#     response_measures_country     Dataframe obteined with the function 'f_get_measures_for_country'
#
#
#   Return value
#   ------------
#     Dataframe 'response_measures_country_dates', in which only the columns containing information
#     about each date are defined
#
#--------------------------------------------------------------------------------------------------
f_get_numeric_info_from_measures <- function( response_measures_country ) {

  #  Return NULL if 'response_measures_country' is NULL or contains no data
  if ( is.null(response_measures_country) || nrow(response_measures_country) == 0 ) {
    return(NULL)
  }
  
  #  Return relevant information
  return( response_measures_country[,-c(1:4)] )
}



#--------------------------------------------------------------------------------------------------
#
#  Function f_get_case_fatality_rate_for_country
#  ---------------------------------------------
#
#  Get case fatality rate (CFR)
#
#
#   Arguments
#   ---------
#     select_country     Country to be analyzed
#
#
#   Return value
#   ------------
#     Value 'case_fatality_rate_country' containing the case fatality rate for the selected country
#
#--------------------------------------------------------------------------------------------------
f_get_case_fatality_rate_for_country <- function( select_country ) {

  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Case fatality rate (CFR)")
  f_log(SEP_LINE_2)

  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start loading case fatality rate (CFR) data") )
  f_log(SEP_LINE_1)
  f_log( paste0(" Read data from file '", basename(input_case_fatality_rate),
                "' for '", select_country, "'") )


  # --- Read data and create list with the information for the selected country ----------------- #

  #     Read file containing case fatality rate (CFR) data for all countries
  case_fatality_rate_in <- read.csv(input_case_fatality_rate)


  #     Assign adequate name to second column
  colnames(case_fatality_rate_in)[2] <- "country"


  #     Select data for country to be processed
  country_indices            <- f_search_country_in_list( select_country,
                                                          case_fatality_rate_in$country ) 
  case_fatality_rate_country <- case_fatality_rate_in$weighted_fatality[country_indices]
  
  
  info_str <- ""
  if (length(case_fatality_rate_country) == 0) {

    case_fatality_rate_country <- mean(case_fatality_rate_in$weighted_fatality)
    info_str <- paste0(" No data found for '", select_country, "'",
                       " / mean value of all countries will be used: ", case_fatality_rate_country)

  } else {

    info_str <- paste0(" Case fatality rate for '", select_country, "': ", case_fatality_rate_country)
  }




  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Loaded data from file '", basename(input_case_fatality_rate),
                "' for '", select_country, "'") )
  f_log(SEP_LINE_1)
  f_log(info_str)

  f_log(" ")
  f_log(" ")


  # --- Return data ----------------------------------------------------------------------------- #
  return(case_fatality_rate_country)
}


#--------------------------------------------------------------------------------------------------
#
#  Function f_get_vaccination_for_country ## ADDED ALBERTO 17/11/2021
#--------------------------------------------------------------------------------------------------
#
#  Get case fatality rate (CFR)
#
#
#   Arguments
#   ---------
#     select_country     Country to be analyzed
#
#
#   Return value
#   ------------
#     Data frame 'vaccination_country' containing the vaccination data from selected country
#
#--------------------------------------------------------------------------------------------------
f_get_vaccination_for_country <- function( select_country ) {
  
  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Vaccination")
  f_log(SEP_LINE_2)
  
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start loading vaccination data") )
  f_log(SEP_LINE_1)
  f_log( paste0(" Read data from file '", basename(input_vaccination),
                "' for '", select_country, "'") )
  
  
  # --- Read data and create list with the information for the selected country ----------------- #
  
  #     Read file containing vaccination data for all countries
  vaccination <- read.csv(input_vaccination)

  #     Select data for country to be processed
  country_indices            <- f_search_country_in_list( select_country,
                                                          vaccination$Country ) 
  vaccination_country  <- vaccination[country_indices,]
  
  
  
  info_str <- ""
  if (nrow(vaccination_country) == 0) {
    
    vaccination_return <- NULL
    info_str <- paste0(" No data found for '", select_country, "'")
                       
    
  } else {
    
    info_str <- paste0(" Data found for '", select_country, "'")
  }
  
  
  
  
  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Loaded data from file '", basename(input_vaccination),
                "' for '", select_country, "'") )
  f_log(SEP_LINE_1)
  f_log(info_str)
  
  f_log(" ")
  f_log(" ")
  
  
  # --- Return data ----------------------------------------------------------------------------- #
  return(vaccination_country)
}

#--------------------------------------------------------------------------------------------------
#
#  Function f_get_first_dose ## ADDED ALBERTO 17/11/2021
#  ---------------------------------------------
#
#  Get first dose vaccination
#
#
#   Arguments
#   ---------
#     data_vaccination     Complete data from vaccination
#
#
#   Return value
#   ------------
#     Columns that  containing the data from first dose
#
#--------------------------------------------------------------------------------------------------
f_get_first_dose <- function( data_vaccination,select_country,pandemic_dates ) {
  
  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Vaccination First Dose")
  f_log(SEP_LINE_2)
  
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start extract first dose vaccination data") )
  f_log(SEP_LINE_1)

  
  

  info_str <- ""
  
  vaccination_dates <- 
    seq(from=(min(pandemic_dates) - 28),to=max(pandemic_dates),by=1)
  if (nrow(data_vaccination) == 0) {
    
    data_first_dose <- NULL
    info_str <- paste0(" No data found for '", "'")
    
    
  } else {
    data_first_dose <- 
      data_vaccination %>% 
      mutate(date=as.Date(paste0(date))) %>% 
      filter(date%in%vaccination_dates) %>% 
      select(Country,
             date,
             contains("first_dose"),
             age_group)
    dates.not.vaccination.found <- 
      vaccination_dates[!vaccination_dates%in%as.Date(data_vaccination$date)]
    levels.age <- 
      levels(factor(data_vaccination$age_group))
    vaccination.dates.not.vaccination.found <- 
                       crossing(date=dates.not.vaccination.found,
                                age_group=levels.age) %>% 
      mutate(Country=select_country,
             n.first_dose=0,
             first_dose=0,
             n.first_dose_perc=0,
             first_dose_perc=0)
    data_first_dose <-
      rbind.data.frame(data_first_dose,vaccination.dates.not.vaccination.found) %>% 
      arrange(date)
    
    if(nrow(data_first_dose)!=0){
    info_str <- paste0(" Data found for '",  "'")
    } 
  }
  
  
  
  
  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Extract first dose from vaccination'", 
                "' for '", select_country, "'") )
  f_log(SEP_LINE_1)
  f_log(info_str)
  
  f_log(" ")
  f_log(" ")
  
  
  # --- Return data ----------------------------------------------------------------------------- #
  return(data_first_dose)
}


#--------------------------------------------------------------------------------------------------
#
#  Function f_get_second_dose ## ADDED ALBERTO 17/11/2021
#  ---------------------------------------------
#
#  Get second dose vaccination
#
#
#   Arguments
#   ---------
#     data_vaccination     Complete data from vaccination
#
#
#   Return value
#   ------------
#     Columns that  containing the data from first dose
#
#--------------------------------------------------------------------------------------------------
f_get_second_dose <- function( data_vaccination,select_country,pandemic_dates ) {
  
  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Vaccination Second Dose")
  f_log(SEP_LINE_2)
  
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start extract second dose vaccination data") )
  f_log(SEP_LINE_1)
  
  
  
  
  info_str <- ""
  vaccination_dates <- 
    seq(from=(min(pandemic_dates) - 28),to=max(pandemic_dates),by=1)
  if (nrow(data_vaccination) == 0) {
    
    data_second_dose <- NULL
    info_str <- paste0(" No data found for '", "'")
    
    
  } else {
    data_second_dose <- 
      data_vaccination %>%  
      mutate(date=as.Date(paste0(date))) %>% 
      filter(date%in%vaccination_dates) %>% 
      select(Country,
             date,
             contains("second_dose"),
             age_group)
    
    dates.not.vaccination.found <- 
      vaccination_dates[!vaccination_dates%in%as.Date(data_vaccination$date)]
    levels.age <- 
      levels(factor(data_vaccination$age_group))
    vaccination.dates.not.vaccination.found <- 
      crossing(date=dates.not.vaccination.found,
               age_group=levels.age) %>% 
      mutate(Country=select_country,
             n.second_dose=0,
             second_dose=0,
             n.second_dose_perc=0,
             second_dose_perc=0)
    data_second_dose <-
      rbind.data.frame(data_second_dose,vaccination.dates.not.vaccination.found) %>% 
      arrange(date)
    
    if(nrow(data_second_dose)!=0){
      info_str <- paste0(" Data found for '",  "'")
    } else{
      info_str <- paste0("Not Data found for those dates'",  "'")
      }
  }
  
  
  
  
  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Extract second dose from vaccination'", 
                "' for '", select_country, "'") )
  f_log(SEP_LINE_1)
  f_log(info_str)
  
  f_log(" ")
  f_log(" ")
  
  
  # --- Return data ----------------------------------------------------------------------------- #
  return(data_second_dose)
}

# ADDED 17 november Alberto 
#--------------------------------------------------------------------------------------------------
#
#  Function f_get_population_by_age
#  ---------------------------------------------
#
#  Get population 
#
#
#   Arguments
#   ---------
#     data_vaccination     Data from Vaccination
#
#
#   Return value
#   ------------
#     Values 'population' containing the population size for the selected country (raw value) and group age 
#
#--------------------------------------------------------------------------------------------------
f_get_population_by_age <- function( data_vaccination,select_country,pandemic_dates ) {
  
  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Population")
  f_log(SEP_LINE_2)
  
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start loading population data") )
  f_log(SEP_LINE_1)
  f_log( paste0(" Load data from vaccination '"))
  
 
  #     Select data for country to be processed
  info_str <- ""
  vaccination_dates <- 
    seq(from=(min(pandemic_dates) - 28),to=max(pandemic_dates),by=1)
  if (nrow(data_vaccination) == 0) {
    
    data_population_age <- NULL
    info_str <- paste0(" No data found for '", "'")
    
    
  } else {
    data_population_age <- 
      data_vaccination %>%  
      mutate(date=as.Date(paste0(date))) %>% 
      filter(date%in%vaccination_dates) %>% 
      select(Country,
             date,
             age_group,
             Poblacion)
    dates.not.vaccination.found <- 
      vaccination_dates[!vaccination_dates%in%as.Date(data_vaccination$date)]
    levels.age <- 
      levels(factor(data_vaccination$age_group))
    population_value <- 
    data_population_age %>% 
      select(age_group,Poblacion) %>% 
      distinct()
    
    population.dates.not.vaccination.found <- 
      crossing(date=dates.not.vaccination.found,
               population_value) %>% 
      mutate(Country=select_country)
    
    
    data_population_age <-
      rbind.data.frame(data_population_age,population.dates.not.vaccination.found) %>% 
      arrange(date)
    info_str <- paste0(" Data found for '",  "'")
  }
  
  
  
  
  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Extract Poblacion and Group Age'", 
                "' for '", select_country, "'") )
  f_log(SEP_LINE_1)
  f_log(info_str)
  
  f_log(" ")
  f_log(" ")
  
  
  # --- Return data ----------------------------------------------------------------------------- #
  return(data_population_age)
}



# ADDED 17 november Alberto 
#--------------------------------------------------------------------------------------------------
#
#  Function f_get_population_inmunizated
#  ---------------------------------------------
#
#  Get population 
#
#
#   Arguments
#   ---------
#     select_country
#     pandemic_dates
#     data_vaccination_first_dose     Data from Vaccination
#     curve_inmunization
#     population_by_age
#
#   Return value
#   ------------
#     List where each element is a data frame by day where record the number of population inmunizated by age group 
#
#--------------------------------------------------------------------------------------------------
f_get_population_immunized <- function(select_country,pandemic_dates,vaccination_country_first_dose,VTIMM_FIRST,population_by_age ) {
  
  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Population")
  f_log(SEP_LINE_2)
  
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start loading population data") )
  f_log(SEP_LINE_1)
  f_log( paste0(" Load data from vaccination '"))
  
  
  #     Select data for country to be processed
  
  population_first_day <-
    population_by_age %>% 
    filter(date==min(date))
  
  population_inmunizated_by_day_by_group_age <-
    lapply(seq_along(pandemic_dates),FUN=function(index.date){
      dates.to.get.vaccination <- 
        seq(from=min(vaccination_country_first_dose$date),
            to=pandemic_dates[index.date]-1,
            by=1)
      if(length(dates.to.get.vaccination) > 
         length(VTIMM_FIRST)){
        n.to.add <- 
          length(dates.to.get.vaccination)- length(VTIMM_FIRST)
        VTIMM_FIRST <- 
          c(VTIMM_FIRST,rep(0.95,n.to.add))
      }
      
      coef_v.t.imm <- 
        VTIMM_FIRST %>% 
        head(length(dates.to.get.vaccination)) %>% 
        rev()
      
        
      vaccination_country_first_dose %>% 
        filter(date%in%dates.to.get.vaccination) %>% 
        left_join(data.frame(date=dates.to.get.vaccination,
                             coef_v.t.imm),by = "date") %>% 
        mutate(pop.inmunizated=n.first_dose*coef_v.t.imm) %>% 
        group_by(age_group) %>% 
        summarise(poblacion.inmunizated=floor(sum(pop.inmunizated)))
        
    })
    
  names(population_inmunizated_by_day_by_group_age) <-
    pandemic_dates
    
  info_str <- paste0(" Data found for '",  "'")
  
  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Extract Poblacion Inmunizated by day'", 
                "' for '", select_country, "'") )
  f_log(SEP_LINE_1)
  f_log(info_str)
  
  f_log(" ")
  f_log(" ")
  
  
  # --- Return data ----------------------------------------------------------------------------- #
  return(population_inmunizated_by_day_by_group_age)
}


# ADDED 17 november Alberto 
#--------------------------------------------------------------------------------------------------
#
#  Function f_get_CFR_by_day
#  ---------------------------------------------
#
#  Calculate CFR each day. We update population in risk with data of population inmunizated
#
#
#   Arguments
#   ---------
#     country_selected
#     pandemic_dates
#     population.inmunizated.by.day
#     population_by_age
#     IFR
#
#   Return value
#   ------------
#     List where each element is a data frame by day where record the number of population inmunizated by age group 
#
#--------------------------------------------------------------------------------------------------

f_get_CFR_by_day <- function(country_selected,pandemic_dates,population.inmunizated.by.day,population_by_age,IFR){
# --- Start of the processing ----------------------------------------------------------------- #
f_log(" ")
f_log(" ")
f_log(SEP_LINE_2)
f_log(" Population")
f_log(SEP_LINE_2)

f_log(SEP_LINE_1)
f_log( paste0(" ", Sys.time(), " - Start loading population data") )
f_log(SEP_LINE_1)
f_log( paste0(" Load data from vaccination '"))


population.first.day <-
  population_by_age %>% 
  filter(date==min(pandemic_dates))
population.inmunizated.by.day <-
  lapply(population.inmunizated.by.day,FUN=function(day.pop.inmunizated){
    day.pop.inmunizated <-
      day.pop.inmunizated %>% 
      left_join(population.first.day,by="age_group") 
    day.pop.inmunizated <-
      day.pop.inmunizated %>% 
      mutate(poblacion.in.riesgo=Poblacion-poblacion.inmunizated)
    day.pop.inmunizated
  })


## OJO esto lo modifico para adaptar los grupos de edad de ITALIA
IFR <- paste0(IFR)
IFR <-
  c(IFR[-1],IFR[(length(IFR))])

## OJO esto lo modifico para adaptar los grupos de edad de ESPAÑA
# IFR <- as.numeric(paste0(IFR))
# IFR <-
#    c(mean(IFR[1:5]),IFR[6],IFR[7],IFR[8],IFR[9])

age.groups.levels <- 
  levels(factor(population_by_age[["age_group"]]))
age.groups.levels <- 
  age.groups.levels[age.groups.levels!="Total"]
# esto para españa
# age.groups.levels <-
#   c("Menos de 49","50-59","60-69","70-79","Mas de 80")

names(IFR) <- age.groups.levels
IFR.df <- 
  data.frame(age_group=names(IFR),
             IFR=IFR)
day.pop.inmunizated <- population.inmunizated.by.day[[1]]
CFR.by.day <- 
   lapply(population.inmunizated.by.day,FUN=function(day.pop.inmunizated){
     day.pop.inmunizated <- 
     day.pop.inmunizated %>% 
       left_join(IFR.df,by="age_group") %>% 
       mutate(IFR=as.numeric(paste0(IFR)))
     day.pop.inmunizated %>% 
       filter(age_group!="Total") %>% 
       summarise(IFR_return=sum(poblacion.in.riesgo*IFR)/sum(poblacion.in.riesgo))
})

return(CFR.by.day)
}
# ADDED 7Mar21 Juan
#--------------------------------------------------------------------------------------------------
#
#  Function f_get_population_for_country
#  ---------------------------------------------
#
#  Get population 
#
#
#   Arguments
#   ---------
#     select_country     Country to be analyzed
#
#
#   Return value
#   ------------
#     Value 'population_K_country' containing the population size for the selected country (in thousands units)
#
#--------------------------------------------------------------------------------------------------
f_get_population_for_country <- function( select_country ) {

  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Population")
  f_log(SEP_LINE_2)

  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start loading population data") )
  f_log(SEP_LINE_1)
  f_log( paste0(" Read data from file '", basename(input_case_fatality_rate),
                "' for '", select_country, "'") )

  # --- Read data and create list with the information for the selected country ----------------- #

  #     Read file containing population data for all countries     # this value (in thousands units) is located in the same file as CFR    7Mar21 Juan
  case_fatality_rate_in <- read.csv(input_case_fatality_rate)


  #     Assign adequate name to second column
  colnames(case_fatality_rate_in)[2] <- "country"


  #     Select data for country to be processed
  country_indices            <- f_search_country_in_list( select_country,
                                                          case_fatality_rate_in$country )
  population_K_country <- case_fatality_rate_in$population[country_indices]


  info_str <- ""
  if (length(population_K_country) == 0) {

    population_K_country <- mean(case_fatality_rate_in$population)
    info_str <- paste0(" No data found for '", select_country, "'",
                       " / mean value of all countries will be used: ", population_K_country)

  } else {

    info_str <- paste0(" Population size for '", select_country, "': ", population_K_country)
  }

  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Loaded data from file '", basename(input_case_fatality_rate),
                "' for '", select_country, "'") )
  f_log(SEP_LINE_1)
  f_log(info_str)

  f_log(" ")
  f_log(" ")


  # --- Return data ----------------------------------------------------------------------------- #
  return(population_K_country)
}
# ADDED 7Mar21 Juan




# 27Mar21 Juan
#--------------------------------------------------------------------------------------------------
#
#  Function f_get_case_hospitalization_rate
#  ---------------------------------------------
#
#  Get case hospitalization rate for selected country 
#
#
#   Arguments
#   ---------
#     select_country     Country to be analyzed
#
#
#   Return value
#   ------------
#     Value 'case_hospitalization_rate' containing the hospitalization rate for infected cases
#
#--------------------------------------------------------------------------------------------------
f_get_case_hospitalization_rate <- function( select_country ) {

  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Case Hospitalization Rate")
  f_log(SEP_LINE_2)

  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start loading case hospitalization rate data") )
  f_log(SEP_LINE_1)
  f_log( paste0(" Read data from file '", basename(input_case_hospitalization_rate),
                "' for '", select_country, "'") )

  # --- Read data and create list with the information for the selected country ----------------- #

  #     Read file containing case hospitalization rate data for all countries     
  case_hospitalization_rate_in <- read.csv(input_case_hospitalization_rate)


  #     Assign adequate name to second column
  colnames(case_hospitalization_rate_in)[2] <- "country"


  #     Select data for country to be processed
  country_indices            <- f_search_country_in_list( select_country,
                                                          case_hospitalization_rate_in$country )
  case_hospitalization_rate <- case_hospitalization_rate_in$weighted_hospitalization[country_indices]

  info_str <- ""
  if (length(case_hospitalization_rate) == 0) {

    case_hospitalization_rate <- mean(case_hospitalization_rate_in$weighted_hospitalization)
    info_str <- paste0(" No data found for '", select_country, "'",
                       " / mean value of all countries will be used: ", case_hospitalization_rate)

  } else {

    info_str <- paste0(" Hospitalization rate for '", select_country, "': ", case_hospitalization_rate)

  }

  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Loaded data from file '", basename(input_case_hospitalization_rate),
                "' for '", select_country, "'") )
  f_log(SEP_LINE_1)
  f_log(info_str)

  f_log(" ")
  f_log(" ")


  # --- Return data ----------------------------------------------------------------------------- #
  return(case_hospitalization_rate)
}

# 27Mar21 Juan






#--------------------------------------------------------------------------------------------------
#
#  Function f_get_serial_interval
#  ------------------------------
#
#  Get serial interval (SI)
#
#
#   Arguments
#   ---------
#     requested_interval_length   Number of values required of the serial interval.
#                                 If this number is higher than the number of records available
#                                 in the input file, the missing records will be set to 0.

#
#   Return value
#   ------------
#     List 'serial_interval'
#
#--------------------------------------------------------------------------------------------------
f_get_serial_interval <- function( requested_interval_length ) {

  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Serial interval (SI)")
  f_log(SEP_LINE_2)

  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start loading serial interval (SI) data") )
  f_log(SEP_LINE_1)
  f_log( paste0(" Read data from file '", basename(input_serial_interval), "'") )


  # --- Read data and create list with the expected length -------------------------------------- #

  #     Read file containing information about the 'serial interval'
  serial_interval_in <- read.csv(input_serial_interval)


  #     Select data from column containing the relevant information
  serial_interval <- serial_interval_in$fit


  #     Adjust length of 'serial_interval' adding data to or removing data from the read records
  n_serial_interval_in = length(serial_interval)
  if (n_serial_interval_in < requested_interval_length) {
    n_missing_data <- requested_interval_length - n_serial_interval_in
    serial_interval[(n_serial_interval_in+1):requested_interval_length] <- rep(0.,n_missing_data)
  }

  serial_interval <- serial_interval[1:requested_interval_length]


  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Loaded data from file '", basename(input_serial_interval), "'") )
  f_log(SEP_LINE_1)

  f_log( paste0(" Number of records in file:        ", nrow(serial_interval_in)) )
  f_log( paste0(" Number of records in return list: ", length(serial_interval)) )

  f_log(" ")
  f_log(" ")


  # --- Return data ----------------------------------------------------------------------------- #
  return(serial_interval)
}


# Added Alberto November 
#--------------------------------------------------------------------------------------------------
#
#  Function f_get_vtimm_first_dose
#  ------------------------------
#
#  Get Vaccination to Immunization
#
#
#   Arguments
#   ---------
#     requested_interval_length   Number of values required of the serial interval.
#                                 If this number is higher than the number of records available
#                                 in the input file, the missing records will be set to 0.95

#
#   Return value
#   ------------
#     List 'vtimm_first_dose'
#
#--------------------------------------------------------------------------------------------------
f_get_vtimm_first_dose <- function( requested_interval_length ) {
  
  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Vaccination to immunization")
  f_log(SEP_LINE_2)
  
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start loading Vaccination Immunization data") )
  f_log(SEP_LINE_1)
  f_log( paste0(" Read data from file '", basename(input_vtimm_first_dose), "'") )
  
  
  # --- Read data and create list with the expected length -------------------------------------- #
  
  #     Read file containing information about the 'serial interval'
  vtimm_first_dose_in <- read.csv(input_vtimm_first_dose)
  
  
  #     Select data from column containing the relevant information
  vtimm_first_dose <- vtimm_first_dose_in$effectivity
  
  
  #     Adjust length of 'serial_interval' adding data to or removing data from the read records
  n_vtimm_first_dose_in = length(vtimm_first_dose)
  if (n_vtimm_first_dose_in < requested_interval_length) {
    n_missing_data <- requested_interval_length - n_vtimm_first_dose_in
    vtimm_first_dose[(n_vtimm_first_dose_in+1):requested_interval_length] <- rep(0.95,n_missing_data)
  }
  
  vtimm_first_dose <- vtimm_first_dose[1:requested_interval_length]
  
  
  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Loaded data from file '", basename(input_vtimm_first_dose), "'") )
  f_log(SEP_LINE_1)
  
  f_log( paste0(" Number of records in file:        ", nrow(vtimm_first_dose_in)) )
  f_log( paste0(" Number of records in return list: ", length(vtimm_first_dose)) )
  
  f_log(" ")
  f_log(" ")
  
  
  # --- Return data ----------------------------------------------------------------------------- #
  return(vtimm_first_dose)
}

# Added November 21
#--------------------------------------------------------------------------------------------------
#
#  Function f_get_vtimm_second_dose
#  ------------------------------
#
#  Get Vaccination to Immunization
#
#
#   Arguments
#   ---------
#     requested_interval_length   Number of values required of the serial interval.
#                                 If this number is higher than the number of records available
#                                 in the input file, the missing records will be set to 0.95

#
#   Return value
#   ------------
#     List 'vtimm_second_dose'
#
#--------------------------------------------------------------------------------------------------
f_get_vtimm_second_dose <- function( requested_interval_length ) {
  
  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Vaccination to immunization")
  f_log(SEP_LINE_2)
  
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start loading Vaccination Immunization data") )
  f_log(SEP_LINE_1)
  f_log( paste0(" Read data from file '", basename(input_vtimm_second_dose), "'") )
  
  
  # --- Read data and create list with the expected length -------------------------------------- #
  
  #     Read file containing information about the 'serial interval'
  vtimm_second_dose_in <- read.csv(input_vtimm_second_dose)
  
  
  #     Select data from column containing the relevant information
  vtimm_second_dose <- vtimm_second_dose_in$effectivity
  
  
  #     Adjust length of 'serial_interval' adding data to or removing data from the read records
  n_vtimm_second_dose_in = length(vtimm_second_dose)
  if (n_vtimm_second_dose_in < requested_interval_length) {
    n_missing_data <- requested_interval_length - n_vtimm_second_dose_in
    vtimm_second_dose[(n_vtimm_second_dose_in+1):requested_interval_length] <- rep(0.95,n_missing_data)
  }
  
  vtimm_second_dose <- vtimm_second_dose[1:requested_interval_length]
  
  
  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Loaded data from file '", basename(input_vtimm_second_dose), "'") )
  f_log(SEP_LINE_1)
  
  f_log( paste0(" Number of records in file:        ", nrow(vtimm_second_dose_in)) )
  f_log( paste0(" Number of records in return list: ", length(vtimm_second_dose)) )
  
  f_log(" ")
  f_log(" ")
  
  
  # --- Return data ----------------------------------------------------------------------------- #
  return(vtimm_second_dose)
}

#--------------------------------------------------------------------------------------------------
#
#  Function f_get_ifr
#  ------------------------------
#
#  Get ifr
#
#
#   Arguments
#   ---------
#
#   Return value
#   ------------
#     List 'serial_interval'
#
#--------------------------------------------------------------------------------------------------
f_get_ifr <- function() {
  
  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" IFR")
  f_log(SEP_LINE_2)
  
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start loading ifr data") )
  f_log(SEP_LINE_1)
  f_log( paste0(" Read data from file '", basename(input_ifr), "'") )
  
  
  # --- Read data and create list with the expected length -------------------------------------- #
  
  #     Read file containing information about the 'serial interval'
  ifr <- read.csv(input_ifr)
 
  
  
  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Loaded data from file '", basename(input_ifr), "'") )
  f_log(SEP_LINE_1)
  
  f_log(" ")
  f_log(" ")
  
  
  # --- Return data ----------------------------------------------------------------------------- #
  return(ifr)
}


#--------------------------------------------------------------------------------------------------
#
#  Function f_get_infection_to_death
#  ---------------------------------
#
#  Get infection to death (ITD)
#
#
#   Arguments
#   ---------
#     requested_itd_length   Number of values required of the infection to death information.
#                            If this number is higher than the number of records available
#                            in the input file, the missing records will be set to 0.
#
#
#   Return value
#   ------------
#     List 'infection_to_death'
#
#--------------------------------------------------------------------------------------------------
f_get_infection_to_death <- function( requested_itd_length ) {

  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Infection to death (ITD)")
  f_log(SEP_LINE_2)

  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start loading infection to death (ITD) data") )
  f_log(SEP_LINE_1)
  f_log( paste0(" Read data from file '", basename(input_infection_to_death), "'") )


  # --- Read data and create list with the expected length -------------------------------------- #

  #     Read file containing information about the 'infection to death'
  infection_to_death_in <- read.csv(input_infection_to_death)


  #     Select data from column containing the relevant information
  infection_to_death <- infection_to_death_in$fit


  #     Adjust length of 'infection_to_death' adding data to or removing data from the read records
  n_infection_to_death_in = length(infection_to_death)
  if (n_infection_to_death_in < requested_itd_length) {
    n_missing_data <- requested_itd_length - n_infection_to_death_in
    infection_to_death[(n_infection_to_death_in+1):requested_itd_length] <- rep(0.,n_missing_data)
  }

  infection_to_death <- infection_to_death[1:requested_itd_length]


  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Loaded data from file '", basename(input_infection_to_death), "'") )
  f_log(SEP_LINE_1)

  f_log( paste0(" Number of records in file:        ", nrow(infection_to_death_in)) )
  f_log( paste0(" Number of records in return list: ", length(infection_to_death)) )

  f_log(" ")
  f_log(" ")


  # --- Return data ----------------------------------------------------------------------------- #
  return(infection_to_death)
}




# 27Mar21 Juan
#--------------------------------------------------------------------------------------------------
#
#  Function f_get_infection_to_hospitalization
#  ------------------------------
#
#  Get infection-to-hospitalization (ITH)
#
#
#   Arguments
#   ---------
#     requested_interval_length   Number of values required of the serial interval.
#                                 If this number is higher than the number of records available
#                                 in the input file, the missing records will be set to 0.

#
#   Return value
#   ------------
#     List 'infection_to_hospitalization'
#
#--------------------------------------------------------------------------------------------------
f_get_infection_to_hospitalization <- function( requested_interval_length ) {

  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Infection to Hospitalization (ITH)")
  f_log(SEP_LINE_2)

  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start calculating infection-to-hospitalization (ITH)") )
  f_log(SEP_LINE_1)


  # --- calculate infection-to-hospitalization (including general wards and ICUs) and create a list with the expected length -------------------------------------- #

    library(mixdist)    # package installed by: install.packages("mixdist")
    # weibullpar(9.1, 4., loc=0)   # to find shape and scale parameters for weibull distribution from mean and sd in NewEnglMed20 (not sure of sd)

    h_ITH=rep(0.,40)
    mean1_ITO = 5.1; cv1_ITO = 0.55; # infection to onset, with more reasonable cv... need to find parameters in literature
    shape_OTH = 2.426088; scale_OTH = 10.26315;

    x1_ITO = rlnorm(5e6,log(mean1_ITO),cv1_ITO)
    x2_OTH = rweibull(5e6,shape_OTH,scale_OTH) #

    fc_ITH = ecdf(x1_ITO+x2_OTH)
    convolution_ITH = function(u) (fc_ITH(u))

    h_ITH[1] = (convolution_ITH(1.5) - convolution_ITH(0))
    for(i_ITH in 2:length(h_ITH)) {
      h_ITH[i_ITH] = (convolution_ITH(i_ITH+.5) - convolution_ITH(i_ITH-.5)) / (1-convolution_ITH(i_ITH-.5))
    }
    s_ITH = rep(0,40)
    s_ITH[1] = 1
    for(i_ITH in 2:40) {
      s_ITH[i_ITH] = s_ITH[i_ITH-1]*(1-h_ITH[i_ITH-1])
    }
    infection_to_hospitalization =  s_ITH * h_ITH



  #     Adjust length of 'infection_to_hospitalization' adding data to or removing data from the read records
  n_infection_to_hospitalization_in = length(infection_to_hospitalization)
  if (n_infection_to_hospitalization_in < requested_interval_length) {
    n_missing_data <- requested_interval_length - n_infection_to_hospitalization_in
    infection_to_hospitalization[(n_infection_to_hospitalization_in+1):requested_interval_length] <- rep(0.,n_missing_data)
  }

  infection_to_hospitalization <- infection_to_hospitalization[1:requested_interval_length]


  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - calculated data ") )
  f_log(SEP_LINE_1)

  f_log( paste0(" Number of records in return list: ", length(infection_to_hospitalization)) )

  f_log(" ")
  f_log(" ")


  # --- Return data ----------------------------------------------------------------------------- #
  return(infection_to_hospitalization)
}
# 27Mar21 Juan







# 27Mar21 Juan
#--------------------------------------------------------------------------------------------------
#
#  Function f_get_hospitalization_to_recovery
#  -----------------------------------------
#
#  Get hospitalization-to-recovery (HTR)
#
#
#   Arguments
#   ---------
#     requested_interval_length   Number of values required of the serial interval.
#                                 If this number is higher than the number of records available
#                                 in the input file, the missing records will be set to 0.

#
#   Return value
#   ------------
#     List 'hospitalization_to_recoery'
#
#--------------------------------------------------------------------------------------------------
f_get_hospitalization_to_recovery <- function( requested_interval_length ) {

  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Hospitalization to Recovery (HTR)")
  f_log(SEP_LINE_2)

  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start calculating hospitalization-to-recovery (HTR)") )
  f_log(SEP_LINE_1)


  # --- calculate hospitalization-to-recovery and create a list with the expected length -------------------------------------- #
  #     (recovered means those who exit from hospital, including general wards and UCIs; for now we assume that 100% of hospitalized recover or exit)

    h_HTR=rep(0.,50)
    shape_HTR = 1.7173; scale_HTR = 11.1293;    # from data in La Rioja, 24Sep-5Nov, male >65 years old not requiring ICU (Angel from CIBIR)... try other averaged values?
    x1_HTR = rweibull(5e6,shape_HTR,scale_HTR)

    x2_HTR = c(rgamma(0.81*5e6,2,2/13.11), rgamma(0.19*5e6,2,2/17.57))  # from data in France, Science 2020, considering proportion of hospital and ICUs in France

    x2b_HTR = c(rgamma(0.86*5e6,2,2/13.11), rgamma(0.14*5e6,2,2/17.57))  # from data in France, Science 2020, considering proportion of hospital and ICUs in Rioja 


    x3_HTR = c(rweibull(0.4325*0.5688*5e6,1.7012,7.2232),rweibull((1-0.4325)*0.5688*5e6,1.7173,11.1293),
               rweibull(0.3717*0.4312*5e6,1.5814,7.7372),rweibull((1-0.3717)*0.4312*5e6,1.4108,10.4412)  )  # from Angel CIBIR, patients not requiring ICU


    x4_HTR = c(rweibull(0.4325*0.5688*(1-0.1284)*5e6,1.7012,7.2232),rweibull((1-0.4325)*0.5688*(1-0.0559)*5e6,1.7173,11.1293),
               rweibull(0.3717*0.4312*(1-0.0086)*5e6,1.5814,7.7372),rweibull((1-0.3717)*0.4312*(1-0.0114)*5e6,1.4108,10.4412),
               rweibull(0.4325*0.5688*(0.1284)*5e6,2.1745,25.3969),rweibull((1-0.4325)*0.5688*(0.0559)*5e6,2.055,11.6422),
               rweibull(0.3717*0.4312*(0.0086)*5e6,2.2624,25.0061),rweibull((1-0.3717)*0.4312*(0.0114)*5e6,2.1359,11.7313)
              )  # from Angel CIBIR, hospitalized not requiring ICU, and direct ICU  (not considered ICU from hospitalized)

    x5_HTR = c(rweibull(0.4325*0.5688*(1-0.1284)*(1-0.1684)*5e6,1.7012,7.2232),rweibull((1-0.4325)*0.5688*(1-0.0559)*(1-0.1111)*5e6,1.7173,11.1293),
               rweibull(0.3717*0.4312*(1-0.0086)*(1-0.0391)*5e6,1.5814,7.7372),rweibull((1-0.3717)*0.4312*(1-0.0114)*(1-0.026)*5e6,1.4108,10.4412),
               rweibull(0.4325*0.5688*(0.1284+0.1684*(1-0.1284))*5e6,2.1745,25.3969),rweibull((1-0.4325)*0.5688*(0.0559+0.1111*(1-0.0559))*5e6,2.055,11.6422),
               rweibull(0.3717*0.4312*(0.0086+0.0391*(1-0.0086))*5e6,2.2624,25.0061),rweibull((1-0.3717)*0.4312*(0.0114+0.026*(1-0.0114))*5e6,2.1359,11.7313)
              )  # from Angel CIBIR, considering direct ICUs and those requiring ICU from ward, but for them only ICU time is considered (not further time in hospital)

    fc_HTR = ecdf(x2_HTR)
    convolution_HTR = function(u) (fc_HTR(u))

    h_HTR[1] = (convolution_HTR(1.5) - convolution_HTR(0))
    for(i_HTR in 2:length(h_HTR)) {
      h_HTR[i_HTR] = (convolution_HTR(i_HTR+.5) - convolution_HTR(i_HTR-.5)) / (1-convolution_HTR(i_HTR-.5))
    }
    s_HTR = rep(0,50)
    s_HTR[1] = 1
    for(i_HTR in 2:50) {
      s_HTR[i_HTR] = s_HTR[i_HTR-1]*(1-h_HTR[i_HTR-1])
    }
    hospitalization_to_recovery =  s_HTR * h_HTR



  #     Adjust length of 'hospitalization_to_recovery' adding data to or removing data from the read records
  n_hospitalization_to_recovery_in = length(hospitalization_to_recovery)
  if (n_hospitalization_to_recovery_in < requested_interval_length) {
    n_missing_data <- requested_interval_length - n_hospitalization_to_recovery_in
    hospitalization_to_recovery[(n_hospitalization_to_recovery_in+1):requested_interval_length] <- rep(0.,n_missing_data)
  }

  hospitalization_to_recovery <- hospitalization_to_recovery[1:requested_interval_length]

  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - calculated data ") )
  f_log(SEP_LINE_1)

  f_log( paste0(" Number of records in return list: ", length(hospitalization_to_recovery)) )

  f_log(" ")
  f_log(" ")


  # --- Return data ----------------------------------------------------------------------------- #
  return(hospitalization_to_recovery)
}
# 27Mar21 Juan





# 16Apr21 Juan
#--------------------------------------------------------------------------------------------------
#
#  Function f_get_infection_to_detection
#  ------------------------------
#
#  Get infection-to-detection (ITDet)
#
#
#   Arguments
#   ---------
#     requested_interval_length   Number of values required of the serial interval.
#                                 If this number is higher than the number of records available
#                                 in the input file, the missing records will be set to 0.
#
#   Return value
#   ------------
#     List 'infection_to_detection'
#
#--------------------------------------------------------------------------------------------------
f_get_infection_to_detection <- function( requested_interval_length ) {

  # --- Start of the processing ----------------------------------------------------------------- #
  f_log(" ")
  f_log(" ")
  f_log(SEP_LINE_2)
  f_log(" Infection to Detecttion (ITDet)")
  f_log(SEP_LINE_2)

  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - Start calculating infection-to-detection (ITDet)") )
  f_log(SEP_LINE_1)


  # --- calculate infection-to-detection and create a list with the expected length -------------------------------------- #

    # infection-to-detection will be same as infection-to-first_medical_visit

    library(mixdist)    # package installed by: install.packages("mixdist")
    # weibullpar(4.6, 4., loc=0)   # to find shape and scale parameters for weibull distribution from mean and sd in NewEnglMed20 (not sure of sd)

    h_ITM=rep(0.,40)
    mean1_ITO = 5.1; cv1_ITO = 0.55; # infection to onset, with more reasonable cv... need to find parmeters in literature
    mean2_OTM = 2.7; cv2_OTM = 1.; # onset to first medical visit (medRxiv20 Linton) ... cv not sure
    shape_OTM = 1.153083; scale_OTM = 4.837206;

    x1_ITO = rlnorm(5e6,log(mean1_ITO),cv1_ITO)
    #x2_OTM = rgammaAlt(5e6,mean2_OTM,cv2_OTM)
    x2_OTM = rweibull(5e6,shape_OTM,scale_OTM) #

    fc_ITM = ecdf(x1_ITO+x2_OTM)
    convolution_ITM = function(u) (fc_ITM(u))

    h_ITM[1] = (convolution_ITM(1.5) - convolution_ITM(0))
    for(i_ITM in 2:length(h_ITM)) {
      h_ITM[i_ITM] = (convolution_ITM(i_ITM+.5) - convolution_ITM(i_ITM-.5)) / (1-convolution_ITM(i_ITM-.5))
    }
    s_ITM = rep(0,40)
    s_ITM[1] = 1
    for(i_ITM in 2:40) {
      s_ITM[i_ITM] = s_ITM[i_ITM-1]*(1-h_ITM[i_ITM-1])
    }
    ITM =  s_ITM * h_ITM


   # infection-to-detection curve (cases infected on i-j day are detected on i day with a given probability) 
   # at the beginning of the outbreak, we can assume that this distribution will follow the same as infection-to-medical
   # asumming that most tests are done when they go to first medical visit

   infection_to_detection=ITM


  #     Adjust length of 'infection_to_detecton' adding data to or removing data from the read records
  n_infection_to_detection_in = length(infection_to_detection)
  if (n_infection_to_detection_in < requested_interval_length) {
    n_missing_data <- requested_interval_length - n_infection_to_detection_in
    infection_to_detection[(n_infection_to_detection_in+1):requested_interval_length] <- rep(0.,n_missing_data)
  }

  infection_to_detection <- infection_to_detection[1:requested_interval_length]


  # --- Summary of data read -------------------------------------------------------------------- #
  f_log(" ")
  f_log(SEP_LINE_1)
  f_log( paste0(" ", Sys.time(), " - calculated data ") )
  f_log(SEP_LINE_1)

  f_log( paste0(" Number of records in return list: ", length(infection_to_detection)) )

  f_log(" ")
  f_log(" ")


  # --- Return data ----------------------------------------------------------------------------- #
  return(infection_to_detection)
}
# 16Apr21 Juan






