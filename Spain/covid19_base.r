#**************************************************************************************************
#
#  Script to estimate the evolution of COVID-19 cases and deaths based on 
#  existing data and on a Stan models
#
#**************************************************************************************************
suppressWarnings(suppressMessages({
  library(rstan)
  library(data.table)
  library(lubridate)
  library(gdata)
  library(EnvStats)
}))

rstan::stan_version()
#**************************************************************************************************
#  Read country to process from first script argument
#**************************************************************************************************

args <- commandArgs(trailingOnly = TRUE)
country_selected <- trim( args[1] )

n_forecast = 90   # number of days after last day with data, for forecast;   ADDED 6Mar21 Juan   PENDING: read from arguments?
n_prior = 30      # number of days added before first day with data (e.g. use when first day of data already have deaths)     ADDED 26May21

# --- Stop execution if no country provided ----------------------------------------------------- #
if ( is.na(country_selected) || country_selected=="" ) {
  
  error_str <- " Script executed without arguments - country to be processed required as argument"
  print(error_str)
  print(" --> Stop execution")
  print(" ")
  print(" Example of correct script execution: 'Rscript.exe covid19_base.r Spain'")
  
  stop(error_str)
}


#**************************************************************************************************
#  Definition of the different directories and files 
#**************************************************************************************************

# --- Directories ------------------------------------------------------------------------------- #
INPUT_DATA_DIR       <- paste0(getwd(), "/data")
INPUT_STAN_MODEL_DIR <- paste0(getwd(), "/stan-models")
OUTPUT_RESULTS_DIR   <- paste0(getwd(), "/results")


# --- Stan model -------------------------------------------------------------------------------- #

# --- If we want analyze evaluate the impact of vaccination
STAN_MODEL           <- "covid19_base_mod"

# --- If we want analyze only response measures 
# STAN_MODEL           <- "covid19_base"

STAN_MODEL_FILE      <- paste0(INPUT_STAN_MODEL_DIR, "/", STAN_MODEL, ".stan")


# --- Unique ID to identify each run of the script (will be used as part of the output file names) #
RUN_ID_STR           <- paste0(country_selected, "_", STAN_MODEL, "_", Sys.time())
RUN_ID_STR           <- trim(RUN_ID_STR)
RUN_ID_STR           <- gsub(" ", "_", RUN_ID_STR)
RUN_ID_STR           <- gsub(":", "-", RUN_ID_STR)


# --- Output files ------------------------------------------------------------------------------ #
RDATA_FILE           <- paste0(OUTPUT_RESULTS_DIR, "/", RUN_ID_STR, '.Rdata')
GRAPHICS_FILE        <- paste0(OUTPUT_RESULTS_DIR, "/", RUN_ID_STR, '.pdf')
LOG_FILE             <- paste0(OUTPUT_RESULTS_DIR, "/", RUN_ID_STR, '.log')


#--------------------------------------------------------------------------------------------------
#  Functions used to get data from different files, create graphics, etc.
#  These are defined in a separate file to make this main script more readable.
#--------------------------------------------------------------------------------------------------
source("covid19_file_functions.r")


#--------------------------------------------------------------------------------------------------
#  Functions to create graphics
#--------------------------------------------------------------------------------------------------
source("covid19_graphic_functions.r")


#**************************************************************************************************
#  Start execution
#**************************************************************************************************
start_time <- Sys.time()

f_log_heading("START EXECUTION")
f_log( paste0(" Processing for country '", country_selected, "'") )
f_log( paste0(" Execution identified by the string '", RUN_ID_STR, "'") )


#--------------------------------------------------------------------------------------------------
#  Read data required for the processing
#--------------------------------------------------------------------------------------------------
f_log_heading("READ DATA")

# --- COVID-19 data (cases, deaths, ...) for the country to be analyzed ------------------------- #
covid19_country_data <- f_get_covid_data_for_country(country_selected,n_deaths_0=10)     # TEST 18Apr21

#     Stop execution if no data available for the country
if ( is.null(covid19_country_data) ) {
  
  error_str <- paste0(" No data found for country '", country_selected, "' --> stop execution") 
  f_log(error_str)
  f_log(" ")
  
  stop(error_str)
}


pandemic_dates   <- covid19_country_data$date
n_pandemic_dates <- length(pandemic_dates) 

pandemic_dates_forecast    <- pandemic_dates[n_pandemic_dates] + 1:n_forecast   # ADDED 6Mar21 Juan
pandemic_dates_prior       <- (pandemic_dates[1] - (n_prior:0) )[0:n_prior]       # ADDED 26May21
pandemic_dates             <- c(pandemic_dates_prior,pandemic_dates,pandemic_dates_forecast)      # ADDED 26May21
n_dates                    <- length(pandemic_dates)     # ADDED 6Mar21 Juan
n_pandemic_dates           <- n_pandemic_dates + n_prior      # ADDED 26May21

reported_cases           <- c(rep(0,n_prior),as.vector(as.numeric(covid19_country_data$Cases)),rep(0,n_forecast))   # ADDED 6Mar21 Juan
cum_reported_cases       <- c(rep(0,n_prior),as.vector(as.numeric(cumsum(covid19_country_data$Cases))),rep(0,n_forecast))   # 21Mar21
reported_deaths          <- c(rep(0,n_prior),as.vector(as.numeric(covid19_country_data$Deaths)),rep(0,n_forecast))  # ADDED 6Mar21 Juan
cum_reported_deaths      <- c(rep(0,n_prior),as.vector(as.numeric(cumsum(covid19_country_data$Deaths))),rep(0,n_forecast))   # ADDED 6Mar21 Juan


if (any(names(covid19_country_data) == "NewHospitalized")){
  reported_newhospitalized    <- c(rep(0,n_prior),as.vector(as.numeric(covid19_country_data$NewHospitalized)),rep(0,n_forecast))   # ADDED 28Mar21 Juan
}else{
  reported_newhospitalized  <- rep(0.,n_dates)  
}

if (any(names(covid19_country_data) == "Hospitalized")){
  reported_hospitalized     <- c(rep(0,n_prior),as.vector(as.numeric(covid19_country_data$Hospitalized)),rep(0,n_forecast))   # ADDED 28Mar21 Juan
}else{
  reported_hospitalized     <- rep(0.,n_dates)   
}


# --- Measures (interventions) defined for the country ------------------------------------------ #
response_measures_country <- f_get_measures_for_country(country_selected, pandemic_dates)


# --- Case fatality rate ------------------------------------------------------------------------ #
CFR <- f_get_case_fatality_rate_for_country(country_selected)

# --- Vaccination Complete ------------------------------------------------------------------------ #
vaccination_country_complete <- f_get_vaccination_for_country(country_selected)

# --- Data first Dose ------------------------------------------------------------------------ #

vaccination_country_first_dose <- f_get_first_dose(vaccination_country_complete,country_selected,pandemic_dates)

# --- Data second Dose ------------------------------------------------------------------------ #
vaccination_country_second_dose <- f_get_second_dose(vaccination_country_complete,country_selected,pandemic_dates)

# --- Data Population and group age with Dates ----------------------------------------------------------------------- #
population_by_age <- f_get_population_by_age(vaccination_country_complete,country_selected,pandemic_dates)

# --- Data Population and group age with Dates ----------------------------------------------------------------------- #
vaccination_total_by_day_acum <- 
  vaccination_country_first_dose %>% 
  filter(date%in%pandemic_dates) %>% 
  filter(age_group=="Total") %>% 
  select(first_dose)
  
vaccination_total_by_day_acum <-
  as.numeric(paste0(vaccination_total_by_day_acum$first_dose))


# --- Population ------------------------------------------------------------------------ #           # ADDED 7Mar21 Juan
# population_K <- f_get_population_for_country(country_selected)                                      # ADDED 7Mar21 Juan
# population <- 1000*population_K    
population <- 
  population_by_age %>% 
  filter(age_group=="Total") %>% 
  select(Poblacion) %>% 
  distinct()
population <- 
  as.numeric(paste0(population))
# PENDING: here we should distinguish initial population for IA14 purposes, from susceptible population (initial population - infected before series start) 17Apr21


# ---- Cumulative Incidence ----------------------------------------------------------------#   1Apr21
reported_IA14 <- rep(0,n_pandemic_dates)     #1Apr21
for (i in (1:n_pandemic_dates)){
  k = i -13
  k = max(k,1)
  reported_IA14[i] = as.integer(round( 100000.*sum(reported_cases[k:i])/population ))   
  # PENDING: here population should be the entire population, without removing number of infected before series start 
}
reported_IA14 <- c(reported_IA14,rep(0,n_forecast))     #1Apr21
# PENDING: perhaps this function should be defined in covid19_file_functions.r, but since it is very simple, we can leave it here?   17Apr21


# --- Serial interval --------------------------------------------------------------------------- #
SI  <- f_get_serial_interval(n_dates)                                                                 # 6Mar21 Juan


# --- Infection to death ------------------------------------------------------------------------ #
ITD <- f_get_infection_to_death(n_dates)                                                              # 6Mar21 Juan

# --- Vaccination to Inmunization first dose----------------------------------------------------- #
VTIMM_FIRST  <- f_get_vtimm_first_dose(n_dates)                                                                 # 22Nov21 Alberto

# --- Population to Inmunization first dose----------------------------------------------------- #

population.inmunizated.by.day <-
  f_get_population_immunized(country_selected,pandemic_dates,vaccination_country_first_dose,VTIMM_FIRST,population_by_age)# 22Nov21 Alberto

# --- CFR by day----------------------------------------------------- #

IFR <- f_get_ifr() # 22Nov21 Alberto

# --- CFR by day----------------------------------------------------- #
CFR_by_day <-
  f_get_CFR_by_day(country_selected,pandemic_dates,population.inmunizated.by.day,population_by_age,IFR) # 22Nov21 Alberto

CFR_by_day <-
  unlist(CFR_by_day)

# --- Calculate fatality rate using data obtained in previous steps ----------------------------- #
fatality_rate <- CFR * ITD


# --- Case hospitalization rate ------------------------------------------------------------------------ #
CHR <- f_get_case_hospitalization_rate(country_selected)                              # 27Mar21 Juan


# --- Infection to hospitalization ------------------------------------------------------------------------ #
ITH <- f_get_infection_to_hospitalization(n_dates)                                                          # 27Mar21 Juan


# --- Calculate hospitalization rate using data obtained in previous steps ----------------------------- #
hospitalization_rate <- CHR * ITH


# --- Hospitalization to recovery ------------------------------------------------------------------------ #
HTR <- f_get_hospitalization_to_recovery(n_dates)                                                          # 27Mar21 Juan   for now it assummes that 100% hospitalized recover


# --- Infection to detection ------------------------------------------------------------------------ #
ITDet <- f_get_infection_to_detection(n_dates)                                                          # 16Apr21 Juan




#--------------------------------------------------------------------------------------------------
#  Start model creation using Stan model
#--------------------------------------------------------------------------------------------------
f_log_heading("CREATE MODEL")

# --- Set options ------------------------------------------------------------------------------- #
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# --- Data for Stan model ----------------------------------------------------------------------- #
stan_data = list()

stan_data$N_DATE_0                  <- 6
stan_data$N_PANDEMIC_DATES          <- n_pandemic_dates   # ADDED 6Mar21 Juan
stan_data$N_DATES                   <- n_dates            # ADDED 6Mar21 Juan
stan_data$SI                        <- SI
stan_data$N_MEASURES                <- nrow(response_measures_country)
stan_data$reported_deaths           <- reported_deaths         # ADDED Juan 6Mar21
stan_data$cum_reported_deaths       <- cum_reported_deaths     # ADDED Juan 6Mar21
stan_data$fatality_rate             <- fatality_rate
stan_data$hospitalization_rate      <- hospitalization_rate    # ADDED Juan 27Mar21
stan_data$HTR                       <- HTR                     # ADDED Juan 27Mar21
stan_data$population                <- population              # ADDED Juan 7Mar21
stan_data$response_measures_country <- f_get_numeric_info_from_measures(response_measures_country)
stan_data$y_cases                   <- as.vector(as.numeric(covid19_country_data$Cases))
stan_data$reported_cases             <- reported_cases         # ADDED Juan 9Apr21
stan_data$reported_hospitalized      <- reported_hospitalized         # ADDED Juan 9Apr21
stan_data$ITDet                     <- ITDet                   # ADDED Juan 16Apr21
stan_data$reported_IA14             <- reported_IA14           # 17Apr21
stan_data$cum_reported_cases       <- cum_reported_cases     # ADDED Juan 18Apr21
stan_data$ITD                       <- ITD                # ADDED Alberto 25Nov21
stan_data$CFR_by_day                <- CFR_by_day                # ADDED Alberto 25Nov21
stan_data$vaccination_total_by_day_acum <- vaccination_total_by_day_acum # ADDED Alberto 25Nov21

# --- Create Stan model ------------------------------------------------------------------------- #
f_log( paste0(" ", Sys.time(), " Compile Stan model '", basename(STAN_MODEL_FILE), "'") )
f_log(" ")

covid19_stan_model <- stan_model(STAN_MODEL_FILE)


# --- Fit model (sampling) ---------------------------------------------------------------------- #

#------------------#---------------------#--------------------#
# TEST 1           #  TEST 2             #  Real execution    #
#------------------#---------------------#--------------------#
# n_iter   <- 70   #  n_iter   <- 200    #  n_iter   <- 4000  #
# n_warmup <- 30   #  n_warmup <- 100    #  n_warmup <- 2000  #
# n_chains <- 1    #  n_chains <- 4      #  n_chains <- 4     #
#------------------#---------------------#--------------------#
n_iter   <- 200
n_warmup <- 100
n_chains <- 4


f_log( paste0(" ", Sys.time(), " Fit model - start sampling (iter=", n_iter,
                               " / warmup=", n_warmup, " / chains=", n_chains, ")") )
f_log(" ")

covid19_model_fit = sampling( covid19_stan_model, data=stan_data,
                                                  iter=n_iter, warmup=n_warmup, chains=n_chains,
                                                  thin=4,
                                                  control=list(adapt_delta=0.90, max_treedepth=10) )

f_log( paste0(" ", Sys.time(), " Fit model - sampling finished") )
f_log(" ")


# --- Save data created in a '.Rdata' file ------------------------------------------------------ #     #21Mar21
save.image(RDATA_FILE)

f_log( paste0(" ", Sys.time(), " Results of the processing saved in '", basename(RDATA_FILE), "'") )
f_log(" ")



#--------------------------------------------------------------------------------------------------
#  Create graphics
#--------------------------------------------------------------------------------------------------
f_log_heading("CREATE GRAPHICS")
f_graph_create_plots(GRAPHICS_FILE)



#**************************************************************************************************
#  Execution finished - summary
#**************************************************************************************************
f_log_heading("END EXECUTION - SUMMARY")

f_log( paste0(" Data read from input directory       '", INPUT_DATA_DIR,            "'") )
f_log(" ")                                                                          
                                                                                    
f_log( paste0(" Stan model read from directory       '", INPUT_STAN_MODEL_DIR,      "'") )
f_log( paste0(" Stan model                           '", basename(STAN_MODEL_FILE), "'") )
f_log(" ")                                                                          
                                                                                    
f_log( paste0(" Procesing results saved in directory '", OUTPUT_RESULTS_DIR,        "'") )
f_log( paste0(" Processing data                      '", basename(RDATA_FILE),      "'") )
f_log( paste0(" Graphics file                        '", basename(GRAPHICS_FILE),   "'") )
f_log( paste0(" Log file                             '", basename(LOG_FILE),        "'") )
f_log(" ")

processing_time     <- Sys.time() - start_time
processing_time_str <- paste0( round(processing_time[[1]], 2), " ", units(processing_time) )
f_log( paste0(" Processing time                       ", processing_time_str) )
f_log(" ")
f_log(" ")


#**************************************************************************************************
#**************************************************************************************************
