//***************************************************************************************************//
//
//   Stan model to estimate the evolution of COVID-19 cases and deaths
//
//***************************************************************************************************//

//***************************************************************************************************//
//  data
//***************************************************************************************************//
data {

    int<lower=1>                N_DATE_0;                       // number of days for which to impute infections
    int<lower=1>                N_PANDEMIC_DATES;               // days of observed data for country since start of the pandemic     // ADDED 6Mar21 Juan
    int<lower=1>                N_DATES;                        // days of observed + forecast data      // ADDED 6Mar21 Juan
    real                        SI[N_DATES];                    // fixed pre-calculated SI using emprical data from Neil
    int <lower=1>               N_MEASURES;                     // number of non pharmacological measures
    int                         reported_deaths[N_DATES];       // reported deaths
    int                         cum_reported_deaths[N_DATES];   // cumulative deaths       ADDED Juan 6Mar21
    real<lower=0>               fatality_rate[N_DATES];         // fatality rate
    matrix[N_MEASURES, N_DATES] response_measures_country;      // information about the active measures in each date
    int<lower=1>                population;                     // to impose that number of cumulative infected should be lower than population    ADDED 7Mar21 Juan
    real<lower=0>               hospitalization_rate[N_DATES];  // hospitalization rate
    real<lower=0>               HTR[N_DATES];                   // hospitalization to recovery
    int                         reported_cases[N_DATES];        // reported cases   ADDED Juan 9Apr21
    int                         reported_hospitalized[N_DATES]; // reported hospitalized   ADDED Juan 9Apr21
    real<lower=0>               ITDet[N_DATES];                 // infection to detection
    int                         reported_IA14[N_DATES];        // IA14   ADDED Juan 17Apr21    
    int                         cum_reported_cases[N_DATES];       // cumulative reported cases
    real                        ITD[N_DATES];                           // ITD ADDED Alberto 25Nov21
    real                        CFR_by_day[N_DATES];           // CFR_by_DAY
    int                         vaccination_total_by_day_acum[N_DATES]; // ADDED Alberto 25Nov21
}


//***************************************************************************************************//
//  parameters
//***************************************************************************************************//
parameters {

    real<lower=0>              mu;                 // intercept for Rt
    real<lower=0>              kappa;
    real<lower=1e-9>           relR[N_MEASURES];  // array with the contribution of each measure to the reduction or increment of Rt
    real<lower=0>              y_cases;            // cases
    real<lower=0>              phi;
    real<lower=0>              tau;
    real<lower=1e-9, upper=1>  detRatio;       // detection ratio
}


//***************************************************************************************************//
//  transformed parameters
//***************************************************************************************************//
transformed parameters {

    real Rt [N_DATES]                            = rep_array(0.0,N_DATES);    // Rt estimated for a specific date
    real Rt_eff [N_DATES]                        = rep_array(0.0,N_DATES);    // Effective Rt considering susceptible population   // 21Mar21
    real prediction_cases [N_DATES]              = rep_array(1e-9,N_DATES);    // New cases predicted for a specific date
    real prediction_deaths [N_DATES]             = rep_array(1e-9,N_DATES);    // Deaths predicted for a specific date
    real cum_prediction_cases [N_DATES]          = rep_array(1e-9,N_DATES);    // Cumulative number of cases predicted until a specific date
    real cum_prediction_deaths [N_DATES]         = rep_array(1e-9,N_DATES);    // Cumulative number of deaths predicted until a specific date
    real prediction_newhospitalized [N_DATES]    = rep_array(1e-9,N_DATES);    // Daily new hospitalized (including general wards and ICUs)     
    real cum_prediction_newhospitalized [N_DATES]    = rep_array(1e-9,N_DATES);    // Cumulative number of new hospitalized (including general wards and ICUs) until a specific date
    real prediction_recovered_hosp [N_DATES]     = rep_array(1e-9,N_DATES);    // Daily recovered from hospitalization (including general wards and ICUs)
    real cum_prediction_recovered_hosp [N_DATES] = rep_array(1e-9,N_DATES);    // Cumulative number of recovere from hospitalization (including general wards and ICUs)
    real prediction_hospitalized [N_DATES]       = rep_array(1e-9,N_DATES);    // Total patients hospitalized in a given day (including general wards and ICUs)
    real prediction_IA14 [N_DATES]               = rep_array(1e-9,N_DATES);    // 14-day cumulative incidence per 100K inhabitants
    real prediction_detected [N_DATES]           = rep_array(1e-9,N_DATES);    // New cases predicted to be detected for a specific date



    // --- Rt -------------------------------------------------------------------------------------- //
    for (ind_date in 1:N_DATES) {
      real added_measure_effect = 1;

      for (ind_measure in 1:N_MEASURES) {
        if (response_measures_country[ind_measure, ind_date] >= 1){      // 21Mar21  even if there are overlapping periods with same measure, it counts only once
           added_measure_effect = added_measure_effect*response_measures_country[ind_measure, ind_date]*relR[ind_measure];    // 8Mar21 Juan
        }
      }

      Rt[ind_date] = mu * added_measure_effect;
    }
    Rt_eff = Rt;    // 21Mar21  to initialize effective Rt values for the first days

    prediction_cases[1:N_DATE_0] = rep_array(y_cases,N_DATE_0); // learn the number of cases in the first N_DATE_0 days
    cum_prediction_cases[1] = prediction_cases[1];
    for (i in 2:N_DATES){         // 28Mar21
      Rt_eff[i] = Rt[i] * (1-((cum_prediction_cases[i-1]+vaccination_total_by_day_acum[i-1])/population)); # habra que añadir en el numerador los inmunizados de ese dia por vacunacion
       // Rt diminishes with susceptible   21Mar21
      for(j in 1:(i-1)) {
        if (i > N_DATE_0){
           
           prediction_cases[i] += Rt_eff[i] * prediction_cases[j]*SI[i-j];   
           prediction_detected[i] += prediction_cases[j]*ITDet[i-j];       
        }
        //prediction_deaths[i] += prediction_cases[j]*fatality_rate[i-j];
        prediction_deaths[i] += prediction_cases[j]*CFR_by_day[i-j]*ITD[i-j]; //Added Alberto
        prediction_newhospitalized[i] += prediction_cases[j]*hospitalization_rate[i-j];
        prediction_recovered_hosp[i] += prediction_newhospitalized[j]*HTR[i-j];
      }
      if ((cum_prediction_cases[i-1] + prediction_cases[i]) > population){
            prediction_cases[i] = population - cum_prediction_cases[i-1];
      }
      cum_prediction_cases[i] = cum_prediction_cases[i-1] + prediction_cases[i];
      cum_prediction_deaths[i] = cum_prediction_deaths[i-1] + prediction_deaths[i];
      cum_prediction_newhospitalized[i] = cum_prediction_newhospitalized[i-1] + prediction_newhospitalized[i];
      cum_prediction_recovered_hosp[i] = cum_prediction_recovered_hosp[i-1] + prediction_recovered_hosp[i];
      prediction_hospitalized[i] = cum_prediction_newhospitalized[i] - cum_prediction_recovered_hosp[i];
      for (k in 1:14){
        if (i >= k){
          // prediction_IA14[i] += 100000.*prediction_cases[i-k+1]/population;   // OLD 16Apr21
          prediction_IA14[i] += 100000.*prediction_detected[i-k+1]/population;   // NEW 16Apr21
        }
      }
    }


}

//***************************************************************************************************//
//  Model
//***************************************************************************************************//
model {

    tau     ~ exponential(0.03);
    y_cases ~ exponential(1.0/tau);

    kappa   ~ normal(0,0.5);
    mu      ~ normal(2.4, kappa);

    relR   ~ gamma(1.,1);  // 6Mar21 Juan
    detRatio ~ normal(0.75, kappa);    // 17Apr21
    
    phi     ~ normal(0,5);
    for(i in 1:N_PANDEMIC_DATES){        // ADDED Juan 6Mar21
       reported_deaths[i] ~ neg_binomial_2(prediction_deaths[i],phi);
       cum_reported_deaths[i] ~ neg_binomial_2(cum_prediction_deaths[i],phi);    // ADDED Juan 6Mar21
       // reported_cases[i] ~ neg_binomial_2(detRatio*prediction_detected[i],phi);   // 17Apr21   TEST
       // cum_reported_cases[i] ~ neg_binomial_2(cum_prediction_cases[i],phi);    // 18Apr21
       // reported_IA14[i] ~ neg_binomial_2(detRatio*prediction_IA14[i],phi);   // 17Apr21   TEST
       // reported_hospitalized[i] ~ neg_binomial_2(prediction_hospitalized[i],phi);  // 9Apr21
    }
}

//***************************************************************************************************//
