library(truncdist)
library(data.table)




sim_quarantine_lnorm <- function(dt_incubation_dists_lnorm,
                                 params,
                                 n_sympt = 6000, n_asympt = 1000,
                                 seed = 91, rand_u = TRUE){
  #Returns 'expected # asymptomatic-infectious days in community' for asymptomatic infectious
  #  travelers and 'expected # presymtomatic-infectious and symptomatic-infectious days
  #  in community' for the average traveller with infection that produces symptoms
  #  for 1000 duration parameter sets at specified quarantine duration and testing
  #  policy
  
  set.seed(seed)
  n_iters = nrow(dt_incubation_dists_lnorm)
  #For each bootstrapped incubation distribution and for a set of 
  #   parameters and policies, calculates the outcomes for 
  #   symptomatic and asymptomatic individuals
  
  #Sample asympt duration gamma parameters
  dt_d_asympt <- data.table(
    mean = runif(n_iters, min = params$dur_asympt_mean_lb, max = params$dur_asympt_mean_ub),
    var = runif(n_iters, min = params$dur_asympt_var_lb, max = params$dur_asympt_var_ub)
  )
  dt_d_asympt[ , par1 := mean^2/var]
  dt_d_asympt[ , par2 := var/mean]
  
  #Sample sympt duration gamma parameters
  dt_d_sympt <- data.table(
    mean = runif(n_iters, min = params$dur_sympt_mean_lb, max = params$dur_sympt_mean_ub),
    var = runif(n_iters, min = params$dur_sympt_var_lb, max = params$dur_sympt_var_ub)
  )
  dt_d_sympt[ , par1 := mean^2/var]
  dt_d_sympt[ , par2 := var/mean]
  
  #Sample presympt duration gamma parameters (truncation implemented in iteration forloop below)
  dt_d_presympt <- data.table(
    mean = runif(n_iters, min = params$dur_presympt_mean_lb, max = params$dur_presympt_mean_ub),
    var = runif(n_iters, min = params$dur_presympt_var_lb, max = params$dur_presympt_var_ub)
  )
  dt_d_presympt[ , par1 := mean^2/var]
  dt_d_presympt[ , par2 := var/mean]
  
  #table to store output
  dt_raw_metrics <- data.table(
    iter = numeric(),
    dur_q = numeric(),
    testing = character(),
    d_presympt = numeric(),
    d_sympt = numeric(),
    d_asympt = numeric()
  )
  
  # Iterate over sets of duration parameters
  for (row in 1:nrow(dt_incubation_dists_lnorm)){
    #[PRE]SYMPOMATIC INFECTION
    #Sample all durations
    dt_sympt <-
      data.table(
        dur_presympt = rtrunc(n_sympt, "gamma", 0.8, 3, 
                              shape=dt_d_presympt[row, par1], scale=dt_d_presympt[row, par2]),
        dur_sympt = rgamma(n_sympt, 
                           shape=dt_d_sympt[row, par1], scale = dt_d_sympt[row, par2]),
        dur_incubation = rlnorm(n_sympt, 
                                meanlog=dt_incubation_dists_lnorm[row, par1], 
                                sdlog = dt_incubation_dists_lnorm[row, par2])
      )
    
    #Calc times. t=0 is time of infection
    dt_sympt[ , t_infectious_start := pmax(0, dur_incubation - dur_presympt)]
    dt_sympt[ , t_sympt_start := dur_incubation]
    dt_sympt[ , t_recovery := dur_incubation + dur_sympt]
    if(rand_u == TRUE){
      dt_sympt[ , dur_infected_prequarantine := runif(n_sympt)*t_recovery]
    } else {
      dt_sympt[ , dur_infected_prequarantine := 0]
    }
    
    
    
    #ASYMPTOMATIC INFECTION
    #Sample all durations
    dt_asympt <-
      data.table(
        dur_presympt = rtrunc(n_asympt, "gamma", 0.8, 3, 
                              shape=dt_d_presympt[row, par1], scale=dt_d_presympt[row, par2]),
        dur_infectious = rgamma(n_asympt, 
                                shape=dt_d_asympt[row, par1], scale = dt_d_asympt[row, par2]),
        dur_incubation = rlnorm(n_asympt, 
                                meanlog=dt_incubation_dists_lnorm[row, par1], 
                                sdlog = dt_incubation_dists_lnorm[row, par2])
      )
    
    #Calc times
    dt_asympt[ , t_infectious_start := pmax(0, dur_incubation - dur_presympt)]
    dt_asympt[ , t_recovery := t_infectious_start + dur_infectious]
    if(rand_u == TRUE){
      dt_asympt[ , dur_infected_prequarantine := runif(n_asympt)*t_recovery]
    } else {
      dt_asympt[ , dur_infected_prequarantine := 0]
    }
    #BOTH
    #Calculate outcomes
    dt_metrics_this_iter <-
      rbind(
        data.table(
          iter = row,
          dur_q = pols_dur_quarantine,
          testing = 0,
          d_presympt = unlist(lapply(X = pols_dur_quarantine, FUN = calc_presympt_days_in_community, 
                                     dt_sympt = dt_sympt, params=params, test = 0)),
          d_sympt = unlist(lapply(X = pols_dur_quarantine, FUN = calc_sympt_days_in_community, 
                                  dt_sympt = dt_sympt, params=params, test = 0)),
          d_asympt = unlist(lapply(X = pols_dur_quarantine, FUN = calc_asympt_days_in_community, 
                                   dt_asympt = dt_asympt, params=params, test = 0))),
        data.table(
          iter = row,
          dur_q = pols_dur_quarantine,
          testing = 1,
          d_presympt = unlist(lapply(X = pols_dur_quarantine, FUN = calc_presympt_days_in_community, 
                                     dt_sympt = dt_sympt, params=params, test = 1)),
          d_sympt = unlist(lapply(X = pols_dur_quarantine, FUN = calc_sympt_days_in_community, 
                                  dt_sympt = dt_sympt, params=params, test = 1)),
          d_asympt = unlist(lapply(X = pols_dur_quarantine, FUN = calc_asympt_days_in_community, 
                                   dt_asympt = dt_asympt, params=params, test = 1)))
      )
    
    dt_raw_metrics <- rbind(dt_raw_metrics,
                            dt_metrics_this_iter)
  }
  return(dt_raw_metrics)
}


calc_presympt_days_in_community <- function(dur_q, dt_sympt, params, test = 0){
  if(test == 0){
    return(
      dt_sympt[ , sum(
        pmax(0,(1-params$prob_quarantine_compliance)*
               (t_sympt_start - pmax(t_infectious_start, dur_infected_prequarantine)))+
          pmax(0,params$prob_quarantine_compliance*
                 (t_sympt_start - pmax(t_infectious_start, dur_q+dur_infected_prequarantine))
          )
      )/nrow(dt_sympt)]
    )
  } else {
    return(
      dt_sympt[ , sum(
        #Noncompliant with quarantine, not tested
        pmax(0,(1-params$prob_quarantine_compliance)*
               (t_sympt_start - pmax(t_infectious_start, dur_infected_prequarantine)))+
          params$prob_quarantine_compliance*
          pmax(0,
               ifelse(dur_q+dur_infected_prequarantine-1 < t_infectious_start,
                      #tested before infectious
                      t_sympt_start - pmax(dur_infected_prequarantine+dur_q, t_infectious_start),
                      #tested while presymtomatic or after
                      ((1-params$sn_presympt)+params$sn_presympt*(1-params$prob_isolate_test))* #prob test neg or test pos & refuse to isolate
                        (t_sympt_start - (dur_q+dur_infected_prequarantine)))
          )
      )/nrow(dt_sympt)]
    )
    
  }
}


calc_sympt_days_in_community <- function(dur_q, dt_sympt, params, test = 0){
  if(test == 0){
    return(
      dt_sympt[ , sum(
        #probability no isolation with symptoms
        (1-params$prob_isolate_sympt)*
          (#noncompliant with quarantine
            pmax(0,(1-params$prob_quarantine_compliance)*(t_recovery - pmax(t_sympt_start, dur_infected_prequarantine)))+
              #compliant with quarantine
              params$prob_quarantine_compliance*
              pmax(0,
                   t_recovery - pmax(t_sympt_start, dur_infected_prequarantine+dur_q)
              )
          )
      )/nrow(dt_sympt)]
      
      
    )
  } else {
    return(
      dt_sympt[ , sum(
        #noncompliant
        (1-params$prob_quarantine_compliance)*(1-params$prob_isolate_sympt)*
          pmax(0, t_recovery - pmax(t_sympt_start, dur_infected_prequarantine))+
          params$prob_quarantine_compliance*
          pmax(0,
               ifelse((dur_q+dur_infected_prequarantine - 1 < t_infectious_start),
                      #tested before infectious therefore neg test
                      (t_recovery - pmax(t_sympt_start, dur_q+dur_infected_prequarantine))*
                        (1-params$prob_isolate_sympt),
                      ifelse((dur_q+dur_infected_prequarantine - 1 < t_sympt_start),
                             #tested during presymptomatic
                             #  tested negative
                             (1-params$sn_presympt)*(1-params$prob_isolate_sympt)*
                               (t_recovery - pmax(t_sympt_start, dur_q+dur_infected_prequarantine)) +
                               #  tested positive
                               params$sn_presympt*(1-params$prob_isolate_both)*
                               (t_recovery - pmax(t_sympt_start, dur_q+dur_infected_prequarantine)),
                             #tested while symptomatic or after
                             #  tested negative
                             (1-params$sn_sympt)*(1-params$prob_isolate_sympt)*
                               (t_recovery - (dur_q + dur_infected_prequarantine))+
                               #  tested positive
                               params$sn_sympt*(1-params$prob_isolate_both)*
                               (t_recovery - (dur_q + dur_infected_prequarantine))
                      )
               )
          )
        
      )/nrow(dt_sympt)]
    )
  }
}


calc_asympt_days_in_community <- function(dur_q, dt_asympt, params, test = 0){
  if(test == 0){
    return(
      dt_asympt[ , sum(
        pmax(0,(1-params$prob_quarantine_compliance)*
               (t_recovery - pmax(dur_infected_prequarantine, t_infectious_start)))+
          pmax(0, params$prob_quarantine_compliance*
                 (t_recovery - pmax(t_infectious_start, dur_q+dur_infected_prequarantine)))
      )/nrow(dt_asympt)]
    )
  } else {
    return(
      dt_asympt[ , sum(
        #No compliance & no test
        pmax(0, (1-params$prob_quarantine_compliance)*
               (t_recovery - pmax(dur_infected_prequarantine, t_infectious_start)))+
          params$prob_quarantine_compliance*pmax(0, 
                                                 ifelse(dur_q+dur_infected_prequarantine - 1 < t_infectious_start,
                                                        #tested before infectious/detectable
                                                        t_recovery - pmax(t_infectious_start, dur_infected_prequarantine+dur_q),
                                                        #tested while or after infectious
                                                        ((1-params$sn_asympt)+params$sn_asympt*(1-params$prob_isolate_test))* #Prob test neg or test pos & refuse to quarantine
                                                          (t_recovery - (dur_q+dur_infected_prequarantine))
                                                 )
          )
      )/nrow(dt_asympt)]
    )
    
  }
}  


sim_calc_metrics <- function(dt_raw_metrics,
                             params,
                             prev_levels,
                             scenario = "basecase"){
  
  dt_raw_metrics[ , d_tot := (d_asympt*params$prob_asympt +
                                (1-params$prob_asympt)*(d_sympt+d_presympt))
                  ]
  
  dt_metrics <- data.table(
    scenario = character(),
    prev = numeric(),
    metric = character(),
    dur_q = numeric(),
    testing = numeric(),
    quantile = character(),
    value = numeric()
  )
  
  
  quant_probs = c(.01, .5, .99)
  for(prev in prev_levels){
    dt_metrics <- rbind(dt_metrics,
                        dt_raw_metrics[ , list(value= prev*quantile(d_tot, probs = quant_probs), 
                                               scenario = scenario,
                                               prev = prev,
                                               metric = "Days at risk per traveler",
                                               quantile = paste0("q",quant_probs)), by = c("dur_q", "testing")]
    )
  }
  return(dt_metrics)
}
