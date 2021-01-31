incubation_dist_fit_lnorm <- readRDS("../incubation_dists/ncov_inc_fit_boot.rds")
dt_incubation_dists_lnorm <- data.table(incubation_dist_fit_lnorm@samples)



params <- list(
  prob_asympt = 0.4,
  prob_isolate_test = 0.9,
  prob_isolate_sympt = 0.8,
  prob_isolate_both = 1.0,
  sn_presympt = 0.7,
  sn_sympt = 0.7,
  sn_asympt = 0.6,
  prob_quarantine_compliance = 0.8,
  dur_presympt_mean_lb = 1.8,
  dur_presympt_mean_ub = 2.8,
  dur_presympt_var_lb = 4.0,
  dur_presympt_var_ub = 6.0,
  dur_presympt_var_ub = 6.0,
  dur_presympt_truncmin_lb = .64,
  dur_presympt_truncmin_ub = .96,
  dur_presympt_truncmax_lb = 2.4,
  dur_presympt_truncmax_ub = 3.6,
  dur_sympt_mean_lb = 2.6,
  dur_sympt_mean_ub = 3.9,
  dur_sympt_var_lb = 3.0,
  dur_sympt_var_ub = 4.5,
  dur_asympt_mean_lb = 4.0,
  dur_asympt_mean_ub = 6.0,
  dur_asympt_var_lb = 4.0,
  dur_asympt_var_ub = 6.0,
  dur_latent_min = 0.5,
  n_iters = 1000,
  dur_quarantine_alone = c(0, 2, 5, 7, 14),
  dur_quarantine_endtest = c(2, 5, 7),
  dur_quarantine_arrivetest = c(0, 2, 5, 7),
  dur_quarantine_pretest = c(0, 2, 5, 7),
  seed = 91,
  infection_timing = "rand_presympt",
  rr_asympt = 1 #Relative risk for transmission from asympt vs. sympt person
)


n_sympt = 6000
n_asympt = 1000
#print("run_sim")

set.seed(params$seed)


#n_iters = params$n_iters

#Sample asympt duration gamma parameters
dt_d_asympt <- data.table(
  mean = runif(params$n_iters, min = params$dur_asympt_mean_lb, max = params$dur_asympt_mean_ub),
  var = runif(params$n_iters, min = params$dur_asympt_var_lb, max = params$dur_asympt_var_ub)
)
dt_d_asympt[ , par1 := mean^2/var]
dt_d_asympt[ , par2 := var/mean]

#Sample sympt duration gamma parameters
dt_d_sympt <- data.table(
  mean = runif(params$n_iters, min = params$dur_sympt_mean_lb, max = params$dur_sympt_mean_ub),
  var = runif(params$n_iters, min = params$dur_sympt_var_lb, max = params$dur_sympt_var_ub)
)
dt_d_sympt[ , par1 := mean^2/var]
dt_d_sympt[ , par2 := var/mean]

#Sample presympt duration gamma parameters (truncation implemented in iteration forloop below)
dt_d_presympt <- data.table(
  mean = runif(params$n_iters, min = params$dur_presympt_mean_lb, max = params$dur_presympt_mean_ub),
  var = runif(params$n_iters, min = params$dur_presympt_var_lb, max = params$dur_presympt_var_ub),
  trunc_min = runif(params$n_iters, min = params$dur_presympt_truncmin_lb, max = params$dur_presympt_truncmin_ub),
  trunc_max = runif(params$n_iters, min = params$dur_presympt_truncmax_lb, max = params$dur_presympt_truncmax_ub)
)
dt_d_presympt[ , par1 := mean^2/var]
dt_d_presympt[ , par2 := var/mean]









d_sympt_100_iters <- data.table(
  iter=numeric(),
  dur_presympt=numeric(),
  dur_sympt=numeric(),
  dur_incubation=numeric(),
  t_infectious_start=numeric(),
  t_sympt_start=numeric(),
  t_recovery=numeric(),
  dur_infected_prequarantine=numeric()
)

d_asympt_100_iters <- data.table(
  iter=numeric(),
  dur_presympt=numeric(),
  dur_infectious=numeric(),
  dur_incubation=numeric(),
  t_infectious_start=numeric(),
  t_recovery=numeric(),
  dur_infected_prequarantine=numeric()
)

n_sympt = 100
n_asympt = 100

for (row in 1:300){
  #print(row)
  #[PRE]SYMPOMATIC INFECTION
  #Sample all durations using distribution parameters from outer loop
  dt_sympt <-
    data.table(
      dur_presympt = rtrunc(n_sympt, "gamma", dt_d_presympt[row, trunc_min], dt_d_presympt[row, trunc_max], 
                            shape=dt_d_presympt[row, par1], scale=dt_d_presympt[row, par2]),
      dur_sympt = rgamma(n_sympt, 
                         shape=dt_d_sympt[row, par1], scale = dt_d_sympt[row, par2]),
      dur_incubation = rlnorm(n_sympt, 
                              meanlog=dt_incubation_dists_lnorm[row, par1], 
                              sdlog = dt_incubation_dists_lnorm[row, par2])
    )
  
  #Calc times from durations. t=0 is time of infection
  dt_sympt[ , t_infectious_start := pmax(params$dur_latent_min, dur_incubation - dur_presympt)]
  dt_sympt[ , t_sympt_start := dur_incubation]
  dt_sympt[ , t_recovery := dur_incubation + dur_sympt]
  #Infection is at random point of presymptomatic phase on arrival
  #dt_sympt[ , dur_infected_prequarantine := runif(n_sympt)*t_sympt_start]
  #Infection is at random point of presymptomatic phase or first 24h of symptomatic phase on arrival
  dt_sympt[ , dur_infected_prequarantine := runif(n_sympt)*pmin(t_sympt_start+1, t_recovery)]
  
  
  
  #ASYMPTOMATIC INFECTION
  #Sample all durations
  dt_asympt <-
    data.table(
      dur_presympt = rtrunc(n_asympt, "gamma", dt_d_presympt[row, trunc_min], dt_d_presympt[row, trunc_max], 
                            shape=dt_d_presympt[row, par1], scale=dt_d_presympt[row, par2]),
      dur_infectious = rgamma(n_asympt, 
                              shape=dt_d_asympt[row, par1], scale = dt_d_asympt[row, par2]),
      dur_incubation = rlnorm(n_asympt, 
                              meanlog=dt_incubation_dists_lnorm[row, par1], 
                              sdlog = dt_incubation_dists_lnorm[row, par2])
    )
  
  #Calc times
  dt_asympt[ , t_infectious_start := pmax(params$dur_latent_min, dur_incubation - dur_presympt)]
  dt_asympt[ , t_recovery := t_infectious_start + dur_infectious]
  if(params$infection_timing == "rand_presympt" | params$infection_timing == "rand_incl_sympt"){
    dt_asympt[ , dur_infected_prequarantine := runif(n_asympt)*t_recovery]
  } else {
    dt_asympt[ , dur_infected_prequarantine := 0]
  }

  
  d_sympt_100_iters <- rbind(d_sympt_100_iters,
                             cbind(iter=row,
                                   dt_sympt))
  
  d_asympt_100_iters <- rbind(d_asympt_100_iters,
                             cbind(iter=row,
                                   dt_asympt))
}



ggplot()+
  geom_point(aes(x=dur_infected_prequarantine, y = t_infectious_start), alpha = 0.2, data=d_sympt_100_iters)+
  geom_abline(intercept = -3, slope = 1, color = "red")+
  geom_label(aes(x = 14, y = 3, label = "Detectable >=3 days before arrival"), color ="red", fontface="bold")+
  ggtitle("Traveller with symtomatic infection")+
  xlim(0,20)+ylim(0,20)+
  xlab("Days from infection to arrival")+
  ylab("Days from infection to infectiousness+detectable")

ggplot()+
  geom_point(aes(x=dur_infected_prequarantine, y = t_infectious_start), alpha = 0.2, data=d_asympt_100_iters)+
  geom_abline(intercept = -3, slope = 1, color = "red")+
  geom_label(aes(x = 14, y = 3, label = "Detectable >=3 days before arrival"), color ="red", fontface="bold")+
  ggtitle("Traveller with asymtomatic infection")+
  xlim(0,20)+ylim(0,20)+
  xlab("Days from infection to arrival")+
  ylab("Days from infection to infectiousness+detectable")

