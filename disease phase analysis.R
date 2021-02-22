## Code for analyzing and plotting disease phase over time under
# different assumptions of infection timing WRT arrival time.



library(data.table)
library(ggplot2)
library(scales)
theme_set(theme_bw())

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
  dur_presympt_shift = 0.5,
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
  dur_quarantine_endtest = c(2, 5, 7, 14),
  dur_quarantine_arrivetest = c(0, 2, 5, 7, 14),
  dur_quarantine_pretest = c(0, 2, 5, 7, 14),
  seed = 91,
  infection_timing = "rand_nosympt_24hr_before_arrival",
  rr_asympt = 0.49 #Relative risk for transmission from asympt vs. sympt person
)


asympt_infection_state_at_time <- function(dt_asympt, time_wrt_arrival){
  dt_asympt[, t_ref := dur_infected_prequarantine+time_wrt_arrival]
  dt_asympt[ , infection_state := ifelse(t_ref < 0, "Not yet infected",
                                         ifelse(t_ref < t_infectious_start, "Pre-infectious (undetectable)",
                                                ifelse(t_ref < t_recovery, "Asymptomatic-infectious",
                                                       "No longer infectious")))]
  
  x<-dt_asympt[, list("percent"=.N/nrow(dt_asympt)), by=infection_state]
  ids <- c("Not yet infected", "Pre-infectious (undetectable)", "Asymptomatic-infectious",
           "No longer infectious")
  id2 <- setdiff(ids, x$infection_state)
  x2 <- data.table(infection_state = id2, percent = 0)
  x <- rbind(x, x2)
  
  return(cbind(time_wrt_arrival = time_wrt_arrival, x))
}


sympt_infection_state_at_time <- function(dt_sympt, time_wrt_arrival){
  dt_sympt[, t_ref := dur_infected_prequarantine+time_wrt_arrival]
  dt_sympt[ , infection_state := ifelse(t_ref < 0, "Not yet infected",
                                        ifelse(t_ref < t_infectious_start, "Pre-infectious (undetectable)",
                                               ifelse(t_ref < t_sympt_start, "Presymptomatic-infectious",
                                                      ifelse(t_ref < t_recovery, "Symptomatic infectious",
                                                             "No longer infectious"))))]
  
  x<-dt_sympt[, list("percent"=.N/nrow(dt_sympt)), by=infection_state]
  ids <- c("Not yet infected", "Pre-infectious (undetectable)", "Presymptomatic-infectious", 
           "Symptomatic infectious", "No longer infectious")
  id2 <- setdiff(ids, x$infection_state)
  x2 <- data.table(infection_state = id2, percent = 0)
  x <- rbind(x, x2)
  
  return(cbind(time_wrt_arrival = time_wrt_arrival, x))
}















dgamma_shifted <- function(x, shape, scale, shift){
  return(dgamma(x=x, shape=shape, scale=scale))
}



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

#Sample presympt duration gamma parameters
dt_d_presympt <- data.table(
  mean = runif(params$n_iters, min = params$dur_presympt_mean_lb, max = params$dur_presympt_mean_ub),
  var = runif(params$n_iters, min = params$dur_presympt_var_lb, max = params$dur_presympt_var_ub),
  shift = params$dur_presympt_shift
)

dt_d_presympt[ , par1 := (mean-shift)^2/var]
dt_d_presympt[ , par2 := var/(mean-shift)]



p_dist_incubation <-
  ggplot() + 
  mapply(function(par1, par2){
    stat_function(fun=dlnorm, args = list(meanlog=par1, sdlog=par2))},
    par1 = dt_incubation_dists_lnorm[1:100, ]$par1,
    par2 = dt_incubation_dists_lnorm[1:100, ]$par2
  )+
  scale_x_continuous(limits = c(0, 13), breaks = seq(0,12,2))+
  #xlab("Days from infection to symptoms")+
  ggtitle("Incubation time\n(lognormal)")+
  xlab("Duration in days")+
  ylab("Density of 100 random distributions")

p_dist_presympt <-
  ggplot() + 
  mapply(function(par1, par2, shift){
    stat_function(fun=dgamma, args = list(shape=par1, scale=par2))},
    par1 = dt_d_presympt[1:100, ]$par1,
    par2 = dt_d_presympt[1:100, ]$par2
  )+
  scale_x_continuous(limits = c(-.5, 12.5), breaks = seq(-.5,12.5,2), labels = seq(0,12,2))+
  #xlab("Days pre-symptomatic\n(truncated gamma; subtracted from incubation)")+
  ylab("Density of 100 random distributions")+
  xlab("Duration in days")+
  ggtitle("Pre-symptomatic phase\n(gamma shifted 0.5)")

p_dist_sympt <- 
  ggplot() + 
  mapply(function(par1, par2){
    stat_function(fun=dgamma, args = list(shape=par1, scale=par2))},
    par1 = dt_d_sympt[1:100, ]$par1,
    par2 = dt_d_sympt[1:100, ]$par2
  )+
  scale_x_continuous(limits = c(0, 13), breaks = seq(0,12,2))+
  xlab("Duration in days")+
  ylab("Density of 100 random distributions")+
  ggtitle("Symptomatic phase\n(gamma)")

p_dist_asympt <-
  ggplot() + 
  mapply(function(par1, par2){
    stat_function(fun=dgamma, args = list(shape=par1, scale=par2))},
    par1 = dt_d_asympt[1:100, ]$par1,
    par2 = dt_d_asympt[1:100, ]$par2
  )+
  scale_x_continuous(limits = c(0, 13), breaks = seq(0,12,2))+
  xlab("Duration in days")+
  ylab("Density of 100 random distributions")+
  ggtitle("Asymptomatic infectious phase\n(gamma)")



ggsave("./manuscript/figs/fig_distribution_sampling.png",
       grid.arrange(p_dist_incubation,
                    p_dist_presympt,
                    p_dist_sympt,
                    p_dist_asympt
       ),
       width = 6.5, height = 6.5, units = "in")

n_sympt = 5000
n_asympt = 5000


#Table to store output
dt_state_wrt_arrival <-data.table(
  iter=numeric(),
  infection_timing = character(),
  infection_type = character(),
  time_wrt_arrival = numeric(),
  infection_state = character(),
  percent = numeric()
)

lst_infection_timing <- c("rand_incl_sympt",
                          "rand_nosympt_before_arrival",
                          "rand_nosympt_24hr_before_arrival",
                          "on_arrival")

for (row in 447:params$n_iters){
  #print(row)
  #[PRE]SYMPOMATIC INFECTION
  #Sample all durations using distribution parameters from outer loop
  dt_sympt <-
    data.table(
      dur_presympt = dt_d_presympt[row, shift]+rgamma(n_sympt, 
                                                      shape=dt_d_presympt[row, par1], scale = dt_d_presympt[row, par2]),
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
  
  
  
  #ASYMPTOMATIC INFECTION
  #Sample all durations
  dt_asympt <-
    data.table(
      dur_presympt = dt_d_presympt[row, shift]+rgamma(n_sympt, 
                                                      shape=dt_d_presympt[row, par1], scale = dt_d_presympt[row, par2]),
      dur_infectious = rgamma(n_asympt, 
                              shape=dt_d_asympt[row, par1], scale = dt_d_asympt[row, par2]),
      dur_incubation = rlnorm(n_asympt, 
                              meanlog=dt_incubation_dists_lnorm[row, par1], 
                              sdlog = dt_incubation_dists_lnorm[row, par2])
    )
  
  #Calc times
  dt_asympt[ , t_infectious_start := pmax(params$dur_latent_min, dur_incubation - dur_presympt)]
  dt_asympt[ , t_recovery := t_infectious_start + dur_infectious]
  
  
  for(infection_timing in lst_infection_timing){
    if(infection_timing == "rand_incl_sympt"){
      dt_sympt[ , dur_infected_prequarantine := runif(n_sympt)*t_recovery]
    } else if(infection_timing == "rand_nosympt_before_arrival"){
      dt_sympt[ , dur_infected_prequarantine := runif(n_sympt)*t_sympt_start]
    } else if(infection_timing == "rand_nosympt_24hr_before_arrival"){
      dt_sympt[ , dur_infected_prequarantine := runif(n_sympt)*pmin(t_sympt_start+1, t_recovery)]
    }else{
      dt_sympt[ , dur_infected_prequarantine := 0]
    }
    
    
    if(infection_timing == "on_arrival"){
      dt_asympt[ , dur_infected_prequarantine := 0]
    } else {
      dt_asympt[ , dur_infected_prequarantine := runif(n_asympt)*t_recovery]
    }
    
    
    ###Calculate at status at days
    for(time_wrt_arrival in seq(-3, 14, 1)){
      dt_state_wrt_arrival <- rbind(dt_state_wrt_arrival,
                                    cbind(iter=row,
                                          infection_timing=infection_timing,
                                          infection_type="Asymptomatic infection",
                                          asympt_infection_state_at_time(dt_asympt, time_wrt_arrival)))
      dt_state_wrt_arrival <- rbind(dt_state_wrt_arrival,
                                    cbind(iter=row,
                                          infection_timing=infection_timing,
                                          infection_type="Symptomatic infection",
                                          sympt_infection_state_at_time(dt_sympt, time_wrt_arrival)))
    }
  }
}



#dt_state_wrt_arrival_plt[time_wrt_arrival==-3 & infection_state %in% c("Not yet infected", "Pre-infectious (undetectable)"), sum(q0.5), by=c("infection_timing", "infection_type")]
temp<- dt_state_wrt_arrival[infection_state %in% c("Not yet infected", "Pre-infectious (undetectable)"), list(percent = sum(percent)), by=c("iter", "infection_timing", "infection_type", "time_wrt_arrival")]
temp2 <- temp[ , list(quantile = paste0("q", c(0.01, 0.5, 0.99)),
                                                         percent = quantile(percent, probs = c(0.01, 0.5, 0.99))), 
                                                  by=c("infection_timing","infection_type", "time_wrt_arrival")]

temp2[time_wrt_arrival == -3 & quantile == "q0.5"]
temp2[time_wrt_arrival == 0& quantile == "q0.5"]
temp2[time_wrt_arrival == -3, min(percent), by=c("quantile")]
temp2[time_wrt_arrival == 0, min(percent), by=c("quantile")]

infection_timing_states<-c("rand_incl_sympt",
                           "rand_nosympt_24hr_before_arrival",
                           "rand_nosympt_before_arrival",
                           "on_arrival")

lbls_infection_timing <- c("Traveling with\nsymptoms", "<24h of symptoms\non arrival", 
                           "No symptoms\nbefore arrival", "Infected on arrival")
names(lbls_infection_timing) <- infection_timing_states



dt_state_wrt_arrival_plt <- dt_state_wrt_arrival[ , list(quantile = paste0("q", c(0.01, 0.5, 0.99)),
                                                         percent = quantile(percent, probs = c(0.01, 0.5, 0.99))), 
                                                  by=c("infection_timing","infection_type", "time_wrt_arrival", "infection_state")]


dt_state_wrt_arrival_plt <- dcast(dt_state_wrt_arrival_plt, 
                                  infection_timing+infection_type+time_wrt_arrival+infection_state ~ quantile, 
                                  value.var = "percent")


dt_state_wrt_arrival_plt[,infection_state := factor(infection_state,
                                                   levels = c("Not yet infected", "Pre-infectious (undetectable)", "Presymptomatic-infectious", 
                                                              "Symptomatic infectious", "Asymptomatic-infectious","No longer infectious"))]

dt_state_wrt_arrival_plt[,infection_timing:=factor(infection_timing,
                                                   levels=infection_timing_states)]



ggplot(dt_state_wrt_arrival_plt)+
  geom_vline(xintercept = 0,  color="black")+
  geom_pointrange(aes(x=time_wrt_arrival, y = q0.5, ymin=q0.01, ymax=q0.99, color = infection_state), shape=20)+
  geom_line(aes(x=time_wrt_arrival, y = q0.5, color = infection_state))+
  facet_grid(cols = vars(infection_type), rows = vars(infection_timing),
             labeller = labeller(infection_timing = lbls_infection_timing))+
  scale_y_continuous(labels = percent_format())+
  scale_color_discrete(name="States")+
  scale_x_continuous(breaks = c(-3, 0, 5, 10, 14), minor_breaks = -3:14)+
  labs(x = "Time with respect to arrival (in days)",
       y = "Percent of infected travelers in state (median and 98% credible interval)")+
    theme(legend.position = "bottom")

ggsave("./manuscript/figs/disease_state_distributions.png",
       width=6.5, height = 7, unit = "in")
  
