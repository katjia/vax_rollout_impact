################################################################################
# @Project - Defining and Estimating Outcomes Directly Averted by a Vaccination 
# Program when Rollout Occurs Over Time
# @Description - Estimating rho and averted outcomes based on obs data aggregate
# d by vax status
# @Author - Katherine Jia
################################################################################
# CONFIGURATION ================================================================
SIMULATION_MODE = 1
strategy = 1

# Initialize configuration
config <- configure_simulation(SIMULATION_MODE, strategy)
N_t <- config$N_t
mod_vax_id <- config$mod_vax_id
RHO_vec <- config$RHO_vec
mod_time_id <- config$mod_time_id
vax_length <- config$vax_length
SIRD_mod_dims <- config$SIRD_mod_dims
SIRD_parmdims <- config$SIRD_parmdims
SIRD_parms <- config$SIRD_parms
scenarios <- config$scenarios
scenario_labels <- config$scenario_labels

# Run model
SCENARIO = 1
SIRD_mod <- set_SIRD_parms_and_run_model(scenario = SCENARIO) 

# ESTIMATING RHO ===============================================================
survival_hazard_by_vax_status <- get_survival_hazard_by_vax_status(SIRD_mod) 

survival_hazard_by_vax_status <- survival_hazard_by_vax_status %>%
  pivot_wider(values_from = c(aggr_survivors_death, aggr_survivors_inf, aggr_hazards_death, aggr_hazards_inf), names_from = status) 

# TABLE S11
obs_dataset <- survival_hazard_by_vax_status %>%
  select(k_first_day,
         aggr_survivors_death_Vax,
         aggr_hazards_death_Vax,
         aggr_survivors_death_Unvax,
         aggr_hazards_death_Unvax,
         kplus1_first_day)

# OUTPUT TABLE S11
# use 9 digits 
(tableS11_9digits <- round(obs_dataset, 9))
if(SIMULATION_MODE==1 & strategy==1 & SCENARIO == 1) write.csv(tableS11_9digits, "~/Documents/GitHub/direct_impact_estimands/4_tables/tableS11.csv")

# add survivors and hazard in the next interval as a separate column
copy_obs_dataset <- obs_dataset %>%
  select(k_first_day,
         aggr_survivors_death_Vax_kplus1 = aggr_survivors_death_Vax,
         aggr_hazards_death_Vax_kplus1 = aggr_hazards_death_Vax,
         aggr_hazards_death_Unvax_kplus1 = aggr_hazards_death_Unvax)

survival_hazard_by_vax_status_combined <- obs_dataset %>%
  left_join(copy_obs_dataset, 
            by=c("kplus1_first_day"="k_first_day"))

# gather the components of eq S17
survival_hazard_by_vax_status_combined <- survival_hazard_by_vax_status_combined %>%
  mutate(aggr_survival_death_Unvax = cumprod((1-aggr_hazards_death_Unvax))) %>%
  select(kplus1_first_day,
         aggr_survivors_death_Vax_kplus1,
         aggr_survivors_death_Vax,
         aggr_hazards_death_Vax,
         aggr_survival_death_Unvax)

# estimate rho using equation S17
survival_hazard_by_vax_status_combined <- survival_hazard_by_vax_status_combined  %>%
  mutate(kplus1_rho=(aggr_survivors_death_Vax_kplus1 - aggr_survivors_death_Vax * (1-aggr_hazards_death_Vax))/(3e5*aggr_survival_death_Unvax)) 

estimated_rho <- survival_hazard_by_vax_status_combined %>%
  select(kplus1_first_day, kplus1_rho)

# add rho of the first interval (which we know) 
# Note: since we add one more row, k starts from -1 to q, instead of 0 to q
if(strategy==1) estimated_rho <- rbind(c(1, 0.2), estimated_rho)
if(strategy==2) estimated_rho <- rbind(c(1, 0.01), estimated_rho)

# calculate rho of the last interval (as we know rho of the first interval)
estimated_rho[nrow(estimated_rho), "kplus1_rho"] <- 1-sum(estimated_rho$kplus1_rho, na.rm = T)

# ESTIMATING AVERTED DEATHS USING OBSERVATIONAL DATA ===========================
obs_dataset %>%
  left_join(estimated_rho %>% 
            select(k_first_day = kplus1_first_day,
                   k_rho = kplus1_rho), 
            by="k_first_day") %>%
  mutate(odds = cumsum(k_rho)/(1-cumsum(k_rho)),
         k_averted_deaths = odds * aggr_survivors_death_Unvax * aggr_hazards_death_Unvax - aggr_survivors_death_Vax * aggr_hazards_death_Vax) %>%
  pull(k_averted_deaths) %>%
  sum()

# ESTIMATING AVERTED DEATHS USING RCT ==========================================
cum_Y <- get_Y(SIRD_mod)

tableS12_9digits <- cum_Y %>% 
  mutate_if(is.numeric, round, digits=9)

if(SIMULATION_MODE == 1 & strategy == 1 & SCENARIO == 1) write.csv(tableS12_9digits, "~/Documents/GitHub/direct_impact_estimands/4_tables/tableS12.csv")

unbiased_estimator <- causal_CID(cum_Y) 
unbiased_estimator["averted_deaths"]
