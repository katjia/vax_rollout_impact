################################################################################
# @Project - Defining and Estimating Outcomes Directly Averted by a Vaccination 
# Program when Rollout Occurs Over Time
# @Description - A snippet to examine why HD estimator underestimates vaccine-av
# ertible deaths
# @Author - Katherine Jia
################################################################################
# CONFIGURATION ================================================================
SIMULATION_MODE = 3
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
SIRD_mod <- set_SIRD_parms_and_run_model(scenario = 3) # covid

# RUN MODEL ====================================================================
d <- get_Delta_Y(SIRD_mod) %>%
  select(vax_id, kplus1 = kplus1_first_day, d_Y = kplus1_prob_death)

# avertible deaths (causal estimand)
3e5*((0.3*d$d_Y[3]+0.5*d$d_Y[5])-0.8*d$d_Y[1]+0.3*(d$d_Y[4]-d$d_Y[2])+0.5*(d$d_Y[6]-d$d_Y[2]))
# avertible deaths (hazard difference estimator)
3e5*((0.3*d$d_Y[3]+0.5*d$d_Y[5])-0.8*d$d_Y[1]+0.5*d$d_Y[6]-0.5*(0.5*(1-d$d_Y[5])/(0.2*(1-d$d_Y[1])+0.3*(1-d$d_Y[3]))*(0.2*d$d_Y[2]+0.3*d$d_Y[4])/0.5))

# bias factors
(bf1 <- 0.3*(d$d_Y[4]-d$d_Y[2])) # plus factor
(bf2 <- 0.5*d$d_Y[2] - 0.5*(0.5*(1-d$d_Y[5])/(0.2*(1-d$d_Y[1])+0.3*(1-d$d_Y[3]))*(0.2*d$d_Y[2]+0.3*d$d_Y[4])/0.5)) # minus factor
