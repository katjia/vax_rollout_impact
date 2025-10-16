################################################################################
# @Project - Defining and Estimating Outcomes Directly Averted by a Vaccination 
# Program when Rollout Occurs Over Time
# @Description - Simulations + bias calculation
# @Author - Katherine Jia
################################################################################

# Load required libraries
library(tidyverse)
library(odin)
library(patchwork)
library(stringr)
library(latex2exp)
library(ggplot2)
library(ggh4x)
options(digits = 15)

# CONFIGURATION SECTION ========================================================

# Simulation modes:
# 1 = Base simulations (varying mu, fixed beta)
# 2 = Varying beta simulations (varying beta, fixed mu) 
# 3 = Realistic parameters (influenza, measles, SARS-CoV-2)
SIMULATION_MODE <- 1

# Analysis mode:
# TRUE = Calculate and output absolute/percentage differences
CALCULATE_DIFFERENCES <- TRUE

# MODEL PARAMETERS AND INDICES =================================================
# beta: The number of effective contacts made by a typical infectious individual per day
# mu: Probability of death due to infection 
# theta: 1 – vaccine efficacy against infection (VE_infection)/100%
# kappa: 1 – vaccine efficacy against death given infection (VE_(death|infection))/100%
# gamma: Recovery rate per day

# Model components
comp_id <- c("S", "I", "R", "D", "Cum_inf")
dt <- 1

# Define vaccination strategies
# Strategy 1 (default): Vaccinate at two time points 
#   - Day 0: 20%, Day 60: 30%, remaining 50% unvaccinated
#   - Note: Vaccination takes effect one day after administration 
#           (e.g., vaccination on Day 60 is effective starting Day 61)
#
# Strategy 2: Vaccinate at 99 weekly time points
#   - 1% vaccinated at the beginning of each week (Weeks 0–98), with 1% remaining unvaccinated

strategy <- 1 
# strategy <- 2

# HELPER FUNCTIONS =============================================================

# Function to create vaccine efficacy matrix
create_ve_matrix <- function(vax_length, N_t, mod_vax_id, ve_value) {
  m2 <- matrix(1, vax_length, N_t + 1)
  for(j in 1:length(mod_vax_id)){
    m2[j, (mod_vax_id[j]):(N_t+1)] <- ve_value
  }
  return(m2)
}

# CONFIGURATION SETUP ==========================================================

# Set up simulation parameters based on mode and strategy
configure_simulation <- function(mode, strategy) {
  if(mode == 3) { # Realistic parameters
    N_t <- 120 
    mod_vax_id <- c(0, 60, N_t) + 1   
    RHO_vec <- c(0.2, 0.3, 0.5)   
    scenarios <- 1:3  # 3 scenarios correspond to 3 pathogens
    scenario_labels <- list(
      "1" = "Influenza",
      "2" = "Measles", 
      "3" = "SARS-CoV-2 WT"
    )
  } else {
    # Base and varying beta simulations
    if(strategy == 1) {
      N_t <- 120
      mod_vax_id <- c(0, 60, N_t) + 1   
      RHO_vec <- c(0.2, 0.3, 0.5)   
    } else {
      N_t <- 693
      mod_vax_id <- seq(0, N_t, 7) + 1
      RHO_vec <- rep(0.01, 100)
    }
    scenarios <- 1:9  # 9 scenarios for base/varying beta
    
    # Set scenario labels based on simulation mode
    if(mode == 2) {
      # Mode 2: Varying beta simulations 
      scenario_labels <- list(
        "1"=TeX('I: $\\beta$=0.15; VE$_{inf}$=90%; VE$_{death}$=0%'), 
        "2"=TeX('II: $\\beta$=0.2; VE$_{inf}$=90%; VE$_{death}$=0%'),
        "3"=TeX('III: $\\beta$=0.25; VE$_{inf}$=90%; VE$_{death}$=0%'),
        "4"=TeX('IV: $\\beta$=0.15; VE$_{inf}$=0; VE$_{death}$=90%'),
        "5"=TeX('V: $\\beta$=0.2; VE$_{inf}$=0; VE$_{death}$=90%'),
        "6"=TeX('VI: $\\beta$=0.25; VE$_{inf}$=0; VE$_{death}$=90%'),
        "7"=TeX('VII: $\\beta$=0.15; VE$_{inf}$=90%; VE$_{death}$=90%'),
        "8"=TeX('VIII: $\\beta$=0.2; VE$_{inf}$=90%; VE$_{death}$=90%'),
        "9"=TeX('IX: $\\beta$=0.25; VE$_{inf}$=90%; VE$_{death}$=90%')
        )
    } else {
      # Mode 1: Base simulations - use IFR values
      scenario_labels <- list(
        "1"=TeX('1: IFR=1%; VE$_{inf}$=90%; VE$_{death}$=0%'),
        "2"=TeX('2: IFR=10%; VE$_{inf}$=90%; VE$_{death}$=0%'),
        "3"=TeX('3: IFR=100%; VE$_{inf}$=90%; VE$_{death}$=0%'),
        "4"=TeX('4: IFR=1%; VE$_{inf}$=0; VE$_{death}$=90%'),
        "5"=TeX('5: IFR=10%; VE$_{inf}$=0; VE$_{death}$=90%'),
        "6"=TeX('6: IFR=100%; VE$_{inf}$=0; VE$_{death}$=90%'),
        "7"=TeX('7: IFR=1%; VE$_{inf}$=90%; VE$_{death}$=90%'),
        "8"=TeX('8: IFR=10%; VE$_{inf}$=90%; VE$_{death}$=90%'),
        "9"=TeX('9: IFR=100%; VE$_{inf}$=90%; VE$_{death}$=90%')
      )
    }
  }
  
  mod_time_id <- seq(1, N_t + 1, dt)
  vax_length <- length(mod_vax_id)
  
  # Create model dimensions
  SIRD_mod_dims <- list(comp_id = comp_id, 
                        vax_id = mod_vax_id, 
                        time_id = mod_time_id)
  
  # Create parameter matrix
  SIRD_parmdims <- list(from = comp_id, 
                        to = comp_id, 
                        vax_id = mod_vax_id, 
                        time = mod_time_id)
  SIRD_parms <- array(0, lengths(SIRD_parmdims), SIRD_parmdims)
  
  # Create results matrix
  res_dims <- list(scenarios = scenarios, 
                   type = c("E[ hazard \n difference \n estimator ]", "causal \n estimand"),
                   measures = c("averted_deaths", "avertible_deaths", "averted_infections", "avertible_infections"))
  res <- array(0, lengths(res_dims), res_dims)
  
  return(list(
    N_t = N_t, mod_vax_id = mod_vax_id, RHO_vec = RHO_vec,
    mod_time_id = mod_time_id, vax_length = vax_length,
    SIRD_mod_dims = SIRD_mod_dims, SIRD_parmdims = SIRD_parmdims,
    SIRD_parms = SIRD_parms, res = res, scenarios = scenarios,
    scenario_labels = scenario_labels
  ))
}

# Initialize configuration (GLOBAL PARAMETERS)
config <- configure_simulation(SIMULATION_MODE, strategy)
N_t <- config$N_t
mod_vax_id <- config$mod_vax_id
RHO_vec <- config$RHO_vec
mod_time_id <- config$mod_time_id
vax_length <- config$vax_length
SIRD_mod_dims <- config$SIRD_mod_dims
SIRD_parmdims <- config$SIRD_parmdims
SIRD_parms <- config$SIRD_parms
res <- config$res
scenarios <- config$scenarios
scenario_labels <- config$scenario_labels

# Initialize placeholders for outputs
prob_outcome_all <- prob_outcome_vax_status_all <- dynamics_scenario_all <- hazards_outcome_by_vax_status_all <- hazards_outcome_by_vax_time_all <- NULL

# MAIN SIMULATION FUNCTION =====================================================
set_SIRD_parms_and_run_model <- function(scenario=1){
  
  # Define base parameters
  gamma <- 0.07 # may be overridden for SIMULATION_MODE==3
  beta <- 0.25  # Default beta, may be overridden for SIMULATION_MODE == 2
  
  # Set parameters based on simulation mode
  if(SIMULATION_MODE == 1) {
    # Base simulations: varying mu, fixed beta
    
    # Set vaccine efficacies based on scenario groups
    if(scenario %in% c(1, 2, 3)){
      theta <- create_ve_matrix(vax_length, N_t, mod_vax_id, 0.1)  # VE inf = 90%
      kappa <- create_ve_matrix(vax_length, N_t, mod_vax_id, 1)    # VE death = 0%
    } else if(scenario %in% c(4, 5, 6)){
      theta <- create_ve_matrix(vax_length, N_t, mod_vax_id, 1)    # VE inf = 0%
      kappa <- create_ve_matrix(vax_length, N_t, mod_vax_id, 0.1)  # VE death = 90%
    } else if(scenario %in% c(7, 8, 9)){
      theta <- create_ve_matrix(vax_length, N_t, mod_vax_id, 0.1)  # VE inf = 90%
      kappa <- create_ve_matrix(vax_length, N_t, mod_vax_id, 0.1)  # VE death = 90%
    }
    
    # Set mortality rates
    if(scenario %in% c(1, 4, 7)) mu <- 0.01
    else if(scenario %in% c(2, 5, 8)) mu <- 0.1
    else if(scenario %in% c(3, 6, 9)) mu <- 1
    
  } else if(SIMULATION_MODE == 2) {
    # Varying beta simulations: varying beta, fixed mu
    mu <- 0.1  # Fixed mu
    
    # Set vaccine efficacies (same as base)
    if(scenario %in% c(1, 2, 3)){
      theta <- create_ve_matrix(vax_length, N_t, mod_vax_id, 0.1)
      kappa <- create_ve_matrix(vax_length, N_t, mod_vax_id, 1)
    } else if(scenario %in% c(4, 5, 6)){
      theta <- create_ve_matrix(vax_length, N_t, mod_vax_id, 1)
      kappa <- create_ve_matrix(vax_length, N_t, mod_vax_id, 0.1)
    } else if(scenario %in% c(7, 8, 9)){
      theta <- create_ve_matrix(vax_length, N_t, mod_vax_id, 0.1)
      kappa <- create_ve_matrix(vax_length, N_t, mod_vax_id, 0.1)
    }
    
    # Set varying beta values
    if(scenario %in% c(1, 4, 7)) beta <- 0.15
    else if(scenario %in% c(2, 5, 8)) beta <- 0.2
    else if(scenario %in% c(3, 6, 9)) beta <- 0.25
    
  } else if(SIMULATION_MODE == 3) {
    # Realistic pathogen simulations
    
    if(scenario == 1) { # Influenza
      beta <- round(1.95 * 0.21, 2)
      theta <- create_ve_matrix(vax_length, N_t, mod_vax_id, 0.66)  # VE inf = 34%
      kappa <- create_ve_matrix(vax_length, N_t, mod_vax_id, 0.69)  # VE death = 31%
      mu <- 0.03
      gamma <- 0.21
    } else if(scenario == 2) { # Measles
      beta <- round(18 * 0.07, 2)
      theta <- create_ve_matrix(vax_length, N_t, mod_vax_id, 0.05)  # VE inf = 95%
      kappa <- create_ve_matrix(vax_length, N_t, mod_vax_id, 0)     # VE death = 100%
      mu <- 0.013
    } else if(scenario == 3) { # SARS-CoV-2 WT
      beta <- round(2.2 * 0.07,2)
      theta <- create_ve_matrix(vax_length, N_t, mod_vax_id, 0.05)  # VE inf = 95%
      kappa <- create_ve_matrix(vax_length, N_t, mod_vax_id, 0.04)  # VE death = 96%
      mu <- round(10.73/10000, 4)
    }
  }
  
  # Build transition rates
  S_to_I <- array(as.vector(theta)*beta, lengths(SIRD_mod_dims[c(2,3)]), SIRD_mod_dims[c(2,3)])
  I_to_R <- array((1 - as.vector(kappa) * mu) * gamma, lengths(SIRD_mod_dims[c(2,3)]), SIRD_mod_dims[c(2,3)])
  I_to_D <- array(as.vector(kappa) * mu * gamma, lengths(SIRD_mod_dims[c(2,3)]), SIRD_mod_dims[c(2,3)])
  
  # Initialize parameter arrays
  for (t in 1:length(mod_time_id)) {
    for (j in 1:length(mod_vax_id)) {
      SIRD_parms[1,1,j,t] <- (1 - S_to_I[j,t]) 
      SIRD_parms[1,2,j,t] <- S_to_I[j,t] 
      SIRD_parms[2,2,j,t] <- (1 - I_to_R[j,t] - I_to_D[j,t])
      SIRD_parms[2,3,j,t] <- I_to_R[j,t]    
      SIRD_parms[2,4,j,t] <- I_to_D[j,t]    
      SIRD_parms[3,3,j,t] <- 1
      SIRD_parms[4,4,j,t] <- 1 
      SIRD_parms[1,5,j,t] <- S_to_I[j,t] 
      SIRD_parms[5,5,j,t] <- 1 
    }
  }
  
  # Initialize model
  SIRD_mod <- array(0,lengths(SIRD_mod_dims),SIRD_mod_dims)  
  SIRD_mod["S",,1] <- 3e5 * (1-0.001) * RHO_vec
  SIRD_mod["I",,1] <- 3e5 * 0.001 * RHO_vec
  SIRD_mod["Cum_inf",,1] <- 3e5 * 0.001 * RHO_vec
  
  # Run simulation
  for(t in 1:(length(mod_time_id)-1)){
    
      # Update force of infection dynamically
      lambda_t <- beta * sum(SIRD_mod["I",,t])/sum(SIRD_mod[c("S","I","R"),,t])
      S_to_I[,t] <- theta[,t]*lambda_t
    
    for(j in 1:(length(mod_vax_id))){
      # Update transmission matrix since lambda is dynamics
        SIRD_parms[1,1,j,t] <- (1 - S_to_I[j,t]) 
        SIRD_parms[1,2,j,t] <- S_to_I[j,t] 
        SIRD_parms[1,5,j,t] <- S_to_I[j,t] 
      
      # Update compartments
      SIRD_mod[,j,(t+1)] <- t(SIRD_parms[,,j,t]) %*%  SIRD_mod[,j,t]
    }
  }
  
  return(SIRD_mod)
}

# <III> GET OUTPUTS ------------------------------------------------------------
# 3.1 A function to get disease dynamics for plotting --------------------------
get_dynamics <- function(SIRD_mod){
  dynamics_scenario <- apply(SIRD_mod, MARGIN = c(1,3), function(x) x / (3e5 * RHO_vec)) %>%
    as.data.frame.table() 
}

# 3.2 Get cumulative incidence -------------------------------------------------
get_Y <- function(SIRD_mod){
  SIRD_mod %>%
    as.data.frame.table() %>%
    pivot_wider(names_from = comp_id, values_from = Freq) %>%
    group_by(vax_id) %>%
    summarise(
      cum_death = last(D),
      cum_inf = last(Cum_inf)) %>%
    mutate(prob_vax = RHO_vec) %>%
    mutate(prob_death = cum_death/(3e5*prob_vax),
           prob_inf = cum_inf/(3e5*prob_vax)) %>%
    select(vax_id, prob_vax, prob_death, prob_inf)
}

# 3.3 Get period incidence -----------------------------------------------------
get_Delta_Y <- function(SIRD_mod){
  SIRD_mod %>%
    as.data.frame.table() %>%
    as.data.frame() %>%
    filter(time_id %in% mod_vax_id) %>%
    pivot_wider(names_from = comp_id, values_from = Freq) %>%
    mutate(prob_vax = rep(RHO_vec, length(mod_vax_id))) %>%
    group_by(vax_id, prob_vax) %>%
    summarise(
      kplus1_diff_death = diff(D),
      kplus1_diff_inf = diff(Cum_inf)) %>%
    mutate(kplus1_first_day=mod_vax_id[-1],
           k_first_day = mod_vax_id[-length(mod_vax_id)]) %>%                   # `kplus1_first_day`: the first DAY of the NEXT interval
    select(vax_id, prob_vax, kplus1_first_day, k_first_day, kplus1_diff_death, kplus1_diff_inf) %>%             
    mutate(kplus1_prob_death = kplus1_diff_death/(3e5*prob_vax),
           kplus1_prob_inf = kplus1_diff_inf/(3e5*prob_vax)) %>%
    select(vax_id, kplus1_first_day, k_first_day,  prob_vax, kplus1_diff_death, kplus1_prob_death, kplus1_diff_inf, kplus1_prob_inf) %>%
    as.data.frame()
}

# 3.4 Get hazards and survival aggregated by vax status ------------------------
get_survival_hazard_by_vax_status <- function(SIRD_mod){

disc_outcome <- get_Delta_Y(SIRD_mod)

# calculating hazards using # survivors as denom
disc_survivors <- SIRD_mod %>%
  as.data.frame.table() %>%
  as.data.frame() %>%
  filter(time_id %in% mod_vax_id) %>%
  pivot_wider(names_from = comp_id, values_from = Freq) %>%
  mutate(prob_vax = rep(RHO_vec, length(mod_vax_id))) %>%
  mutate(k_survivors_death = 3e5*prob_vax-D,
         k_survivors_inf = 3e5*prob_vax-Cum_inf) %>%
  select(vax_id, time_id, k_survivors_death, k_survivors_inf)                   # time_id is essentially k_first_day

disc_survivors$time_id <- as.numeric(disc_survivors$time_id)

# merge incident cases with survival (survival is one interval behind outcomes)
disc_incident_cases_surv <- disc_outcome %>%
  left_join(disc_survivors, by= c("vax_id", "k_first_day" = "time_id"))

disc_incident_cases_surv %>%
  mutate(status = case_when(as.numeric(as.character(vax_id)) < kplus1_first_day ~ "Vax",
                            as.numeric(as.character(vax_id)) >= kplus1_first_day ~ "Unvax")) %>%
  group_by(kplus1_first_day, k_first_day, status) %>%
  summarise(aggr_hazards_death = sum(kplus1_diff_death)/sum(k_survivors_death),
            aggr_hazards_inf = sum(kplus1_diff_inf)/sum(k_survivors_inf),
            aggr_survivors_death = sum(k_survivors_death),
            aggr_survivors_inf = sum(k_survivors_inf)) %>%
  select(kplus1_first_day, k_first_day, status, aggr_survivors_death, aggr_survivors_inf, aggr_hazards_death, aggr_hazards_inf) %>%
  ungroup()

}

# <IV> CALCULATE AVERTED & AVERTIBLE OUTCOMES ----------------------------------
# 4.1 Calculate CID estimand using cumulative incidence difference
causal_CID <- function(cum_Y){
  cum_Y %>%
    group_by() %>%
    mutate(avertible_deaths_diff=prob_death-first(prob_death),
           avertible_deaths = 3e5 * prob_vax * avertible_deaths_diff,
           averted_deaths_diff = last(prob_death)-prob_death,
           averted_deaths = 3e5 * prob_vax * averted_deaths_diff,
           avertible_infections_diff=prob_inf-first(prob_inf),
           avertible_infections = 3e5 * prob_vax * avertible_infections_diff,
           averted_infections_diff = last(prob_inf)-prob_inf,
           averted_infections = 3e5 * prob_vax * averted_infections_diff) %>%
    select(averted_deaths, avertible_deaths, averted_infections, avertible_infections) %>%
    colSums()
}

# 4.2 Calculate HD estimator 
naive_method <- function(SIRD_mod){
  
  survival_hazard_by_vax_status <- get_survival_hazard_by_vax_status(SIRD_mod)
  
  survival_hazard_by_vax_status %>%
    pivot_wider(values_from = c(aggr_survivors_death, aggr_survivors_inf, aggr_hazards_death, aggr_hazards_inf), names_from = status) %>%
    mutate(averted_deaths = aggr_survivors_death_Vax*(aggr_hazards_death_Unvax-aggr_hazards_death_Vax),
           averted_infections = aggr_survivors_inf_Vax*(aggr_hazards_inf_Unvax-aggr_hazards_inf_Vax),
           avertible_deaths = aggr_survivors_death_Unvax*(aggr_hazards_death_Unvax-aggr_hazards_death_Vax),
           avertible_infections = aggr_survivors_inf_Unvax*(aggr_hazards_inf_Unvax-aggr_hazards_inf_Vax)) %>%
    select(averted_deaths, avertible_deaths, averted_infections, avertible_infections) %>%
    colSums() 
  
}

# MAIN SIMULATION EXECUTION ===================================================

# Function to create scenario labels based on mode
create_labels <- function(mode) {
  if(mode == 3) {
    list("1" = "Influenza", "2" = "Measles", "3" = "SARS-CoV-2 WT")
  } else if(mode == 2) {
    # Mode 2: Varying beta simulations - use beta values instead of IFR
    list(
      "1"=TeX('I: $\\beta$=0.15; VE$_{inf}$=90%; VE$_{death}$=0%'), 
      "2"=TeX('II: $\\beta$=0.2; VE$_{inf}$=90%; VE$_{death}$=0%'),
      "3"=TeX('III: $\\beta$=0.25; VE$_{inf}$=90%; VE$_{death}$=0%'),
      "4"=TeX('IV: $\\beta$=0.15; VE$_{inf}$=0; VE$_{death}$=90%'),
      "5"=TeX('V: $\\beta$=0.2; VE$_{inf}$=0; VE$_{death}$=90%'),
      "6"=TeX('VI: $\\beta$=0.25; VE$_{inf}$=0; VE$_{death}$=90%'),
      "7"=TeX('VII: $\\beta$=0.15; VE$_{inf}$=90%; VE$_{death}$=90%'),
      "8"=TeX('VIII: $\\beta$=0.2; VE$_{inf}$=90%; VE$_{death}$=90%'),
      "9"=TeX('IX: $\\beta$=0.25; VE$_{inf}$=90%; VE$_{death}$=90%')
    )
  } else {
    # Mode 1: Base simulations - use IFR values
    list(
      "1"=TeX('1: IFR=1%; VE$_{inf}$=90%; VE$_{death}$=0%'), 
      "2"=TeX('2: IFR=10%; VE$_{inf}$=90%; VE$_{death}$=0%'),
      "3"=TeX('3: IFR=100%; VE$_{inf}$=90%; VE$_{death}$=0%'),
      "4"=TeX('4: IFR=1%; VE$_{inf}$=0; VE$_{death}$=90%'),
      "5"=TeX('5: IFR=10%; VE$_{inf}$=0; VE$_{death}$=90%'),
      "6"=TeX('6: IFR=100%; VE$_{inf}$=0; VE$_{death}$=90%'),
      "7"=TeX('7: IFR=1%; VE$_{inf}$=90%; VE$_{death}$=90%'),
      "8"=TeX('8: IFR=10%; VE$_{inf}$=90%; VE$_{death}$=90%'),
      "9"=TeX('9: IFR=100%; VE$_{inf}$=90%; VE$_{death}$=90%')
    )
  }
}

# Run simulations
cat("Running simulations in mode", SIMULATION_MODE, "...n")
for(scenario in scenarios){
  cat("Processing scenario", scenario, "\n")
  
  SIRD_mod <- set_SIRD_parms_and_run_model(scenario = scenario)
  cum_Y <- get_Y(SIRD_mod)
  Delta_Y <- get_Delta_Y(SIRD_mod)
  causal_method_res <- causal_CID(cum_Y)
  naive_method_res <- naive_method(SIRD_mod)
  
  res[scenario, 1,] <- round(naive_method_res, 9)
  res[scenario, 2,] <- round(causal_method_res, 9)
  
  # Store dynamics for plotting
  dynamics_scenario <- get_dynamics(SIRD_mod)
  dynamics_scenario_all <- rbind(dynamics_scenario %>% mutate(scenarios=scenario), dynamics_scenario_all)
}

# Create labels for current mode
arr <- create_labels(SIMULATION_MODE)
mylabel <- function(val) { return(lapply(val, function(x) arr[x])) }

# CALCULATE DIFFERENCES (from abs_pct_diff.R) ==================================
if(CALCULATE_DIFFERENCES) {
  cat("Calculating absolute and percentage differences...n")
  
  # Calculate absolute and percentage differences
  # Note (1): Outcome measures are rounded to 9 decimal places to avoid displaying 
  # uninterpretable, extremely small non-zero values (below e^-9) for scenarios 4–6,
  # here averted or avertible infections are practically zero.
  # Note (2): Negative sign is to set the causal estimand as the ref group
  abs_diff <- - apply(res, c(1,3), diff)
  pct_diff <- abs_diff/res[,2,] * 100
  
  # Display the results (rounded to 2 d.p.)
  abs_diff <- t(round(abs_diff,2))
  pct_diff <- t(round(pct_diff,2))
  
  # Create output tables
  averted_outcomes <- rbind(paste0(abs_diff[3,], paste0("\n (", pct_diff[3,], "%)")),
                            paste0(abs_diff[1,], paste0("\n (", pct_diff[1,], "%)")))
  rownames(averted_outcomes) <- c("averted_infections", "averted_deaths")
    
  avertible_outcomes <- rbind(paste0(abs_diff[4,], paste0("\n (", pct_diff[4,], "%)")),
                              paste0(abs_diff[2,], paste0("\n (", pct_diff[2,], "%)")))
  rownames(avertible_outcomes) <- c("avertible_infections","avertible_deaths")
  
  # Output to CSV files
  mode_suffix <- switch(as.character(SIMULATION_MODE),
                       "1" = "base",
                       "2" = "varying_beta", 
                       "3" = "realistic")
  
  write.csv(averted_outcomes, 
            paste0("3_tables/averted_outcomes_", mode_suffix, "_strategy", strategy, ".csv"), 
            row.names = FALSE)
  if(SIMULATION_MODE!=2){
    write.csv(avertible_outcomes, 
            paste0("3_tables/avertible_outcomes_", mode_suffix, "_strategy", strategy, ".csv"), 
            row.names = FALSE)
  }
  
  cat("Difference tables saved to 3_tables/n")
  
  # CREATE HEATMAPS ===========================================================
  cat("Creating heatmaps for abs_diff and pct_diff...n")
  
  # Prepare data for heatmaps
  # abs_diff and pct_diff are currently transposed matrices with dimensions:
  # Rows: averted_infections, avertible_infections, averted_deaths, avertible_deaths (indices 3,4,1,2 respectively)
  # Columns: scenarios 1-9
  
  # Function to reshape data for heatmap (3x3 grid)
  prepare_heatmap_data <- function(data_matrix, row_index, outcome_type, sim_mode) {
    # Extract the relevant row
    values <- data_matrix[row_index, ]
    
    # Reshape into 3x3 matrix
    # Rows represent IFR/beta levels, Columns represent VE combinations
    # Scenarios 1,2,3 = first VE combination (VE_inf=90%, VE_death=0%)
    # Scenarios 4,5,6 = second VE combination (VE_inf=0%, VE_death=90%)
    # Scenarios 7,8,9 = third VE combination (VE_inf=90%, VE_death=90%)
    
    heatmap_matrix <- matrix(c(values[1], values[2], values[3],
                               values[4], values[5], values[6],
                               values[7], values[8], values[9]),
                             nrow = 3, ncol = 3, byrow = TRUE)
    
    # Create data frame for ggplot
    df <- data.frame(
      Row = rep(1:3, each = 3),
      Col = rep(1:3, times = 3),
      Value = as.vector(t(heatmap_matrix))
    )
    
    # Set row and column labels based on simulation mode
    if(sim_mode == 1) {
      col_labs <- c("IFR=1%", "IFR=10%", "IFR=100%")
    } else if(sim_mode == 2) {
      col_labs <- c("β=0.15", "β=0.2", "β=0.25")
    }
    
    row_labs <- c("VE inf=90%, \n  VE death=0",
                  "VE inf=0, \n  VE death=90%",
                  "VE inf=90%, \n  VE death=90%")
    
    df$Row_Label <- factor(df$Row, levels = 1:3, labels = row_labs)
    df$Col_Label <- factor(df$Col, levels = 1:3, labels = col_labs)
    
    return(df)
  }
  
  # Create heatmaps for each outcome type
  create_heatmap <- function(df, value_type) {
    
    if(SIMULATION_MODE == 1) {
      upper_limit <- ifelse(value_type == "Abs. Diff.", 42000, 120)
    } else if(SIMULATION_MODE == 2) {
      upper_limit <- ifelse(value_type == "Abs. Diff.", 30000, 50)
    }
    p <- ggplot(df, aes(x = Col_Label, y = Row_Label, fill = Value)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", Value)), color = "white", size = 5) +
      scale_fill_viridis_c(option = "viridis", name = value_type, limits = c(0, upper_limit)) +
      scale_x_discrete(position = "top") +
      scale_y_discrete(labels = function(x) as.character(x), limits = rev(levels(df$Row_Label))) +
      labs(x = "", y = "") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.grid = element_blank()
      )
    return(p)
  }
  
  # Prepare data for all four heatmaps
  # For AVERTED outcomes
  averted_inf_abs_df <- prepare_heatmap_data(abs_diff, 3, "infections", SIMULATION_MODE)  # row 3 = averted_infections
  averted_death_abs_df <- prepare_heatmap_data(abs_diff, 1, "deaths", SIMULATION_MODE)    # row 1 = averted_deaths
  averted_inf_pct_df <- prepare_heatmap_data(pct_diff, 3, "infections", SIMULATION_MODE)
  averted_death_pct_df <- prepare_heatmap_data(pct_diff, 1, "deaths", SIMULATION_MODE)
  
  # Create averted outcome heatmaps
  h1 <- create_heatmap(averted_inf_abs_df, "Abs. Diff.")
  h2 <- create_heatmap(averted_death_abs_df, "Abs. Diff.")
  h3 <- create_heatmap(averted_inf_pct_df, "% Diff.")
  h4 <- create_heatmap(averted_death_pct_df, "% Diff.")
  
  # Arrange plots 
  col1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Infection", size = 5, hjust=1) + theme_void() 
  col2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Death", size = 5, hjust=1) + theme_void() 
  
  layoutplot <- "
aaaaaaaabbbbbbbb
ggggggggeeeeeeee
ggggggggeeeeeeee
ggggggggeeeeeeee
ggggggggeeeeeeee
ggggggggeeeeeeee
ffffffffhhhhhhhh
ccccccccdddddddd
ccccccccdddddddd
ccccccccdddddddd
ccccccccdddddddd
ccccccccdddddddd
"
  
  plotlist_averted <-
    list(
      a = col1,
      b = col2,
      g = h1,
      e = h2,
      f = col1,
      h = col2,
      c = h3,
      d = h4)
  
  averted <- wrap_plots(plotlist_averted, design = layoutplot) 
  
  if(SIMULATION_MODE==1 | SIMULATION_MODE == 2) ggsave(paste0("2_figures/heatmap_averted_outcomes_", mode_suffix, "_strategy", strategy, ".png"), averted, width = 12, height = 7, dpi=300, units="in")

# PLOTTING SECTION =============================================================

# Helper function to create y-axis scales
create_y_scale <- function(upper_limit) {
  scale_y_continuous(limits = c(0, upper_limit))
}

# Set fixed y-axis limits based on simulation mode
if(SIMULATION_MODE == 1) {
  # Mode 1: Base simulations (9 scenarios) - uses fixed limits for all outcomes
  if(strategy==1) AVERTED_DEATH_LIMITS <- c(1000, 10000, 1.5e5, 1000, 10000, 1.5e5, 1000, 10000, 1.5e5)
  else if(strategy==2) AVERTED_DEATH_LIMITS <- c(2000, 20000, 2e5, 2000, 20000, 2e5, 2000, 20000, 2e5)
  if(strategy==1) AVERTIBLE_DEATH_LIMITS <- c(2000, 20000, 2e5, 2000, 20000, 2e5, 2000, 20000, 2e5)
  else if(strategy==2) AVERTIBLE_DEATH_LIMITS <- c(3000, 30000, 3e5, 3000, 30000, 3e5, 3000, 30000, 3e5)
  AVERTED_INFECTION_LIMITS <- c(1e5, 1e5, 1e5, 1, 1, 1, 1e5, 1e5, 1e5)
  if(strategy==1) AVERTIBLE_INFECTION_LIMITS <- c(2e5, 2e5, 2e5, 1, 1, 1, 2e5, 2e5, 2e5)
  else if(strategy==2) AVERTIBLE_INFECTION_LIMITS <- c(3e5, 3e5, 3e5, 1, 1, 1, 3e5, 3e5, 3e5)
  
  custom_y_averted_death <- map(AVERTED_DEATH_LIMITS, create_y_scale)
  custom_y_avertible_death <- map(AVERTIBLE_DEATH_LIMITS, create_y_scale)
  custom_y_averted_infections <- map(AVERTED_INFECTION_LIMITS, create_y_scale)
  custom_y_avertible_infections <- map(AVERTIBLE_INFECTION_LIMITS, create_y_scale)
  
} else if(SIMULATION_MODE == 2) {
  # Mode 2: Varying beta simulations (9 scenarios) - uses fixed limits only for averted infections
  AVERTED_INFECTION_LIMITS <- c(1e5, 1e5, 1e5, 1, 1, 1, 1e5, 1e5, 1e5)
  custom_y_averted_infections <- map(AVERTED_INFECTION_LIMITS, create_y_scale)
  
  # Death plots and avertible infections use free scales (no custom limits)
  AVERTED_DEATH_LIMITS <- rep(1e4, 9)
  custom_y_averted_death <- map(AVERTED_DEATH_LIMITS, create_y_scale)
  custom_y_avertible_death <- NULL
  custom_y_avertible_infections <- NULL
  
} else if(SIMULATION_MODE == 3) {
  # Mode 3: Realistic pathogen simulations (3 scenarios) - all plots use free scales
  custom_y_averted_death <- NULL
  custom_y_avertible_death <- NULL
  custom_y_averted_infections <- NULL
  custom_y_avertible_infections <- NULL
}

# Generate plots
averted_deaths <- ggplot() +
  geom_bar(data=res[,,1] %>% as.data.frame.table(), aes(y=Freq, x=type, group=scenarios), stat = "identity") +
  theme_minimal() +
  labs(y="Averted deaths", x="") +
  facet_wrap(.~scenarios, labeller=mylabel, scales = "free_y") +
  {if(!is.null(custom_y_averted_death)) facetted_pos_scales(y = custom_y_averted_death) else NULL} +
  theme(text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

avertible_deaths <- ggplot() +
  geom_bar(data=res[,,2] %>% as.data.frame.table(), aes(y=Freq, x=type, group=scenarios), stat = "identity") +
  theme_minimal() +
  labs(y="Avertible deaths", x="") +
  facet_wrap(.~scenarios, labeller=mylabel, scales = "free_y") +
  {if(!is.null(custom_y_avertible_death)) facetted_pos_scales(y = custom_y_avertible_death) else NULL} +
  theme(text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

averted_infections <- ggplot() +
  geom_bar(data=res[,,3] %>% as.data.frame.table(), aes(y=Freq, x=type, group=scenarios), stat = "identity") +
  theme_minimal() +
  labs(y="Averted infections", x="") +
  facet_wrap(.~scenarios, labeller=mylabel, scales = "free_y") +
  {if(!is.null(custom_y_averted_infections)) facetted_pos_scales(y = custom_y_averted_infections) else NULL} +
  theme(text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) 

avertible_infections <- ggplot() +
  geom_bar(data=res[,,4] %>% as.data.frame.table(), aes(y=Freq, x=type, group=scenarios), stat = "identity") +
  theme_minimal() +
  labs(y="Avertible infections", x="") +
  facet_wrap(.~scenarios, labeller=mylabel, scales = "free_y") +
  {if(!is.null(custom_y_avertible_infections)) facetted_pos_scales(y = custom_y_avertible_infections) else NULL} +
  theme(text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

# <VI> ARRANGE PLOTS -----------------------------------------------------------
# 6.1 Averted outcomes 
col1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Infection", size = 5, hjust=1) + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Death", size = 5, hjust=1) + theme_void() 

layoutplot <- "
aaaaaaaabbbbbbbb
ggggggggeeeeeeee
ggggggggeeeeeeee
ggggggggeeeeeeee
ggggggggeeeeeeee
ggggggggeeeeeeee
ggggggggeeeeeeee
ggggggggeeeeeeee
ggggggggeeeeeeee
ggggggggeeeeeeee
"

plotlist_averted <-
  list(
    a = col1,
    b = col2,
    g = averted_infections,
    e = averted_deaths)

averted <- wrap_plots(plotlist_averted, guides = 'collect', design = layoutplot) 

if(SIMULATION_MODE==1 | SIMULATION_MODE == 2) ggsave(paste0("2_figures/averted_outcomes_", mode_suffix, "_strategy", strategy, ".png"), averted, width = 22, height = 12, dpi=300, units="in")

# 6.2 Avertible outcomes 
plotlist_avertible <-
  list(
    a = col1,
    b = col2,
    g = avertible_infections,
    e = avertible_deaths)

avertible <- wrap_plots(plotlist_avertible, guides = 'collect', design = layoutplot) 

if(SIMULATION_MODE==1) ggsave(paste0("2_figures/avertible_outcomes_", mode_suffix, "_strategy", strategy, ".png"), avertible, width = 22, height = 12, dpi=300, units="in")

combined <- averted / avertible
if(SIMULATION_MODE==3) ggsave(paste0("2_figures/averted_avertible_outcomes_", mode_suffix, "_strategy", strategy, ".png"), combined, width = 22, height = 10, dpi=300, units="in")

# <VII> SUPPL PLOTS ------------------------------------------------------------
# 7.1 Dynamics for vax at 2 points
if(SIMULATION_MODE==1 & strategy==1){
  dynamics_all_scenarios <- dynamics_scenario_all %>%
    # filter(scenarios==1,
    #        comp_id=="D") %>%
    mutate(comp_label = case_when(comp_id == "S"~"1. Susceptible",
                                  comp_id == "I"~"2. Infectious",
                                  comp_id == "R"~"3. Recovered",
                                  comp_id == "D"~"4. Death",
                                  comp_id == "Cum_inf"~"5. Cumulative \n infections"),
           scenario_label = case_when(scenarios == "1"~"1: IFR=1%; \n VE inf=90%; \n  VE death=0",
                                      scenarios == "2"~"2: IFR=10%; \n VE inf=90%; \n  VE death=0",
                                      scenarios == "3"~"3: IFR=100%; \n VE inf=90%; \n  VE death=0",
                                      scenarios == "4"~"4: IFR=1%; \n VE inf=0; \n  VE death=90%",
                                      scenarios == "5"~"5: IFR=10%; \n VE inf=0; \n  VE death=90%",
                                      scenarios == "6"~"6: IFR=100%; \n VE inf=0; \n  VE death=90%",
                                      scenarios == "7"~"7: IFR=1%; \n VE inf=90%; \n  VE death=90%",
                                      scenarios == "8"~"5: IFR=10%; \n VE inf=90%; \n  VE death=90%",
                                      scenarios == "9"~"9: IFR=100%; \n VE inf=90%; \n  VE death=90%")) %>%
    ggplot() +
    geom_line(aes(x=as.numeric(time_id), y=Freq, col=vax_id, group=vax_id)) +
    scale_color_manual(values = c("1" = "forestgreen", 
                                  "61" = "darkblue", 
                                  "121" = "darkred")) +
    facet_grid(scenario_label~comp_label) +
    lims(x=c(0,121)) +
    theme_minimal() +
    labs(y="Probability",
         x="Day",
         col="Vax time") +
    theme(text = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)) +
    geom_vline(xintercept = 60, lty = "dashed", color = "gray50")
  
  ggsave(paste0("2_figures/dynamics_", mode_suffix, "_strategy", strategy, ".png"), dynamics_all_scenarios, width = 20, height = 13, units = "in", dpi = 300)
}


