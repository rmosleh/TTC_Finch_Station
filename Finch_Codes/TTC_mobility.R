# Load the necessary libraries
library(adaptivetau)
library(ggplot2)

# Define the parameters and initial states for each simulation period

# ODEs results
Inf_com_ode_rush_mor = c(2.798147700336436,2.705861847971072,2.630402551005460,
                        2.568701809782487) 

Inf_com_ode_off_mor = c(2.568701809782487,2.653955874306719,2.684282982757323,
                        2.695071552920243,2.698909289099528,2.700274417341106,
                        2.700760105091462) 

Inf_com_ode_rush_eve = c(2.700760105091462,2.625448563413102,2.599566029109999,
                        2.590671155023782, 2.587614151832164) 

Inf_com_ode_off_eve = c( 2.587614151832164, 2.651578971284700,2.692112753125213,
                        2.717798242339989, 2.734075024867917, 2.744389305633340,
                        2.750925390991835,2.755067202203305)
                        

Inf_hub_ode_rush_mor = c(0, 0.092285852365364, 0.167745149330976, 0.229445890553949) 

Inf_hub_ode_off_mor = c(0.229445890553949, 0.144191826029717,0.113864717579113,
                       0.103076147416193, 0.099238411236909,0.097873282995331,
                       0.097387595244974) 

Inf_hub_ode_rush_eve = c(0.097387595244974, 0.172699136923334,0.198581671226437,
                        0.207476545312654,0.210533548504272) 

Inf_hub_ode_off_eve = c( 0.210533548504272, 0.146568729051735, 0.106034947211223,
                        0.080349457996446, 0.064072675468518, 0.053758394703095,
                        0.047222309344599, 0.043080498133130)
# Exposed

Ex_com_ode_rush_mor = c(0,0,0,0) 

Ex_com_ode_off_mor = c(0,0,0,0,0,0,0) 

Ex_com_ode_rush_eve = c(0,0,0, 0, 0) 

Ex_com_ode_off_eve = c( 0, 0,0,0, 0, 0,0,0)


Ex_hub_ode_rush_mor = c(0, 0, 0, 0) 

Ex_hub_ode_off_mor = c(0, 0,0,0, 0,0,0) 

Ex_hub_ode_rush_eve = c(0, 0,0,0,0) 

Ex_hub_ode_off_eve = c( 0, 0, 0,0, 0, 0, 0, 0)

# Morning Rush Hours Simulation
parameters_rush <- c(alpha = 0.0364,  # Transmission rate
                     gamma = 0.1648 
                    ) # Recovery rate

initial_state_rush <- c(S_c = 131756 - 3,
                        L_c = 0, 
                        I_c = 3, 
                        R_c = 0,
                        S_h = 610,
                        L_h = 0,
                        I_h = 0,
                        R_h = 0)

# Define the transition matrix
transitions <- list(
  c(S_c = -1, S_h = +1),  
  c(L_c = -1, L_h = +1),
  c(I_c = -1, I_h = +1),  
  c(R_c = -1, R_h = +1),
  c(S_c = +1, S_h = -1),
  c(L_c = +1, L_h = -1),
  c(I_c = +1, I_h = -1),
  c(R_c = +1, R_h = -1)
)

# Define the rate function
rate_function <- function(state, params, t) {
  with(as.list(c(state, params)), {
    return(c(alpha * S_c,   # Rate of infection
             alpha * L_c,
             alpha * I_c,
             alpha * R_c,
             gamma * S_h,
             gamma * L_h,
             gamma * I_h,
             gamma * R_h))  # Rate of recovery
  })
}

# Run the Morning Rush Hours simulation
set.seed(1)
simulation_result_rush <- ssa.adaptivetau(
  init.values = initial_state_rush,
  transitions = transitions,
  rateFunc = rate_function,
  params = parameters_rush,
  tf = 3
)

# Convert results to data frame
simulation_df_rush <- as.data.frame(simulation_result_rush)

# Extract final states to use as initial state for Morning Off-Peak Hours
last_S_c_rush <- tail(simulation_df_rush$S_c, n = 1)
last_L_c_rush <- tail(simulation_df_rush$L_c, n = 1)
last_I_c_rush <- tail(simulation_df_rush$I_c, n = 1)
last_R_c_rush <- tail(simulation_df_rush$R_c, n = 1)
last_S_h_rush <- tail(simulation_df_rush$S_h, n = 1)
last_L_h_rush <- tail(simulation_df_rush$L_h, n = 1)
last_I_h_rush <- tail(simulation_df_rush$I_h, n = 1)
last_R_h_rush <- tail(simulation_df_rush$R_h, n = 1)

# Morning Off-Peak Hours Simulation
parameters_offpeak <- c(alpha = 0.0358,  # Transmission rate
                        gamma = 0.9977 
                       )

initial_state_offpeak <- c(S_c = last_S_c_rush,
                           L_c = last_L_c_rush,
                           I_c = last_I_c_rush,
                           R_c = last_R_c_rush,
                           S_h = last_S_h_rush,
                           L_h = last_L_h_rush,
                           I_h = last_I_h_rush,
                           R_h = last_R_h_rush)

# Run the Morning Off-Peak Hours simulation
set.seed(1)
simulation_result_offpeak <- ssa.adaptivetau(
  init.values = initial_state_offpeak,
  transitions = transitions,
  rateFunc = rate_function,
  params = parameters_offpeak,
  tf = 6
)

# Convert results to data frame
simulation_df_offpeak <- as.data.frame(simulation_result_offpeak)

# Adjust time in Off-Peak Hours result to start from 3
simulation_df_offpeak$time <- simulation_df_offpeak$time + 3

# Extract final states to use as initial state for Evening Rush Hours
last_S_c_offpeak <- tail(simulation_df_offpeak$S_c, n = 1)
last_L_c_offpeak <- tail(simulation_df_offpeak$L_c, n = 1)
last_I_c_offpeak <- tail(simulation_df_offpeak$I_c, n = 1)
last_R_c_offpeak <- tail(simulation_df_offpeak$R_c, n = 1)
last_S_h_offpeak <- tail(simulation_df_offpeak$S_h, n = 1)
last_L_h_offpeak <- tail(simulation_df_offpeak$L_h, n = 1)
last_I_h_offpeak <- tail(simulation_df_offpeak$I_h, n = 1)
last_R_h_offpeak <- tail(simulation_df_offpeak$R_h, n = 1)

# Evening Rush Hours Simulation
parameters_rush_eve <- c(alpha = 0.0809,  # Transmission rate
                         gamma = 0.9870
                         ) # Recovery rate

initial_state_rush_eve <- c(S_c = last_S_c_offpeak,
                            L_c = last_L_c_offpeak, 
                            I_c = last_I_c_offpeak, 
                            R_c = last_R_c_offpeak,
                            S_h = last_S_h_offpeak,
                            L_h = last_L_h_offpeak,
                            I_h = last_I_h_offpeak,
                            R_h = last_R_h_offpeak)

# Run the Evening Rush Hours simulation
set.seed(1)
simulation_result_rush_eve <- ssa.adaptivetau(
  init.values = initial_state_rush_eve,
  transitions = transitions,
  rateFunc = rate_function,
  params = parameters_rush_eve,
  tf = 4
)

# Convert results to data frame
simulation_df_rush_eve <- as.data.frame(simulation_result_rush_eve)

# Adjust time in Evening Rush Hours result to start from 9
simulation_df_rush_eve$time <- simulation_df_rush_eve$time + 9

# Extract final states to use as initial state for Evening Off-Peak Hours
last_S_c_rush_eve <- tail(simulation_df_rush_eve$S_c, n = 1)
last_L_c_rush_eve <- tail(simulation_df_rush_eve$L_c, n = 1)
last_I_c_rush_eve <- tail(simulation_df_rush_eve$I_c, n = 1)
last_R_c_rush_eve <- tail(simulation_df_rush_eve$R_c, n = 1)
last_S_h_rush_eve <- tail(simulation_df_rush_eve$S_h, n = 1)
last_L_h_rush_eve <- tail(simulation_df_rush_eve$L_h, n = 1)
last_I_h_rush_eve <- tail(simulation_df_rush_eve$I_h, n = 1)
last_R_h_rush_eve <- tail(simulation_df_rush_eve$R_h, n = 1)

# Evening Off-Peak Hours Simulation
parameters_offpeak_eve <- c(alpha = 0.0058,  # Transmission rate
                            gamma = 0.4503 
                            ) # Recovery rate

initial_state_offpeak_eve <- c(S_c = last_S_c_rush_eve,
                               L_c = last_L_c_rush_eve, 
                               I_c = last_I_c_rush_eve, 
                               R_c = last_R_c_rush_eve,
                               S_h = last_S_h_rush_eve,
                               L_h = last_L_h_rush_eve,
                               I_h = last_I_h_rush_eve,
                               R_h = last_R_h_rush_eve)

# Run the Evening Off-Peak Hours simulation
set.seed(1)
simulation_result_offpeak_eve <- ssa.adaptivetau(
  init.values = initial_state_offpeak_eve,
  transitions = transitions,
  rateFunc = rate_function,
  params = parameters_offpeak_eve,
  tf = 7
)

# Convert results to data frame
simulation_df_offpeak_eve <- as.data.frame(simulation_result_offpeak_eve)

# Adjust time in Evening Off-Peak Hours result to start from 13
simulation_df_offpeak_eve$time <- simulation_df_offpeak_eve$time + 13

# Combine all simulation results
combined_df <- rbind(simulation_df_rush, simulation_df_offpeak, simulation_df_rush_eve, simulation_df_offpeak_eve)

# Create a separate data frame for plotting ODEs results
combined_ode_df <- data.frame(
  time = c(seq(0, 3, length.out = length(Inf_com_ode_rush_mor)), 
           seq(3, 9, length.out = length(Inf_com_ode_off_mor)),
           seq(9, 13, length.out = length(Inf_com_ode_rush_eve)),
           seq(13, 20, length.out = length(Inf_com_ode_off_eve))),
  I_c = c(Inf_com_ode_rush_mor, Inf_com_ode_off_mor, Inf_com_ode_rush_eve, Inf_com_ode_off_eve),
  Type = rep("ODE_com", length(c(Inf_com_ode_rush_mor, Inf_com_ode_off_mor, Inf_com_ode_rush_eve, Inf_com_ode_off_eve)))
)

combined_hub_df <- data.frame(
  time = c(seq(0, 3, length.out = length(Inf_hub_ode_rush_mor)), 
           seq(3, 9, length.out = length(Inf_hub_ode_off_mor)),
           seq(9, 13, length.out = length(Inf_hub_ode_rush_eve)),
           seq(13, 20, length.out = length(Inf_hub_ode_off_eve))),
  I_h = c(Inf_hub_ode_rush_mor, Inf_hub_ode_off_mor, Inf_hub_ode_rush_eve, Inf_hub_ode_off_eve),
  Type = rep("ODE_hub", length(c(Inf_hub_ode_rush_mor, Inf_hub_ode_off_mor, Inf_hub_ode_rush_eve, Inf_hub_ode_off_eve)))
)

# Plot the combined results
ggplot() +
 geom_line(data = combined_df, aes(x = time, y = I_c, color = "CTMC_Infctee_Com")) +
  geom_line(data = combined_df, aes(x = time, y = I_h, color = "CTMC_Infectee_Hub")) +
  geom_line(data = combined_ode_df, aes(x = time, y = I_c, color = "ODE_Infectee_Com")) +
  geom_line(data = combined_hub_df, aes(x = time, y = I_h, color = "ODE_Infectee_Hub")) +
  labs(
    x = "Time (Hour)",
    y = "Number of the Exposed Individuals") +
  scale_color_manual(values = c("CTMC_Infctee_Com" = "red", "CTMC_Infectee_Hub" = "blue", "ODE_Infectee_Com" = "black", "ODE_Infectee_Hub" = "green")) +
  theme_minimal()

#combined_ode_df <- data.frame(
#  time = c(seq(0, 3, length.out = length(Ex_com_ode_rush_mor)), 
       #    seq(3, 9, length.out = length(Ex_com_ode_off_mor)),
     #      seq(9, 13, length.out = length(Ex_com_ode_rush_eve)),
      #     seq(13, 20, length.out = length(Ex_com_ode_off_eve))),
#  L_c = c(Ex_com_ode_rush_mor, Ex_com_ode_off_mor, Ex_com_ode_rush_eve, Ex_com_ode_off_eve),
#  Type = rep("ODE_com", length(c(Ex_com_ode_rush_mor, Ex_com_ode_off_mor, Ex_com_ode_rush_eve, Ex_com_ode_off_eve)))
#)

#combined_hub_df <- data.frame(
 # time = c(seq(0, 3, length.out = length(Ex_hub_ode_rush_mor)), 
 #          seq(3, 9, length.out = length(Ex_hub_ode_off_mor)),
 #          seq(9, 13, length.out = length(Ex_hub_ode_rush_eve)),
 #          seq(13, 20, length.out = length(Ex_hub_ode_off_eve))),
 # L_h = c(Ex_hub_ode_rush_mor, Ex_hub_ode_off_mor, Ex_hub_ode_rush_eve, Ex_hub_ode_off_eve),
 # Type = rep("ODE_hub", length(c(Ex_hub_ode_rush_mor, Ex_hub_ode_off_mor, Ex_hub_ode_rush_eve, Ex_hub_ode_off_eve)))
#)

# Plot the combined results
#ggplot() +
 # geom_line(data = combined_df, aes(x = time, y = L_c, color = "CTMC_Infctee_Com")) +
#  geom_line(data = combined_df, aes(x = time, y = L_h, color = "CTMC_Infectee_Hub")) +
#  geom_line(data = combined_ode_df, aes(x = time, y = L_c, color = "ODE_Infectee_Com")) +
#  geom_line(data = combined_hub_df, aes(x = time, y = L_h, color = "ODE_Infectee_Hub")) +
#  labs(
 #   x = "Time (Hour)",
 #   y = "Number of the Exposed Individuals") +
 # scale_color_manual(values = c("CTMC_Infctee_Com" = "red", "CTMC_Infectee_Hub" = "blue", "ODE_Infectee_Com" = "black", "ODE_Infectee_Hub" = "green")) +
 # theme_minimal()