# Load the necessary libraries
library(adaptivetau)
library(ggplot2)
library(dplyr)


# ODEs results
Ex_com_ode_rush_mor = c( 0,
                        4.175579609005524,
                        14.924025213988539,
                        31.620628849789711) 

Ex_com_ode_off_mor = c(0.031620628849790,
                       0.806392133618238,
                       1.380480337287443,
                       1.835501978145210,
                       2.219197846423139,
                       2.559494225246198,
                       2.872730274725067) * 1.0e+03

Ex_com_ode_rush_eve = c(2.872730274725067,
                        3.749990675880330,
                        5.037799489443120,
                        6.419540570682869,
                        7.800374764646471) * 1.0e+03

Ex_com_ode_off_eve = c(7.800374764646471,
                       8.417001310356166,
                       8.823501503801852,
                       9.087264717588353,
                       9.254112606631676,
                       9.355192783780376,
                       9.411683389160315,
                       9.437980268879517) * 1.0e+03

Ex_hub_ode_rush_mor = c(   0,
                           3.256474102614896,
                           6.359772029467082,
                           9.319231651509556) * 1.0e+02

Ex_hub_ode_off_mor = c(9.319231651509556,
                       7.011114272824757,
                       5.735314121974124,
                       5.051092865626721,
                       4.705109174025061,
                       4.551919847132165,
                       4.508140506770838) * 1.0e+02

Ex_hub_ode_rush_eve = c(0.450814050677084,
                        0.978802789632168,
                        1.232430204528012,
                        1.406971986757701,
                        1.557324329020907) * 1.0e+03

Ex_hub_ode_off_eve = c( 1.557324329020907,
                        1.127633021479866,
                        0.835943831635315,
                        0.637883864608475,
                        0.503343503250731,
                        0.411897647302102,
                        0.349687200026037,
                        0.307310689990465) * 1.0e+03



# Parameters for simulation
# Morning Rush Hours Simulation
parameters_rush <- c(alpha = 0.0364,  
                     gamma = 0.0114, 
                     beta_c = 0.0066,
                     k_c = 22598.2129,
                     epsilon_c = 0.0116,
                     delta_c = 0.0019,
                     beta_h = 0.0116, 
                     k_h = 37068.2276,
                     epsilon_h = 0.00116, 
                     delta_h = 0.00019,
                     p = 0.0916) 

# Morning Off-Peak Hours Simulation
parameters_offpeak <- c(alpha = 0.0196,  # Transmission rate
                        gamma = 0.4998, 
                        beta_c = 0.0066,
                        k_c = 121851.9505,
                        epsilon_c = 0.0116,
                        delta_c = 0.0019,
                        beta_h = 0.0023, 
                        k_h = 784.0781,
                        epsilon_h = 0.00116, 
                        delta_h = 0.00019,
                        p = 0.0916)

# Evening Rush Hours Simulation
parameters_rush_eve <- c(alpha = 0.1108,  # Transmission rate
                         gamma = 0.9990, 
                         beta_c = 0.0066,
                         k_c = 121851.9505,
                         epsilon_c = 0.0116,
                         delta_c = 0.0019,
                         beta_h = 0.1990, 
                         k_h = 37068.2276,
                         epsilon_h = 0.00116, 
                         delta_h = 0.00019,
                         p =  0.0916) # Recovery rate

parameters_offpeak_eve <- c(alpha = 0.0046,  # Transmission rate
                            gamma = 0.3776, 
                            beta_c =0.0066,
                            k_c = 380859.6078,
                            epsilon_c = 0.0116,
                            delta_c = 0.0019,
                            beta_h = 0.0017, 
                            k_h = 784.0781,
                            epsilon_h = 0.00116, 
                            delta_h = 0.00019,
                            p = 0.0916) 

initial_state_rush <- c(S_c = 100006 - 10,
                        L_c = 0, 
                        I_c = 10, 
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
  c(R_c = +1, R_h = -1),
  c(S_c = -1, L_c = +1),
  c(S_c = -1, L_h = +1),
  c(L_c = -1, I_c = +1),
  c(I_c = -1, R_c = +1),
  c(S_h = -1, L_h = +1),
  c(S_h = -1, L_c = +1),
  c(L_h = -1, I_h = +1),
  c(I_h = -1, R_h = +1)
)

# Define the rate function
rate_function <- function(state, params, t) {
  with(as.list(c(state, params)), {
    return(c(alpha * (1-p) * S_c,   
             alpha * L_c,
             alpha * I_c,
             alpha * R_c,
             gamma * (1-p) * S_h,
             gamma * L_h,
             gamma * I_h,
             gamma * R_h,
             beta_c * ((I_c + I_h) / (1 + k_c * (I_c + I_h))) * S_c,
             p * alpha * S_c,
             epsilon_c * L_c,
             delta_c * I_c,
             beta_h * ((I_c + I_h) / (1 + k_h * (I_c + I_h))) * S_h,
             p * gamma * S_h,
             epsilon_h * L_h,
             delta_h * I_h)) 
  })
}

# Run multiple simulations and average results
# Fixing the averaging of L_c and R_c
run_simulation <- function(initial_state, transitions, rate_function, params, tf, n_simulations = 50) {
  all_simulations <- replicate(n_simulations, {
    simulation_result <- ssa.adaptivetau(init.values = round(initial_state), transitions = transitions,
                                         rateFunc = rate_function, params = params, tf = tf)
    as.data.frame(simulation_result)
  }, simplify = FALSE)
  
  # Interpolating over a fixed time grid and averaging
  time_grid <- seq(0, tf, length.out = 100)
  
  averaged_simulation_I_c <- lapply(all_simulations, function(df) {
    approx(df$time, df$I_c, time_grid)$y
  })
  averaged_I_c <- rowMeans(do.call(cbind, averaged_simulation_I_c))
  
  averaged_simulation_I_h <- lapply(all_simulations, function(df) {
    approx(df$time, df$I_h, time_grid)$y
  })
  averaged_I_h <- rowMeans(do.call(cbind, averaged_simulation_I_h))
  
  averaged_simulation_S_c <- lapply(all_simulations, function(df) {
    approx(df$time, df$S_c, time_grid)$y
  })
  averaged_S_c <- rowMeans(do.call(cbind, averaged_simulation_S_c))
  
  averaged_simulation_S_h <- lapply(all_simulations, function(df) {
    approx(df$time, df$S_h, time_grid)$y
  })
  averaged_S_h <- rowMeans(do.call(cbind, averaged_simulation_S_h))
  
  averaged_simulation_L_c <- lapply(all_simulations, function(df) {
    approx(df$time, df$L_c, time_grid)$y
  })
  averaged_L_c <- rowMeans(do.call(cbind, averaged_simulation_L_c))  # Corrected from averaging L_h
  
  averaged_simulation_L_h <- lapply(all_simulations, function(df) {
    approx(df$time, df$L_h, time_grid)$y
  })
  averaged_L_h <- rowMeans(do.call(cbind, averaged_simulation_L_h))
  
  averaged_simulation_R_h <- lapply(all_simulations, function(df) {
    approx(df$time, df$R_h, time_grid)$y
  })
  averaged_R_h <- rowMeans(do.call(cbind, averaged_simulation_R_h))
  
  averaged_simulation_R_c <- lapply(all_simulations, function(df) {
    approx(df$time, df$R_c, time_grid)$y
  })
  averaged_R_c <- rowMeans(do.call(cbind, averaged_simulation_R_c))  # Corrected from averaging R_h
  
  return(data.frame(time = time_grid, I_c = averaged_I_c, I_h = averaged_I_h, 
                    S_c = averaged_S_c, S_h = averaged_S_h, 
                    L_h = averaged_L_h, L_c = averaged_L_c, R_h = averaged_R_h, R_c = averaged_R_c))
}


  
  # Running the simulations for each time period and combining the results
  simulation_rush <- run_simulation(initial_state_rush, transitions, rate_function, parameters_rush, tf = 3)
  simulation_offpeak <- run_simulation(c(S_c = tail(simulation_rush$S_c, 1), L_c = tail(simulation_rush$L_c, 1), 
                                         I_c = tail(simulation_rush$I_c, 1), 
                                         R_c =tail(simulation_rush$R_c, 1), S_h = tail(simulation_rush$S_h, 1), 
                                         L_h = tail(simulation_rush$L_h, 1), 
                                         I_h = tail(simulation_rush$I_h, 1), R_h = tail(simulation_rush$R_h, 1)),
                                       transitions, rate_function, parameters_offpeak, tf = 6)
  simulation_rush_eve <- run_simulation(c(S_c = tail(simulation_offpeak$S_c, 1), L_c =  tail(simulation_offpeak$L_c, 1),
                                          I_c = tail(simulation_offpeak$I_c, 1),
                                          R_c =  tail(simulation_offpeak$R_c, 1), S_h = tail(simulation_offpeak$S_h, 1),
                                          L_h =  tail(simulation_offpeak$L_h, 1), 
                                          I_h = tail(simulation_offpeak$I_h, 1), R_h =  tail(simulation_offpeak$R_h, 1)),
                                        transitions, rate_function, parameters_rush_eve, tf = 4)
  simulation_offpeak_eve <- run_simulation(c(S_c = tail(simulation_rush_eve$S_c, 1), L_c =  tail(simulation_rush_eve$L_c, 1), 
                                             I_c = tail(simulation_rush_eve$I_c, 1), 
                                             R_c =  tail(simulation_rush_eve$R_c, 1), S_h = tail(simulation_rush_eve$S_h, 1), 
                                             L_h =  tail(simulation_rush_eve$L_h, 1), 
                                             I_h = tail(simulation_rush_eve$I_h, 1), R_h =  tail(simulation_rush_eve$R_h, 1)),
                                           transitions, rate_function, parameters_offpeak_eve, tf = 7)
  


  # Fix the combined simulation dataframe and plotting
  combined_simulation <- bind_rows(
    mutate(simulation_rush, time = time + 0),
    mutate(simulation_offpeak, time = time + 3),
    mutate(simulation_rush_eve, time = time + 9),
    mutate(simulation_offpeak_eve, time = time + 13)
  )
  
# ODE data (as before)
combined_ode_df <- data.frame(
  time = c(seq(0, 3, length.out = length(Ex_com_ode_rush_mor)), seq(3, 9, length.out = length(Ex_com_ode_off_mor)),
           seq(9, 13, length.out = length(Ex_com_ode_rush_eve)), seq(13, 20, length.out = length(Ex_com_ode_off_eve))),
  L_c = c(Ex_com_ode_rush_mor, Ex_com_ode_off_mor, Ex_com_ode_rush_eve, Ex_com_ode_off_eve),
  Type = rep("ODE_com", length(c(Ex_com_ode_rush_mor, Ex_com_ode_off_mor, Ex_com_ode_rush_eve, Ex_com_ode_off_eve)))
)

combined_hub_df <- data.frame(
  time = c(seq(0, 3, length.out = length(Ex_hub_ode_rush_mor)), seq(3, 9, length.out = length(Ex_hub_ode_off_mor)),
           seq(9, 13, length.out = length(Ex_hub_ode_rush_eve)), seq(13, 20, length.out = length(Ex_hub_ode_off_eve ))),
 L_h = c(Ex_hub_ode_rush_mor, Ex_hub_ode_off_mor, Ex_hub_ode_rush_eve, Ex_hub_ode_off_eve ),
  Type = rep("ODE_hub", length(c(Ex_hub_ode_rush_mor, Ex_hub_ode_off_mor, Ex_hub_ode_rush_eve, Ex_hub_ode_off_eve)))
)
# Plot the combined results with corrected dataframe reference
ggplot() +
  geom_line(data = combined_simulation, aes(x = time, y = L_c, color = "CTMC_Exposed_Com" )) +
  geom_line(data = combined_simulation, aes(x = time, y = L_h, color = "CTMC_Exposed_Hub")) +
  geom_line(data = combined_ode_df, aes(x = time, y = L_c, color = "ODE_Exposed_Com")) +
  geom_line(data = combined_hub_df, aes(x = time, y = L_h, color = "ODE_Exposed_Hub")) +
  labs(
    x = "Time (Hour)",
    y = "Number of the Exposed Individuals") +
  scale_color_manual(values = c("CTMC_Exposed_Com" = "red", "CTMC_Exposed_Hub" = "blue", "ODE_Exposed_Com" = "black", "ODE_Exposed_Hub" = "green")) +
  theme_minimal()