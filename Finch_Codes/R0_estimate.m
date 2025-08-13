function [R0]= R0_estimate(par)

%% this program evaluate R_0 over a single day based on the estimated parameters

beta_c=par(1);        % transmission rate of the community
beta_h_rush_mor=par(2);   % transmission rate of the hub in morning rush hours
beta_h_off_mor=par(3);    % transmission rate of the hub in morning off-peak hours
beta_h_rush_eve=par(4);   % transmission rate of the hub in evening rush hours
beta_h_off_eve=par(5);    % transmission rate of the hub in evening off-peak hours
alpha_rush_mor=par(6);    % inflow rate  from community to the hub in morning rush hours
alpha_off_mor=par(7);     % inflow rate  from community to the hub in morning off-peak hours
alpha_rush_eve=par(8);    % inflow rate  from community to the hub in evening rush hours
alpha_off_eve=par(9);     % inflow rate  from community to the hub in evening off-peak hours
gamma_rush_mor=par(10);    %  outflow rate  from hub to the community in morning rush hours
gamma_off_mor=par(11);     % outflow rate  from hub to the community in morning off-peak hours
gamma_rush_eve=par(12);    %  outflow rate  from hub to the community in evening rush hours
gamma_off_eve=par(13);  
k_off=par(14);         % satuartion number for the force of infection of the community  in off-peak hours
k_rush=par(15);        % satuartion number for the force of infection of the community  in rush hours 
G_off=par(16);        % satuartion number for the force of infection of the hub  in off-peak hours
G_rush=par(17);       % satuartion number for the force of infection of the hub  in rush hours
 N=par(18);           % Number of the total population of the community
 I_0c=par(19);      % Initial data of infectees in the community

epsilon_c= 0.0116;    % propagation rate for community
epsilon_h= 0.00116;   % propagation rate for hub
delta_c=0.0019;       % recovery rate for community
delta_h=0.00019;      % recovery rate for hub

% R_O for the morning rush hours

P_rush_mor=alpha_rush_mor^2*beta_c*epsilon_h+alpha_rush_mor*beta_c*delta_c*epsilon_h+alpha_rush_mor*beta_c*epsilon_c*epsilon_h+...
beta_c*delta_h*epsilon_c*epsilon_h+alpha_rush_mor*beta_c*epsilon_c*gamma_rush_mor+beta_c*delta_c*epsilon_c*gamma_rush_mor+...
alpha_rush_mor*beta_c*epsilon_h*gamma_rush_mor+beta_c*epsilon_c*epsilon_h*gamma_rush_mor+beta_c*epsilon_c*gamma_rush_mor^2;

R_rush_mor=alpha_rush_mor^2*beta_h_rush_mor*epsilon_h+alpha_rush_mor*beta_h_rush_mor*delta_c*epsilon_h+alpha_rush_mor*beta_h_rush_mor*epsilon_c*epsilon_h+...
beta_c*delta_c*epsilon_c*epsilon_h+alpha_rush_mor*beta_h_rush_mor*epsilon_c*gamma_rush_mor+beta_h_rush_mor*delta_h*epsilon_c*gamma_rush_mor+...
alpha_rush_mor*beta_h_rush_mor*epsilon_h*gamma_rush_mor+beta_c*epsilon_c*epsilon_h*gamma_rush_mor+beta_h_rush_mor*epsilon_c*gamma_rush_mor^2;

R_0_rush_mor=(P_rush_mor+(alpha_rush_mor/gamma_rush_mor)*R_rush_mor)/...
((alpha_rush_mor*delta_h+delta_h*delta_c+delta_c*gamma_rush_mor)*(alpha_rush_mor*epsilon_h+epsilon_c*epsilon_h+epsilon_c*gamma_rush_mor))

% R_O for the morning off-peak hours

P_off_mor=(alpha_off_mor)^2*beta_c*epsilon_h+alpha_off_mor*beta_c*delta_c*epsilon_h+alpha_off_mor*beta_c*epsilon_c*epsilon_h+...
beta_c*delta_h*epsilon_c*epsilon_h+alpha_off_mor*beta_c*epsilon_c*gamma_off_mor+beta_c*delta_c*epsilon_c*gamma_off_mor+...
alpha_off_mor*beta_c*epsilon_h*gamma_off_mor+beta_c*epsilon_c*epsilon_h*gamma_off_mor+beta_c*epsilon_c*gamma_off_mor^2;

R_off_mor=alpha_off_mor^2*beta_h_off_mor*epsilon_h+alpha_off_mor*beta_h_off_mor*delta_c*epsilon_h+alpha_off_mor*beta_h_off_mor*epsilon_c*epsilon_h+...
beta_c*delta_c*epsilon_c*epsilon_h+alpha_off_mor*beta_h_off_mor*epsilon_c*gamma_off_mor+beta_h_off_mor*delta_h*epsilon_c*gamma_off_mor+...
alpha_off_mor*beta_h_off_mor*epsilon_h*gamma_off_mor+beta_c*epsilon_c*epsilon_h*gamma_off_mor+beta_h_off_mor*epsilon_c*gamma_off_mor^2;

R_0_off_mor=(P_off_mor+(alpha_off_mor/gamma_off_mor)*R_off_mor)/...
((alpha_off_mor*delta_h+delta_h*delta_c+delta_c*gamma_off_mor)*(alpha_off_mor*epsilon_h+epsilon_c*epsilon_h+epsilon_c*gamma_off_mor))


% R_O for the evening rush hours

P_rush_eve=alpha_rush_eve^2*beta_c*epsilon_h+alpha_rush_eve*beta_c*delta_c*epsilon_h+alpha_rush_eve*beta_c*epsilon_c*epsilon_h+...
beta_c*delta_h*epsilon_c*epsilon_h+alpha_rush_eve*beta_c*epsilon_c*gamma_rush_eve+beta_c*delta_c*epsilon_c*gamma_rush_eve+...
alpha_rush_eve*beta_c*epsilon_h*gamma_rush_eve+beta_c*epsilon_c*epsilon_h*gamma_rush_eve+beta_c*epsilon_c*gamma_rush_eve^2;

R_rush_eve=alpha_rush_eve^2*beta_h_rush_eve*epsilon_h+alpha_rush_eve*beta_h_rush_eve*delta_c*epsilon_h+alpha_rush_eve*beta_h_rush_eve*epsilon_c*epsilon_h+...
beta_c*delta_c*epsilon_c*epsilon_h+alpha_rush_eve*beta_h_rush_eve*epsilon_c*gamma_rush_eve+beta_h_rush_eve*delta_h*epsilon_c*gamma_rush_eve+...
alpha_rush_eve*beta_h_rush_eve*epsilon_h*gamma_rush_eve+beta_c*epsilon_c*epsilon_h*gamma_rush_eve+beta_h_rush_eve*epsilon_c*gamma_rush_eve^2;

R_0_rush_eve=(P_rush_eve+(alpha_rush_eve/gamma_rush_eve)*R_rush_eve)/...
((alpha_rush_eve*delta_h+delta_h*delta_c+delta_c*gamma_rush_eve)*(alpha_rush_eve*epsilon_h+epsilon_c*epsilon_h+epsilon_c*gamma_rush_eve))

% R_O for the evening off-peak hours

P_off_eve=alpha_off_eve^2*beta_c*epsilon_h+alpha_off_eve*beta_c*delta_c*epsilon_h+alpha_off_eve*beta_c*epsilon_c*epsilon_h+...
beta_c*delta_h*epsilon_c*epsilon_h+alpha_off_eve*beta_c*epsilon_c*gamma_off_eve+beta_c*delta_c*epsilon_c*gamma_off_eve+...
alpha_off_eve*beta_c*epsilon_h*gamma_off_eve+beta_c*epsilon_c*epsilon_h*gamma_off_eve+beta_c*epsilon_c*gamma_off_eve^2;

R_off_eve=alpha_off_eve^2*beta_h_off_eve*epsilon_h+alpha_off_eve*beta_h_off_eve*delta_c*epsilon_h+alpha_off_eve*beta_h_off_eve*epsilon_c*epsilon_h+...
beta_c*delta_c*epsilon_c*epsilon_h+alpha_off_eve*beta_h_off_eve*epsilon_c*gamma_off_eve+beta_h_off_eve*delta_h*epsilon_c*gamma_off_eve+...
alpha_off_eve*beta_h_off_eve*epsilon_h*gamma_off_eve+beta_c*epsilon_c*epsilon_h*gamma_off_eve+beta_h_off_eve*epsilon_c*gamma_off_eve^2;

R_0_off_eve=(P_off_eve+(alpha_off_eve/gamma_off_eve)*R_off_eve)/...
((alpha_off_eve*delta_h+delta_h*delta_c+delta_c*gamma_off_eve)*(alpha_off_eve*epsilon_h+epsilon_c*epsilon_h+epsilon_c*gamma_off_eve))
R0=max([R_0_rush_mor,R_0_off_mor,R_0_rush_eve,R_0_off_eve])
