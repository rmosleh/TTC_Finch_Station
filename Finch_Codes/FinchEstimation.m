function [pa_estimate fval history] = FinchEstimation
format long

%% Estimating Process for parameters in the community and hub(Finch station) over a single day

 history=[];
 Data=load('Data_mod.txt');

%-----------Estimate_Parameters_initial value--------

beta_c_in=0.003;           % initial value for transmission rate of the community
beta_h_rush_mor_in=0.82;    % initial value for transmission rate of the hub in morning rush hours
beta_h_off_mor_in=0.0001;   % initial value for transmission rate of the hub in morning off-peak hours
beta_h_rush_eve_in=0.8;     % initial value for transmission rate of the hub in evening rush hours
beta_h_off_eve_in=0.0001;   % initial value for transmission rate of the hub in evening off-peak hours
alpha_rush_mor_in=0.908;    % initial value for inflow rate  from community to the hub in morning rush hours
alpha_off_mor_in=0.733;     %  initial value for inflow rate  from community to the hub in morning off-peak                         
alpha_rush_eve_in=0.098;    % initial value for inflow rate  from community to the hub in evening rush hours
alpha_off_eve_in=0.128;     % initial value for inflow rate  from community to the hub in evening off-peak hours
gamma_rush_mor_in=0.0102;   % initial value for outflow rate  from hub to the community in morning rush hours
gamma_off_mor_in=0.0107;    % initial value for outflow rate  from hub to the community in morning off-peak hours
gamma_rush_eve_in=0.0755;   %  initial value for outflow rate  from hub to the community in evening rush hours
gamma_off_eve_in=0.0703;    % initial value for outflow rate  from hub to the community in evening off-peak hours

k_off_in=195;               % initial value for satuartion number for the force of infection of the community  in off-peak hours
k_rush_in=20;               % initial value for satuartion number for the force of infection of the community  in rush hours  =
G_off_in=8000;              % initial value for satuartion number for the force of infection of the hub  in off-peak hours
G_rush_in=100;              % initial value for satuartion number for the force of infection of the hub  in rush hours
N_in=200000;                % initial value for number of the total population of the community
I_c_in=2;                   % initial value for I_c(0)
p_in=0.01;                  % initial value for percentage of individuals carryig infectious

param_estimate=[beta_c_in,  beta_h_rush_mor_in, beta_h_off_mor_in,  beta_h_rush_eve_in, ...
    beta_h_off_eve_in,alpha_rush_mor_in, alpha_off_mor_in, alpha_rush_eve_in, alpha_off_eve_in,...
gamma_rush_mor_in, gamma_off_mor_in, gamma_rush_eve_in, gamma_off_eve_in,...
k_off_in, k_rush_in, G_off_in, G_rush_in, N_in,I_c_in,p_in ]; 
%-----------------Estimating Function ----------------------------------
LB=[zeros(1,13), ones(1,4), 100000, 0,0 ]; %lower Bound
UB=[ ones(1,13),1000000*ones(1,4), 1000000, 10,1]; % Upper Bound

options = optimset('OutputFcn', @myoutput, 'TolFun', 1e-6);
[pa_estimate,fval]=fminsearchbnd(@(param)fmin_estimatettc_one(param), param_estimate, ...  % Estimating function
    LB,UB,options)
[error_estimate  x_rush_mor x_off_mor x_rush_eve x_off_eve ]=fmin_estimatettc_one(pa_estimate);

%---------------Plotting Process------------------------

t_rush_mor=linspace(0,3,4);
t_off_mor=linspace(3,9,7);
t_rush_eve=linspace(9,13,5);
t_off_eve=linspace(13,20,8);
t_total=linspace(0,20,21);

 error_estimate_rush_mor=zeros(length(t_rush_mor),1);
 error_estimate_off_mor=zeros(length(t_off_mor),1);
 error_estimate_rush_eve=zeros(length(t_rush_eve),1);
 error_estimate_off_eve=zeros(length(t_off_eve),1);


figure(1)
% 
plot (t_rush_mor, x_rush_mor(:,2),'-r',t_rush_mor,x_rush_mor(:,6),'-b' ,t_off_mor, x_off_mor(:,2),'-r', t_off_mor, x_off_mor(:,6),'-b' , ...
    t_rush_eve, x_rush_eve(:,2),'-r' , t_rush_eve, x_rush_eve(:,6), '-b' ,t_off_eve,x_off_eve(:,2),'-r' ,t_off_eve, x_off_eve(:,6),'-b' , 'LineWidth',2)
legend('community','Hub')
xlabel('Time (Hour)')
ylabel(' Number of Exposed Individuals within the Hub')


% ----------- Dataset-----------

function stop = myoutput(pa_estimate,optimvalues,state);
        stop = false;
        if isequal(state,'iter')
          history = [history; pa_estimate];
        end

end
%% Confidence Interval  Process
history(:,11)
 N_1=size(history(:,1),1);
 N_Sample_1=size(history(3108:N_1,1),1);

  N_2=size(history(:,2),1);
 N_Sample_2=size(history(3107:N_2,2),1);

  N_3=size(history(:,3),1);
 N_Sample_3=size(history(3107:N_3,3),1);

 N_4=size(history(:,4),1);
 N_Sample_4=size(history(3107:N_4,4),1);

 N_5=size(history(:,5),1);
 N_Sample_5=size(history(3107:N_5,5),1);

 N_6=size(history(:,6),1);
 N_Sample_6=size(history(3107:N_6,6),1);

 N_7=size(history(:,7),1);
 N_Sample_7=size(history(3107:N_7,7),1);

 N_8=size(history(:,8),1);
 N_Sample_8=size(history(3107:N_8,8),1);

 N_9=size(history(:,9),1);
 N_Sample_9=size(history(3107:N_9,9),1);

 N_10=size(history(:,10),1);
 N_Sample_10=size(history(3107:N_10,10),1);

 N_11=size(history(:,11),1);
 N_Sample_11=size(history(3107:N_11,11),1);

 N_12=size(history(:,12),1);
 N_Sample_12=size(history(3107:N_12,12),1);

 N_13=size(history(:,13),1);
 N_Sample_13=size(history(3107:N_13,13),1);

 N_14=size(history(:,14),1);
 N_Sample_14=size(history(3107:N_14,14),1);

N_15=size(history(:,15),1);
N_Sample_15=size(history(3107:N_15,15),1);


N_16=size(history(:,16),1);
N_Sample_16=size(history(3107:N_16,16),1);

N_17=size(history(:,17),1);
N_Sample_17=size(history(3107:N_17,17),1);


N_18=size(history(:,18),1);
N_Sample_18=size(history(3107:N_18,18),1);


N_19=size(history(:,19),1);
N_Sample_19=size(history(3107:N_19,19),1);

N_20=size(history(:,20),1);
N_Sample_20=size(history(3107:N_20,20),1);






%% Mean
% 
 mu_beta_c=mean(history(3108:N_1,1));
 mu_beta_h_rush_mor=mean(history(3107:N_2,2));
 mu_beta_h_off_mor=mean(history(3107:N_3,3));
 mu_beta_h_rush_eve=mean(history(3107:N_4,4));
 mu_beta_h_off_eve=mean(history(3107:N_5,5));
 mu_alpha_rush_mor=mean(history(3107:N_6,6));
 mu_alpha_off_mor=mean(history(3107:N_7,7));
 mu_alpha_rush_eve=mean(history(3107:N_8,8));
 mu_alpha_off_eve=mean(history(3107:N_9,9));
 mu_gamma_rush_mor=mean(history(3107:N_10,10));
 mu_gamma_off_mor=mean(history(3107:N_11,11));
 mu_gamma_rush_eve=mean(history(3107:N_12,12));
 mu_gamma_off_eve=mean(history(3107:N_13,13));
 mu_k_off=mean(history(3107:N_14,14));
 mu_k_rush=mean(history(3107:N_15,15));
 mu_G_off=mean(history(3107:N_16,16));
 mu_G_rush=mean(history(3107:N_17,17));
 mu_N=mean(history(3107:N_18,18));
 mu_I=mean(history(3107:N_19,19));
 mu_p=mean(history(3107:N_20,20));
 
%% SD
 SD_beta_c=std(history(3108:N_1,1));
 SD_beta_h_rush_mor=std(history(3107:N_2,2));
 SD_beta_h_off_mor=std(history(3107:N_3,3));
 SD_beta_h_rush_eve=std(history(3107:N_4,4));
 SD_beta_h_off_eve=std(history(3107:N_5,5));
 SD_alpha_rush_mor=std(history(3107:N_6,6));
 SD_alpha_off_mor=std(history(3107:N_7,7));
 SD_alpha_rush_eve=std(history(3107:N_8,8));
 SD_alpha_off_eve=std(history(3107:N_9,9));
 SD_gamma_rush_mor=std(history(3107:N_10,10));
 SD_gamma_off_mor=std(history(3107:N_11,11));
 SD_gamma_rush_eve=std(history(3107:N_12,12));
 SD_gamma_off_eve=std(history(3107:N_13,13));
 SD_k_off=std(history(3107:N_14,14));
 SD_k_rush=std(history(3107:N_15,15));
 SD_G_off=std(history(3107:N_16,16));
 SD_G_rush=std(history(3107:N_17,17));
 SD_N=std(history(3107:N_18,18));
 SD_I=std(history(3107:N_19,19));
 SD_p=std(history(3107:N_20,20));
%% CI

CI_beta_c_left=mu_beta_c-1.96*(SD_beta_c/sqrt(N_Sample_1));
CI_beta_c_right=mu_beta_c+1.96*(SD_beta_c/sqrt(N_Sample_1));

CI_beta_h_rush_mor_left=mu_beta_h_rush_mor-1.96*(SD_beta_h_rush_mor/sqrt(N_Sample_2));
 CI_beta_h_rush_mor_right=mu_beta_h_rush_mor+1.96*(SD_beta_h_rush_mor/sqrt(N_Sample_2));

 CI_beta_h_off_mor_left=mu_beta_h_off_mor-1.96*(SD_beta_h_off_mor/sqrt(N_Sample_3));
 CI_beta_h_off_mor_right=mu_beta_h_off_mor+1.96*(SD_beta_h_off_mor/sqrt(N_Sample_3));


CI_beta_h_rush_eve_left=mu_beta_h_rush_eve-1.96*(SD_beta_h_rush_eve/sqrt(N_Sample_4));
CI_beta_h_rush_eve_right=mu_beta_h_rush_eve+1.96*(SD_beta_h_rush_eve/sqrt(N_Sample_4));

CI_beta_h_off_eve_left=mu_beta_h_off_eve-1.96*(SD_beta_h_off_eve/sqrt(N_Sample_5));
CI_beta_h_off_eve_right=mu_beta_h_off_eve+1.96*(SD_beta_h_off_eve/sqrt(N_Sample_5));


CI_alpha_rush_mor_left=mu_alpha_rush_mor-1.96*(SD_alpha_rush_mor/sqrt(N_Sample_6));
CI_alpha_rush_mor_right=mu_alpha_rush_mor+1.96*(SD_alpha_rush_mor/sqrt(N_Sample_6));


CI_alpha_off_mor_left=mu_alpha_off_mor-1.96*(SD_alpha_off_mor/sqrt(N_Sample_7));
CI_alpha_off_mor_right=mu_alpha_off_mor+1.96*(SD_alpha_off_mor/sqrt(N_Sample_7));

CI_alpha_rush_eve_left=mu_alpha_rush_eve-1.96*(SD_alpha_rush_eve/sqrt(N_Sample_8));
CI_alpha_rush_eve_right=mu_alpha_rush_eve+1.96*(SD_alpha_rush_eve/sqrt(N_Sample_8));

CI_alpha_off_eve_left=mu_alpha_off_eve-1.96*(SD_alpha_off_eve/sqrt(N_Sample_9));
CI_alpha_off_eve_right=mu_alpha_off_eve+1.96*(SD_alpha_off_eve/sqrt(N_Sample_9));

CI_gamma_rush_mor_left=mu_gamma_rush_mor-1.96*(SD_gamma_rush_mor/sqrt(N_Sample_10));
CI_gamma_rush_mor_right=mu_gamma_rush_mor+1.96*(SD_gamma_rush_mor/sqrt(N_Sample_10));

CI_gamma_off_mor_left=mu_gamma_off_mor-1.96*(SD_gamma_off_mor/sqrt(N_Sample_11));
CI_gamma_off_mor_right=mu_gamma_off_mor+1.96*(SD_gamma_off_mor/sqrt(N_Sample_11));

CI_gamma_rush_eve_left=mu_gamma_rush_eve-1.96*(SD_gamma_rush_eve/sqrt(N_Sample_12));
CI_gamma_rush_eve_right=mu_gamma_rush_eve+1.96*(SD_gamma_rush_eve/sqrt(N_Sample_12));

CI_gamma_off_eve_left=mu_gamma_off_eve-1.96*(SD_gamma_off_eve/sqrt(N_Sample_13));
CI_gamma_off_eve_right=mu_gamma_off_eve+1.96*(SD_gamma_off_eve/sqrt(N_Sample_13));

CI_k_off_left=mu_k_off-1.96*(SD_k_off/sqrt(N_Sample_14));
CI_k_off_right=mu_k_off+1.96*(SD_k_off/sqrt(N_Sample_14));

CI_k_rush_left=mu_k_rush-1.96*(SD_k_rush/sqrt(N_Sample_15));
CI_k_rush_right=mu_k_rush+1.96*(SD_k_rush/sqrt(N_Sample_15));

CI_G_off_left=mu_G_off-1.96*(SD_G_off/sqrt(N_Sample_16));
CI_G_off_right=mu_G_off+1.96*(SD_G_off/sqrt(N_Sample_16));

CI_G_rush_left=mu_G_rush-1.96*(SD_G_rush/sqrt(N_Sample_17));
CI_G_rush_right=mu_G_rush+1.96*(SD_G_rush/sqrt(N_Sample_17));


CI_N_left=mu_N-1.96*(SD_N/sqrt(N_Sample_18));
CI_N_right=mu_N+1.96*(SD_N/sqrt(N_Sample_18));

CI_I_left=mu_I-1.96*(SD_I/sqrt(N_Sample_19));
CI_I_right=mu_I+1.96*(SD_I/sqrt(N_Sample_19));

CI_p_left=mu_p-1.96*(SD_p/sqrt(N_Sample_20));
CI_p_right=mu_p+1.96*(SD_p/sqrt(N_Sample_20));










%----------- Model Function----------------

    function [error_estimate  x_rush_mor x_off_mor x_rush_eve x_off_eve]=fmin_estimatettc_one(param)

beta_c=param(1);            % transmission rate of the community
beta_h_rush_mor=param(2);   % transmission rate of the hub in morning rush hours
beta_h_off_mor=param(3);    % transmission rate of the hub in morning off-peak hours
beta_h_rush_eve=param(4)    % transmission rate of the hub in evening rush hours
beta_h_off_eve=param(5);    % transmission rate of the hub in evening off-peak hours
alpha_rush_mor=param(6);    % inflow rate  from community to the hub in morning rush hours
alpha_off_mor=param(7);     % inflow rate  from community to the hub in morning off-peak hours
alpha_rush_eve=param(8);    % inflow rate  from community to the hub in evening rush hours
alpha_off_eve=param(9);     % inflow rate  from community to the hub in evening off-peak hours
gamma_rush_mor=param(10);    %  outflow rate  from hub to the community in morning rush hours
gamma_off_mor=param(11);     % outflow rate  from hub to the community in morning off-peak hours
gamma_rush_eve=param(12);    %  outflow rate  from hub to the community in evening rush hours
gamma_off_eve=param(13);  
k_off=param(14);         % satuartion number for the force of infection of the community  in off-peak hours
k_rush=param(15);        % satuartion number for the force of infection of the community  in rush hours 
G_off=param(16);        % satuartion number for the force of infection of the hub  in off-peak hours
G_rush=param(17);       % satuartion number for the force of infection of the hub  in rush hours
 N=param(18);           % Number of the total population of the community
 I_0c=param(19);       % Initial data of infectees in the community
 p=param(20);            % Percentge of getting infectious

epsilon_c= 0.0116;
epsilon_h= 0.00116;
delta_c=0.0041;
delta_h=0.00041;



%-------------------- Initial Data-----

N_h_Finch=61323;    % Number of the indviduals visiting Finch station daily
E_0c=0;              % initial data for number of the exposed individuals in the community
R_0c=0;               % initial data for number of the recovered individuals in the community
S_0c=N-(I_0c+E_0c);  %initial data for number of the scuceptible individuals in the community
S_0h=389;              %initial data for number of the scuceptible individuals in the hub
E_0h=0;              % initial data for number of the exposed individuals in the hub
I_0h=0;              % initial data for number of the infected individuals in the hub
R_0h=0;              % initial data for number of the recovered individuals in the hub
IC_estimate=[S_0c,E_0c,I_0c,R_0c,S_0h,E_0h,I_0h,R_0h]; % Vector of the initial data

%---------------- Data Fitting -------

%% Rush-peak hour Morning: 6AM-9AM
t_rush_mor=linspace(0,3,4);
IC_estimate_rush_mor=IC_estimate;
p_rush_mor=[beta_c,beta_h_rush_mor,alpha_rush_mor,gamma_rush_mor,k_rush G_rush,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_mor]=ode45(@(t,x_rush_mor)ttccase_one(t,x_rush_mor,p_rush_mor),t_rush_mor,IC_estimate_rush_mor,op);


%% Off-peak hour midday: 9AM-3PM
t_off_mor=linspace(3,9,7);

IC_estimate_off_mor=[x_rush_mor(4,1),x_rush_mor(4,2),x_rush_mor(4,3),x_rush_mor(4,4),x_rush_mor(4,5),x_rush_mor(4,6),x_rush_mor(4,7),x_rush_mor(4,8)]; % Vector of the initial data in off-peak hours
p_off_mor=[beta_c,beta_h_off_mor,alpha_off_mor,gamma_off_mor,k_off, G_off,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_mor]=ode45(@(t,x_off_mor)ttccase_one(t,x_off_mor,p_off_mor),t_off_mor,IC_estimate_off_mor,op);

%% Rush hour evening: 3AM-7AM
t_rush_eve=linspace(9,13,5);

IC_estimate_rush_eve=[x_off_mor(7,1),x_off_mor(7,2),x_off_mor(7,3),x_off_mor(7,4),x_off_mor(7,5),x_off_mor(7,6),x_off_mor(7,7),x_off_mor(7,8)]% Vector of the initial data in rush hours
p_rush_eve=[beta_c,beta_h_rush_eve,alpha_rush_eve,gamma_rush_eve,k_rush, G_rush,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_eve]=ode45(@(t,x_rush_eve)ttccase_one(t,x_rush_eve,p_rush_eve),t_rush_eve,IC_estimate_rush_eve,op);

%% Off-peak hour evening: 7AM-2AM
t_off_eve=linspace(13,20,8);
IC_estimate_off_eve=[x_rush_eve(5,1),x_rush_eve(5,2),x_rush_eve(5,3),x_rush_eve(5,4),x_rush_eve(5,5),x_rush_eve(5,6),x_rush_eve(5,7),x_rush_eve(5,8)]; % Vector of the initial data in off-peak hours
p_off_eve=[beta_c,beta_h_off_eve,alpha_off_eve,gamma_off_eve,k_off, G_off,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_eve]=ode45(@(t,x_off_eve)ttccase_one(t,x_off_eve,p_off_eve),t_off_eve,IC_estimate_off_eve,op);
%% Result

% Cumulative numbers of infectees in the community
I_c_cum_rush_mor=cumsum(x_rush_mor(:,3));
I_c_cum_off_mor=cumsum(x_off_mor(:,3));
 I_c_cum_rush_eve=cumsum(x_rush_eve(:,3));
 I_c_cum_off_eve=cumsum(x_off_eve(:,3));

I_c_cum_end=I_c_cum_rush_mor(end)-I_0c+I_c_cum_off_mor(end)-x_off_mor(1,3)+I_c_cum_rush_eve(end)-...
x_rush_eve(1,3)+I_c_cum_off_eve(end)-x_off_eve(1,3) ;  

 % Cumulative numbers of exposed individuals in the community

E_c_cum_rush_mor=cumsum(x_rush_mor(:,2));
E_c_cum_off_mor=cumsum(x_off_mor(:,2));
 E_c_cum_rush_eve=cumsum(x_rush_eve(:,2));
 E_c_cum_off_eve=cumsum(x_off_eve(:,2));

E_c_cum_end=E_c_cum_rush_mor(end)-E_0c+E_c_cum_off_mor(end)-x_off_mor(1,2)+E_c_cum_rush_eve(end)-x_rush_eve(1,2)+...
    E_c_cum_off_eve(end)-x_off_eve(1,2) ; 

 % Cumulative numbers of exposed individuals in the hub
E_h_cum_rush_mor=cumsum(x_rush_mor(:,6));
E_h_cum_off_mor=cumsum(x_off_mor(:,6));
 E_h_cum_rush_eve=cumsum(x_rush_eve(:,6));
 E_h_cum_off_eve=cumsum(x_off_eve(:,6));

E_h_cum_end=E_h_cum_rush_mor(end)+E_h_cum_off_mor(end)-x_off_mor(1,6)+E_h_cum_rush_eve(end)-x_rush_eve(1,6)+...
    E_h_cum_off_eve(end)-x_off_eve(1,6);

% Cumulative numbers of infeted individuals in the hub
 
I_h_cum_rush_mor=cumsum(x_rush_mor(:,7));
I_h_cum_off_mor=cumsum(x_off_mor(:,7));
 I_h_cum_rush_eve=cumsum(x_rush_eve(:,7));
 I_h_cum_off_eve=cumsum(x_off_eve(:,7));

I_h_cum_end=I_h_cum_rush_mor(end)+I_h_cum_off_mor(end)-x_off_mor(1,7)+I_h_cum_rush_eve(end)-x_rush_eve(1,7)+...
    I_h_cum_off_eve(end)-x_off_eve(1,7) ;  

% Cumulative numbers of infeted individuals in the hub

R_h_cum_rush_mor=cumsum(x_rush_mor(:,8));
R_h_cum_off_mor=cumsum(x_off_mor(:,8));
 R_h_cum_rush_eve=cumsum(x_rush_eve(:,8));
 R_h_cum_off_eve=cumsum(x_off_eve(:,8));


R_h_cum_end=R_h_cum_rush_mor(end)+R_h_cum_off_mor(end)-x_off_mor(1,8)+R_h_cum_rush_eve(end)-x_rush_eve(1,8)+...
    R_h_cum_off_eve(end)-x_off_eve(1,8) ;  


 % Cumulative numbers of scucptible individuals in the hub

S_h_cum_rush_mor=cumsum(x_rush_mor(:,5));
S_h_cum_rush_mor(end);
S_h_cum_off_mor=cumsum(x_off_mor(:,5));
 S_h_cum_rush_eve=cumsum(x_rush_eve(:,5));
 S_h_cum_off_eve=cumsum(x_off_eve(:,5));

S_h_cum_end=S_h_cum_rush_mor(end)+S_h_cum_off_mor(end)-x_off_mor(1,5)+S_h_cum_rush_eve(end)-x_rush_eve(1,5)+...
    S_h_cum_off_eve(end)-x_off_eve(1,5);

 % Cumulative numbers of  individuals visiting the hub (model-based)
N_h_cum=S_h_cum_end+E_h_cum_end+I_h_cum_end+R_h_cum_end; 

 
for i=1:length(t_rush_mor)
error_estimate_rush_mor(i)=(x_rush_mor(i,5)+x_rush_mor(i,6)+x_rush_mor(i,7)+x_rush_mor(i,8)-Data(i))^2;
x_rush_total_mor_h(i)=x_rush_mor(i,5)+x_rush_mor(i,6)+x_rush_mor(i,7)+x_rush_mor(i,8);
x_rush_total_mor_c(i)=x_rush_mor(i,1)+x_rush_mor(i,2)+x_rush_mor(i,3)+x_rush_mor(i,4);

end
for i=2:length(t_off_mor)
    x_off_total_mor_h(1)=x_rush_total_mor_h(end);
     x_off_total_mor_c(1)=x_rush_total_mor_c(end);
error_estimate_off_mor(i)=(x_off_mor(i,5)+x_off_mor(i,6)+x_off_mor(i,7)+x_off_mor(i,8)-Data(i+3))^2;
x_off_total_mor_h(i)=x_off_mor(i,5)+x_off_mor(i,6)+x_off_mor(i,7)+x_off_mor(i,8);
x_off_total_mor_c(i)=x_off_mor(i,1)+x_off_mor(i,2)+x_off_mor(i,3)+x_off_mor(i,4);


end
for i=2:length(t_rush_eve)
    x_rush_total_eve_h(1)=x_off_total_mor_h(end);
    x_rush_total_eve_c(1)=x_off_total_mor_c(end);
error_estimate_rush_eve(i)=(x_rush_eve(i,5)+x_rush_eve(i,6)+x_rush_eve(i,7)+x_rush_eve(i,8)-Data(i+9))^2;
x_rush_total_eve_h(i)=x_rush_eve(i,5)+x_rush_eve(i,6)+x_rush_eve(i,7)+x_rush_eve(i,8);
x_rush_total_eve_c(i)=x_rush_eve(i,1)+x_rush_eve(i,2)+x_rush_eve(i,3)+x_rush_eve(i,4);

end
for i=2:length(t_off_eve)
    x_off_total_eve_h(1)=x_rush_total_eve_h(end);
     x_off_total_eve_c(1)=x_rush_total_eve_c(end);
error_estimate_off_eve(i)=(x_off_eve(i,5)+x_off_eve(i,6)+x_off_eve(i,7)+x_off_eve(i,8)-Data(i+13))^2;
x_off_total_eve_h(i)=x_off_eve(i,5)+x_off_eve(i,6)+x_off_eve(i,7)+x_off_eve(i,8);
x_off_total_eve_c(i)=x_off_eve(i,1)+x_off_eve(i,2)+x_off_eve(i,3)+x_off_eve(i,4);

end


EE=[error_estimate_rush_mor'; error_estimate_off_mor';error_estimate_rush_eve';error_estimate_off_eve'];

EE_cum=cumsum(EE);
R0_m=R0_estimate(param)
error_estimate_rush_eve';
error_estimate_off_mor';
N_h_cum;
x_rush_total_mor_h
x_rush_total_mor_c;
x_off_total_mor_h
x_off_total_mor_c;
x_rush_total_eve_h
x_rush_total_eve_c;
x_off_total_eve_h
x_off_total_eve_c;

t_total=linspace(0,20,21);
Data_cum=cumsum(Data);
x_rush_cum_total_mor_h=cumsum(x_rush_total_mor_h);
x_rush_cum_total_mor_c=cumsum(x_rush_total_mor_c);

x_off_cum_total_mor_h=cumsum(x_off_total_mor_h);
x_off_cum_total_mor_c=cumsum(x_off_total_mor_c);

x_rush_cum_total_eve_h=cumsum(x_rush_total_eve_h);
x_rush_cum_total_eve_c=cumsum(x_rush_total_eve_c);

x_off_cum_total_eve_h=cumsum(x_off_total_eve_h);
x_off_cum_total_eve_c=cumsum(x_off_total_eve_c);

N_h_total=x_rush_cum_total_mor_h(end)+x_off_cum_total_mor_h(end)-x_off_cum_total_mor_h(1)+...
   x_rush_cum_total_eve_h(end)- x_rush_cum_total_eve_h(1)+x_off_cum_total_eve_h(end)-x_off_cum_total_eve_h(1)

error_estimate= EE_cum(end)+(R0_m-8.2)^2+(N_h_total-N_h_Finch)^2;

figure(2)
 plot(t_rush_mor,x_rush_total_mor_h,'-b', t_off_mor,x_off_total_mor_h,'-b', ...
     t_rush_eve,x_rush_total_eve_h,'-b',t_off_eve,x_off_total_eve_h,'-b', ...
     t_total, Data,'k .', 'markersize', 12,'LineWidth',2)
 xlabel('Time (Hour)')
ylabel(' Number of Individuals Visiting Finch Station in a Single Day')

    end
Ex_cum_ode_rush_mor=[x_rush_mor(:,2)]
Ex_cum_ode_off_mor=[ x_off_mor(:,2)]
Ex_cum_ode_rush_eve=[x_rush_eve(:,2)]
Ex_cum_ode_off_eve=[ x_off_eve(:,2)]
Ex_hub_ode_rush_mor=[x_rush_mor(:,6)]
Ex_hub_ode_off_mor=[ x_off_mor(:,6)]
Ex_hub_ode_rush_eve=[x_rush_eve(:,6)]
Ex_hub_ode_off_eve=[ x_off_eve(:,6)]

end


