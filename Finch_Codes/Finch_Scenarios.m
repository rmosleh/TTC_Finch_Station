function Finch_Scenarios

format long

%--------------- Parameters-----------

beta_c=0.0050038755;               % transmission rate of the community
beta_h_rush_mor=0.1320834617;      % transmission rate of the hub in morning rush hours 
beta_h_off_mor=0.0476675483   ;    % transmission rate of the hub in morning off-peakh hours
beta_h_rush_eve=0.8197465080;     % transmission rate of the hub in evening rush hours
beta_h_off_eve=0.0968375423  ;    % transmission rate of the hub in evening off-peak hours
alpha_rush_mor=0.0275091431;      % inflow rate  from community to the hub in morning rush hours
alpha_off_mor=0.0196874190;       % inflow rate  from community to the hub in morning off-peak hours
alpha_rush_eve=0.0354639424;    % inflow rate  from community to the hub in evening rush hours
alpha_off_eve=0.0033932538  ;     % inflow rate  from community to the hub in evening off-peak hours
gamma_rush_mor=0.3342648200;    %  outflow rate  from hub to the community in morning rush hours
gamma_off_mor=0.8172129389 ;     % outflow rate  from hub to the community in morning off-peak hours
gamma_rush_eve=0.5930270081;    %  outflow rate  from hub to the community in evening rush hours
gamma_off_eve=0.4444764165;     % outflow rate  from hub to the community in evening ff-peak hours
k_off=62242.7236876310;         % satuartion number for the force of infection of the community  in off-peak hours
k_rush=139.8166578181;          % satuartion number for the force of infection of the community  in rush hours 
G_off=36339.5314556524;        % satuartion number for the force of infection of the hub  in off-peak hours
G_rush=788210.2232412722;       % satuartion number for the force of infection of the hub  in rush hours
 p=0.0265663442;               %rate of carring infectious
 epsilon_c= 0.0116;
epsilon_h= 0.00116;
delta_c=0.0019;
delta_h=0.00019;



%-------------------- Initial Data-----

 N_h_Finch=61323;          % actual number of the indviduals visiting Finch station daily
 N=100676.3918036157  ;           % Estimated number of the total population of the community
E_0c=0;                         % initial data for number of the exposed individuals in the community
I_0c=4.8599130629 ;              % estimated initial data for number of the infected individuals in the community
S_0c=N-(I_0c+E_0c);              %initial data for scuceptible individuals in the community
R_0c=0;                       %initial data for number of the recovered individuals in the community
S_0h=389;                     %initial data for number of the scuceptible individuals in the hub
E_0h=0;                       % initial data for number of the exposed individuals in the hub
I_0h=0;                       %initial data for number of the infectious individuals in the hub
R_0h=0;                       % initial data for number of the recovered individuals in the hub
IC_estimate=[S_0c,E_0c,I_0c,R_0c,S_0h,E_0h,I_0h,R_0h]; % Vector of the initial data

%% Transmission rate for the community 

s=[1, 0.8, 0.5, 0.3, 0.1];
c=['r','b','k','g','m'];
for i=1:length(s)

    %% Rush-peak hour Morning: 6AM-9AM
t_rush_1=linspace(0,3,4);
IC_estimate_rush_1=IC_estimate;
p_rush_1=[beta_c * s(i),beta_h_rush_mor ,alpha_rush_mor ,gamma_rush_mor,k_rush G_rush,p,epsilon_c,epsilon_h,delta_c,delta_h];

op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_1]=ode45(@(t,x_rush_1)ttccase_one(t,x_rush_1,p_rush_1),t_rush_1,IC_estimate_rush_1,op);


%% Off-peak hour midday: 9AM-3PM
t_off_1=linspace(3,9,7);

IC_estimate_off_1=[x_rush_1(4,1),x_rush_1(4,2),x_rush_1(4,3),x_rush_1(4,4),x_rush_1(4,5),x_rush_1(4,6),x_rush_1(4,7),x_rush_1(4,8)]; % Vector of the initial data in off-peak hours
p_off_1=[beta_c * s(i),beta_h_off_mor ,alpha_off_mor ,gamma_off_mor,k_off, G_off,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_1]=ode45(@(t,x_off_1)ttccase_one(t,x_off_1,p_off_1),t_off_1,IC_estimate_off_1,op);

%% Rush hour evening: 3AM-7AM
t_rush_2=linspace(9,13,5);

IC_estimate_rush_2=[x_off_1(7,1),x_off_1(7,2),x_off_1(7,3),x_off_1(7,4),x_off_1(7,5),x_off_1(7,6),x_off_1(7,7),x_off_1(7,8)];% Vector of the initial data in rush hours
p_rush_2=[beta_c * s(i),beta_h_rush_eve ,alpha_rush_eve ,gamma_rush_eve,k_rush, G_rush,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_2]=ode45(@(t,x_rush_2)ttccase_one(t,x_rush_2,p_rush_2),t_rush_2,IC_estimate_rush_2,op);

%% Off-peak hour evening: 7AM-2AM
t_off_2=linspace(13,20,8);
IC_estimate_off_2=[x_rush_2(5,1),x_rush_2(5,2),x_rush_2(5,3),x_rush_2(5,4),x_rush_2(5,5),x_rush_2(5,6),x_rush_2(5,7),x_rush_2(5,8)]; % Vector of the initial data in off-peak hours
p_off_2=[beta_c * s(i),beta_h_off_eve,alpha_off_eve,gamma_off_eve,k_off, G_off,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_2]=ode45(@(t,x_off_2)ttccase_one(t,x_off_2,p_off_2),t_off_2,IC_estimate_off_2,op);

figure(1)


plot (t_rush_1, x_rush_1(:,2), c(i), t_off_1, x_off_1(:,2), c(i),t_rush_2, x_rush_2(:,2),c(i),t_off_2,x_off_2(:,2),c(i),'LineWidth',2)
 hold on


%legend({' Baseline ',' 20% decline ','50% decline ','70% decline ','90% decline '} )
xlabel('Time (Hour)')
ylabel(' Number of Exposed Individuals within the Community')
title(' Impacts of Community Transmission Rates')

figure(2)

plot (t_rush_1, x_rush_1(:,6),c(i), t_off_1, x_off_1(:,6),c(i), t_rush_2, x_rush_2(:,6),c(i),t_off_2,x_off_2(:,6),c(i),'LineWidth',2)
hold on 
%legend({'Baseline','20% decline','50% decline','70% decline','90% decline'} )
xlabel('Time (Hour)')
ylabel(' Number of Exposed Individuals within the Hub')
title('Impacts of Community Transmission Rates')
end

%% Transmission rate for the hub

for i=1:length(s)

    %% Rush-peak hour Morning: 6AM-9AM
t_rush_1=linspace(0,3,4);
IC_estimate_rush_1=IC_estimate;
p_rush_1=[beta_c ,beta_h_rush_mor* s(i)  ,alpha_rush_mor,gamma_rush_mor,k_rush G_rush,p,epsilon_c,epsilon_h,delta_c,delta_h];

op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_1]=ode45(@(t,x_rush_1)ttccase_one(t,x_rush_1,p_rush_1),t_rush_1,IC_estimate_rush_1,op);






%% Off-peak hour midday: 9AM-3PM
t_off_1=linspace(3,9,7);

IC_estimate_off_1=[x_rush_1(4,1),x_rush_1(4,2),x_rush_1(4,3),x_rush_1(4,4),x_rush_1(4,5),x_rush_1(4,6),x_rush_1(4,7),x_rush_1(4,8)]; % Vector of the initial data in off-peak hours
p_off_1=[beta_c,beta_h_off_mor* s(i),alpha_off_mor ,gamma_off_mor,k_off, G_off,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_1]=ode45(@(t,x_off_1)ttccase_one(t,x_off_1,p_off_1),t_off_1,IC_estimate_off_1,op);

%% Rush hour evening: 3AM-7AM
t_rush_2=linspace(9,13,5);

IC_estimate_rush_2=[x_off_1(7,1),x_off_1(7,2),x_off_1(7,3),x_off_1(7,4),x_off_1(7,5),x_off_1(7,6),x_off_1(7,7),x_off_1(7,8)];% Vector of the initial data in rush hours
p_rush_2=[beta_c ,beta_h_rush_eve* s(i),alpha_rush_eve ,gamma_rush_eve,k_rush, G_rush,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_2]=ode45(@(t,x_rush_2)ttccase_one(t,x_rush_2,p_rush_2),t_rush_2,IC_estimate_rush_2,op);

%% Off-peak hour evening: 7AM-2AM
t_off_2=linspace(13,20,8);
IC_estimate_off_2=[x_rush_2(5,1),x_rush_2(5,2),x_rush_2(5,3),x_rush_2(5,4),x_rush_2(5,5),x_rush_2(5,6),x_rush_2(5,7),x_rush_2(5,8)]; % Vector of the initial data in off-peak hours
p_off_2=[beta_c ,beta_h_off_eve* s(i),alpha_off_eve,gamma_off_eve,k_off, G_off,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_2]=ode45(@(t,x_off_2)ttccase_one(t,x_off_2,p_off_2),t_off_2,IC_estimate_off_2,op);

figure(3)

plot (t_rush_1, x_rush_1(:,2), c(i), t_off_1, x_off_1(:,2), c(i),t_rush_2, x_rush_2(:,2),c(i),t_off_2,x_off_2(:,2),c(i),'LineWidth',2)
hold on 
%legend('Baseline','20% decline','50% decline','70% decline','90% decline' )
xlabel('Time (Hour)')
ylabel(' Number of Exposed Individuals within the Community')
title('Impacts of Hub Transmission Rates')

figure(4)

plot (t_rush_1, x_rush_1(:,6),c(i), t_off_1, x_off_1(:,6),c(i), t_rush_2, x_rush_2(:,6),c(i),t_off_2,x_off_2(:,6),c(i),'LineWidth',2)
hold on 
%legend('Baseline','20% decline','50% decline','70% decline','90% decline' )
xlabel('Time (Hour)')
ylabel(' Number of Exposed Individuals within the Hub')
title('Impacts of Hub Transmission Rates')
end
% 
%% Inflow rate

for i=1:length(s)

    %% Rush-peak hour Morning: 6AM-9AM
t_rush_1=linspace(0,3,4);
IC_estimate_rush_1=IC_estimate;
p_rush_1=[beta_c ,beta_h_rush_mor,alpha_rush_mor * s(i) ,gamma_rush_mor,k_rush G_rush,p,epsilon_c,epsilon_h,delta_c,delta_h];

op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_1]=ode45(@(t,x_rush_1)ttccase_one(t,x_rush_1,p_rush_1),t_rush_1,IC_estimate_rush_1,op);






%% Off-peak hour midday: 9AM-3PM
t_off_1=linspace(3,9,7);

IC_estimate_off_1=[x_rush_1(4,1),x_rush_1(4,2),x_rush_1(4,3),x_rush_1(4,4),x_rush_1(4,5),x_rush_1(4,6),x_rush_1(4,7),x_rush_1(4,8)]; % Vector of the initial data in off-peak hours
p_off_1=[beta_c,beta_h_off_mor  ,alpha_off_mor * s(i) ,gamma_off_mor,k_off, G_off,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_1]=ode45(@(t,x_off_1)ttccase_one(t,x_off_1,p_off_1),t_off_1,IC_estimate_off_1,op);

%% Rush hour evening: 3AM-7AM
t_rush_2=linspace(9,13,5);

IC_estimate_rush_2=[x_off_1(7,1),x_off_1(7,2),x_off_1(7,3),x_off_1(7,4),x_off_1(7,5),x_off_1(7,6),x_off_1(7,7),x_off_1(7,8)];% Vector of the initial data in rush hours
p_rush_2=[beta_c ,beta_h_rush_eve,alpha_rush_eve * s(i) ,gamma_rush_eve,k_rush, G_rush,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_2]=ode45(@(t,x_rush_2)ttccase_one(t,x_rush_2,p_rush_2),t_rush_2,IC_estimate_rush_2,op);

%% Off-peak hour evening: 7AM-2AM
t_off_2=linspace(13,20,8);
IC_estimate_off_2=[x_rush_2(5,1),x_rush_2(5,2),x_rush_2(5,3),x_rush_2(5,4),x_rush_2(5,5),x_rush_2(5,6),x_rush_2(5,7),x_rush_2(5,8)]; % Vector of the initial data in off-peak hours
p_off_2=[beta_c,beta_h_off_eve,alpha_off_eve* s(i),gamma_off_eve,k_off, G_off,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_2]=ode45(@(t,x_off_2)ttccase_one(t,x_off_2,p_off_2),t_off_2,IC_estimate_off_2,op);

figure(5)

plot (t_rush_1, x_rush_1(:,2), c(i), t_off_1, x_off_1(:,2), c(i),t_rush_2, x_rush_2(:,2),c(i),t_off_2,x_off_2(:,2),c(i),'LineWidth',2)
hold on 
%legend('Baseline','20% decline','50% decline','70% decline','90% decline' )
xlabel('Time (Hour)')
ylabel(' Number of Exposed Individuals within the Community')
title('Impacts of Inflow Rates')

figure(6)

plot (t_rush_1, x_rush_1(:,6),c(i), t_off_1, x_off_1(:,6),c(i), t_rush_2, x_rush_2(:,6),c(i),t_off_2,x_off_2(:,6),c(i),'LineWidth',2)
hold on 
%legend('Baseline','20% decline','50% decline','70% decline','90% decline' )
xlabel('Time (Hour)')
ylabel(' Number of Exposed Individuals within the Hub')
title('Impacts of Inflow Rates')
end
% 
%% Outflow rate

for i=1:length(s)

    %% Rush-peak hour Morning: 6AM-9AM
t_rush_1=linspace(0,3,4);
IC_estimate_rush_1=IC_estimate;
p_rush_1=[beta_c ,beta_h_rush_mor  ,alpha_rush_mor  ,gamma_rush_mor * s(i),k_rush G_rush,p,epsilon_c,epsilon_h,delta_c,delta_h];

op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_1]=ode45(@(t,x_rush_1)ttccase_one(t,x_rush_1,p_rush_1),t_rush_1,IC_estimate_rush_1,op);






%% Off-peak hour midday: 9AM-3PM
t_off_1=linspace(3,9,7);

IC_estimate_off_1=[x_rush_1(4,1),x_rush_1(4,2),x_rush_1(4,3),x_rush_1(4,4),x_rush_1(4,5),x_rush_1(4,6),x_rush_1(4,7),x_rush_1(4,8)]; % Vector of the initial data in off-peak hours
p_off_1=[beta_c,beta_h_off_mor  ,alpha_off_mor  ,gamma_off_mor * s(i),k_off, G_off,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_1]=ode45(@(t,x_off_1)ttccase_one(t,x_off_1,p_off_1),t_off_1,IC_estimate_off_1,op);

%% Rush hour evening: 3AM-7AM
t_rush_2=linspace(9,13,5);

IC_estimate_rush_2=[x_off_1(7,1),x_off_1(7,2),x_off_1(7,3),x_off_1(7,4),x_off_1(7,5),x_off_1(7,6),x_off_1(7,7),x_off_1(7,8)];% Vector of the initial data in rush hours
p_rush_2=[beta_c ,beta_h_rush_eve,alpha_rush_eve ,gamma_rush_eve * s(i),k_rush, G_rush,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_2]=ode45(@(t,x_rush_2)ttccase_one(t,x_rush_2,p_rush_2),t_rush_2,IC_estimate_rush_2,op);

%% Off-peak hour evening: 7AM-2AM
t_off_2=linspace(13,20,8);
IC_estimate_off_2=[x_rush_2(5,1),x_rush_2(5,2),x_rush_2(5,3),x_rush_2(5,4),x_rush_2(5,5),x_rush_2(5,6),x_rush_2(5,7),x_rush_2(5,8)]; % Vector of the initial data in off-peak hours
p_off_2=[beta_c,beta_h_off_eve,alpha_off_eve ,gamma_off_eve * s(i),k_off, G_off,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_2]=ode45(@(t,x_off_2)ttccase_one(t,x_off_2,p_off_2),t_off_2,IC_estimate_off_2,op);

figure(7)

plot (t_rush_1, x_rush_1(:,2), c(i), t_off_1, x_off_1(:,2), c(i),t_rush_2, x_rush_2(:,2),c(i),t_off_2,x_off_2(:,2),c(i),'LineWidth',2)
hold on 
%legend('Baseline','20% decline','50% decline','70% decline','90% decline' )
xlabel('Time (Hour)')
ylabel(' Number of Exposed Individuals within the Community')
title('Impacts of Outflow Rates')

figure(8)

plot (t_rush_1, x_rush_1(:,6),c(i), t_off_1, x_off_1(:,6),c(i), t_rush_2, x_rush_2(:,6),c(i),t_off_2,x_off_2(:,6),c(i),'LineWidth',2)
hold on 
%legend('Baseline','20% decline','50% decline','70% decline','90% decline' )
xlabel('Time (Hour)')
ylabel(' Number of Exposed Individuals within the Hub')
title('Impacts of Outflow Rates')
end
% % 
%% Inflow and Outflow rates

for i=1:length(s)

    %% Rush-peak hour Morning: 6AM-9AM
t_rush_1=linspace(0,3,4);
IC_estimate_rush_1=IC_estimate;
p_rush_1=[beta_c ,beta_h_rush_mor  ,alpha_rush_mor * s(i) ,gamma_rush_mor * s(i),k_rush G_rush,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_1]=ode45(@(t,x_rush_1)ttccase_one(t,x_rush_1,p_rush_1),t_rush_1,IC_estimate_rush_1,op);






%% Off-peak hour midday: 9AM-3PM
t_off_1=linspace(3,9,7);

IC_estimate_off_1=[x_rush_1(4,1),x_rush_1(4,2),x_rush_1(4,3),x_rush_1(4,4),x_rush_1(4,5),x_rush_1(4,6),x_rush_1(4,7),x_rush_1(4,8)]; % Vector of the initial data in off-peak hours
p_off_1=[beta_c,beta_h_off_mor  ,alpha_off_mor * s(i) ,gamma_off_mor * s(i),k_off, G_off,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_1]=ode45(@(t,x_off_1)ttccase_one(t,x_off_1,p_off_1),t_off_1,IC_estimate_off_1,op);

%% Rush hour evening: 3AM-7AM
t_rush_2=linspace(9,13,5);

IC_estimate_rush_2=[x_off_1(7,1),x_off_1(7,2),x_off_1(7,3),x_off_1(7,4),x_off_1(7,5),x_off_1(7,6),x_off_1(7,7),x_off_1(7,8)];% Vector of the initial data in rush hours
p_rush_2=[beta_c ,beta_h_rush_eve,alpha_rush_eve * s(i) ,gamma_rush_eve  * s(i),k_rush, G_rush,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_2]=ode45(@(t,x_rush_2)ttccase_one(t,x_rush_2,p_rush_2),t_rush_2,IC_estimate_rush_2,op);

%% Off-peak hour evening: 7AM-2AM
t_off_2=linspace(13,20,8);
IC_estimate_off_2=[x_rush_2(5,1),x_rush_2(5,2),x_rush_2(5,3),x_rush_2(5,4),x_rush_2(5,5),x_rush_2(5,6),x_rush_2(5,7),x_rush_2(5,8)]; % Vector of the initial data in off-peak hours
p_off_2=[beta_c,beta_h_off_eve,alpha_off_eve * s(i),gamma_off_eve  * s(i),k_off, G_off,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_2]=ode45(@(t,x_off_2)ttccase_one(t,x_off_2,p_off_2),t_off_2,IC_estimate_off_2,op);

figure(9)

plot (t_rush_1, x_rush_1(:,2), c(i), t_off_1, x_off_1(:,2), c(i),t_rush_2, x_rush_2(:,2),c(i),t_off_2,x_off_2(:,2),c(i),'LineWidth',2)
hold on 
%legend('Baseline','20% decline','50% decline','70% decline','90% decline' )
xlabel('Time (Hour)')
ylabel(' Number of Exposed Individuals within the Community')
title('Impacts of Inflow and Outflow Rates')

figure(10)

plot (t_rush_1, x_rush_1(:,6),c(i), t_off_1, x_off_1(:,6),c(i), t_rush_2, x_rush_2(:,6),c(i),t_off_2,x_off_2(:,6),c(i),'LineWidth',2)
hold on 
%legend('Baseline','20% decline','50% decline','70% decline','90% decline' )
xlabel('Time (Hour)')
ylabel(' Number of Exposed Individuals within the Hub')
title('Impacts of Inflow and Outflow Rates')
end
% 
%% Inflow and Outflow rates-Heatmaps

s_1=0.0:0.01:1;
s_2=0.0:0.01:1;
for k=1:length(s_1)
    for j=1:length(s_2)

    %% Rush-peak hour Morning: 6AM-9AM
t_rush_1=linspace(0,3,4);
IC_estimate_rush_1=IC_estimate;
p_rush_1=[beta_c ,beta_h_rush_mor,alpha_rush_mor * (1-s_1(k)) ,gamma_rush_mor * (1-s_2(j)),k_rush G_rush,p,epsilon_c,epsilon_h,delta_c,delta_h];

op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_1]=ode45(@(t,x_rush_1)ttccase_one(t,x_rush_1,p_rush_1),t_rush_1,IC_estimate_rush_1,op);






%% Off-peak hour midday: 9AM-3PM
t_off_1=linspace(3,9,7);

IC_estimate_off_1=[x_rush_1(4,1),x_rush_1(4,2),x_rush_1(4,3),x_rush_1(4,4),x_rush_1(4,5),x_rush_1(4,6),x_rush_1(4,7),x_rush_1(4,8)]; % Vector of the initial data in off-peak hours
p_off_1=[beta_c,beta_h_off_mor,alpha_off_mor * (1-s_1(k)) ,gamma_off_mor * (1-s_2(j)),k_off, G_off,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_1]=ode45(@(t,x_off_1)ttccase_one(t,x_off_1,p_off_1),t_off_1,IC_estimate_off_1,op);

%% Rush hour evening: 3AM-7AM
t_rush_2=linspace(9,13,5);

IC_estimate_rush_2=[x_off_1(7,1),x_off_1(7,2),x_off_1(7,3),x_off_1(7,4),x_off_1(7,5),x_off_1(7,6),x_off_1(7,7),x_off_1(7,8)];% Vector of the initial data in rush hours
p_rush_2=[beta_c ,beta_h_rush_eve,alpha_rush_eve * (1-s_1(k)) ,gamma_rush_eve  * (1-s_2(j)),k_rush, G_rush,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_2]=ode45(@(t,x_rush_2)ttccase_one(t,x_rush_2,p_rush_2),t_rush_2,IC_estimate_rush_2,op);

%% Off-peak hour evening: 7AM-2AM
t_off_2=linspace(13,20,8);
IC_estimate_off_2=[x_rush_2(5,1),x_rush_2(5,2),x_rush_2(5,3),x_rush_2(5,4),x_rush_2(5,5),x_rush_2(5,6),x_rush_2(5,7),x_rush_2(5,8)]; % Vector of the initial data in off-peak hours
p_off_2=[beta_c,beta_h_off_eve,alpha_off_eve * (1-s_1(k)),gamma_off_eve * (1-s_2(j)),k_off, G_off,p,epsilon_c,epsilon_h,delta_c,delta_h];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_2]=ode45(@(t,x_off_2)ttccase_one (t,x_off_2,p_off_2),t_off_2,IC_estimate_off_2,op);

%%



M_com_rush_1=max(x_rush_1(:,2));
M_com_off_1=max(x_off_1(:,2));
M_com_rush_2=max(x_rush_2(:,2));
M_com_off_2=max(x_off_2(:,2));
M_com_total=[M_com_rush_1, M_com_off_1,M_com_rush_2,M_com_off_2];
MM_com(k,j)=max(M_com_total);


M_hub_rush_1=max(x_rush_1(:,6));
M_hub_off_1=max(x_off_1(:,6));
M_hub_rush_2=max(x_rush_2(:,6));
M_hub_off_2=max(x_off_2(:,6));
M_hub_total_1=[M_hub_rush_1, M_hub_off_1 ];
M_hub_total_2=[M_hub_rush_2,M_hub_off_2];
MM_hub_1(k,j)=max(M_hub_total_1);
MM_hub_2(k,j)=max(M_hub_total_2);




 end
end
figure(11)
[X,Y]=meshgrid(s_1,s_2)
 [C,h]=contourf(X,Y,MM_com,50)
 colormap(hot)
 colorbar
 title('Community')
xlabel('Reduction in Outflow Rates')
ylabel(' Reduction in Inflow Rates')
 figure(12)
[X,Y]=meshgrid(s_1,s_2)
 [C,h]=contourf(X,Y,MM_hub_1,50)
 colormap(hot)
 colorbar
 title('Hub - First Peak')
 xlabel('Reduction in Outflow Rates')
ylabel('Reduction in Inflow Rates')

 figure(13)
[X,Y]=meshgrid(s_1,s_2)
 [C,h]=contourf(X,Y,MM_hub_2,50)
 colormap(hot)
 colorbar
 title('Hub - Second Peak')
 xlabel('Reduction in Outflow Rates')
ylabel(' Reduction in Inflow Rates')

       [MAX_MM, ind_max]=max(MM_com,[],"all")
       [row,col,page] = ind2sub(size(MM_com),ind_max)
       [min_MM, ind_min]=min(MM_com,[],"all")
       [row,col,page] = ind2sub(size(MM_com),ind_min)

       [MAX_MM_hub_1, ind_max_hub_1]=max(MM_hub_1,[],"all")
       [row,col,page] = ind2sub(size(MM_hub_1),ind_max_hub_1)
       [min_MM_hub_1, ind_min_hub_1]=min(MM_hub_1,[],"all")
       [row,col,page] = ind2sub(size(MM_hub_1),ind_min_hub_1)

       [MAX_MM_hub_2, ind_max_hub_2]=max(MM_hub_2,[],"all")
       [row,col,page] = ind2sub(size(MM_hub_2),ind_max_hub_2)
       [min_MM_hub_2, ind_min_hub_2]=min(MM_hub_2,[],"all")
       [row,col,page] = ind2sub(size(MM_hub_2),ind_min_hub_2)


