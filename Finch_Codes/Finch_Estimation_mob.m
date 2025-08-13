function [pa_estimate fval history] = Finch_Estimation_mod
format long

%% Estimating Process for the mobility model over a single day with the least squares method

 history=[];
 Data=load('Data_mod.txt');

%-----------Estimate_Parameters_initial value--------


alpha_rush_mor_in=0.898;    % inflow rate  from community to the hub in morning rush hours
alpha_off_mor_in=0.133;     % inflow rate  from community to the hub in morning off-peak hours
alpha_rush_eve_in=0.198;    % inflow rate  from community to the hub in evening rush hours
alpha_off_eve_in=0.128;     % inflow rate  from community to the hub in evening off-peak hours
gamma_rush_mor_in=0.555;    %  outflow rate  from hub to the community in morning rush hours
gamma_off_mor_in=0.0904;     % outflow rate  from hub to the community in morning off-peak hours
gamma_rush_eve_in=0.0755;    %  outflow rate  from hub to the community in evening rush hours
gamma_off_eve_in=0.307;     % outflow rate  from hub to the community in evening off-peak hours


N_in=100000;           % initial value for number of the total population of the community
I_c_in=2;              % initial value for I_c(0)


param_estimate=[alpha_rush_mor_in, alpha_off_mor_in, alpha_rush_eve_in, alpha_off_eve_in,...
gamma_rush_mor_in, gamma_off_mor_in, gamma_rush_eve_in, gamma_off_eve_in, N_in,I_c_in ]; 
%-----------------Estimating Function ----------------------------------
LB=[zeros(1,8),  100000, 0 ]; %lower Bound
UB=[ones(1,4),1,.9,.9,1, 1000000, 10]; % Upper Bound

options = optimset('OutputFcn', @myoutput, 'TolFun', 1e-6);
[pa_estimate,fval]=fminsearchbnd(@(param)fmin_estimatettc_one(param), param_estimate, ...  % Estimating function
    LB,UB,options)
[error_estimate  x_rush_mor x_off_mor x_rush_eve x_off_eve ]=fmin_estimatettc_one(pa_estimate);
%---------------Plotting Process------------------------

t_rush_mor=linspace(0,3,4);  % time interval for morning rush hours
t_off_mor=linspace(3,9,7);   % time interval for morning off-peak hours
t_rush_eve=linspace(9,13,5); % time interval for evening rush hours
t_off_eve=linspace(13,20,8); % time interval for evening off-peak hours
t_total=linspace(0,20,21);   % total time interval

 % initialize and evaluate the error estimation matrices to zero

 error_estimate_rush_mor=zeros(length(t_rush_mor),1); 
 error_estimate_off_mor=zeros(length(t_off_mor),1);
 error_estimate_rush_eve=zeros(length(t_rush_eve),1);
 error_estimate_off_eve=zeros(length(t_off_eve),1);


 % plot the figure of the number of the exposed individuals over a single day in the
 % community and hub based on the estimated parameters

figure(1)
plot (t_rush_mor, x_rush_mor(:,3),'-r',t_rush_mor,x_rush_mor(:,7),'-b' ,t_off_mor, x_off_mor(:,3),'-r', t_off_mor, x_off_mor(:,7),'-b' , ...
    t_rush_eve, x_rush_eve(:,3),'-r' , t_rush_eve, x_rush_eve(:,7), '-b' ,t_off_eve,x_off_eve(:,3),'-r' ,t_off_eve, x_off_eve(:,7),'-b' , 'LineWidth',2)

legend('community','Hub')
xlabel('Time (Hour)')
ylabel(' Number of Infected Individuals')

% ----------- Dataset-----------

function stop = myoutput(pa_estimate,optimvalues,state);
        stop = false;
        if isequal(state,'iter')
          history = [history; pa_estimate];
        end

end
 %% Confidence Interval  Process
 history(:,10)
 N_size=size(history(:,3),1)
 N_2=size(history(:,2),1);
 N_Sample_1=size(history(1360:N_size,1),1);
 N_Sample_2=size(history(1390:N_size,2),1);
 N_Sample_3=size(history(1405:N_size,3),1);
 N_Sample_4=size(history(1400:N_size,4),1);
 N_Sample_5=size(history(1405:N_size,5),1);
 N_Sample_6=size(history(1405:N_size,6),1);
 N_Sample_7=size(history(1400:N_size,7),1);
 N_Sample_8=size(history(1400:N_size,8),1);
 N_Sample_9=size(history(1405:N_size,9),1);
 N_Sample_10=size(history(1405:N_size,10),1);
%% Mean
 
 mu_alpha_rush_mor=mean(history(1360:N_size,1));
 mu_alpha_off_mor=mean(history(1390:N_size,2));
 mu_alpha_rush_eve=mean(history(1405:N_size,3));
 mu_alpha_off_eve=mean(history(1400:N_size,4));
 mu_gamma_rush_mor=mean(history(1405:N_size,5));
 mu_gamma_off_mor=mean(history(1405:N_size,6));
 mu_gamma_rush_eve=mean(history(1400:N_size,7));
 mu_gamma_off_eve=mean(history(1400:N_size,8));
 mu_N_com=mean(history(1405:N_size,9));
 mu_I=mean(history(1405:N_size,10));

%% SD
 SD_alpha_rush_mor=std(history(1360:N_size,1));
 SD_alpha_off_mor=std(history(1390:N_size,2));
 SD_alpha_rush_eve=std(history(1405:N_size,3));
 SD_alpha_off_eve=std(history(1400:N_size,4));
 SD_gamma_rush_mor=std(history(1405:N_size,5));
 SD_gamma_off_mor=std(history(1405:N_size,6));
 SD_gamma_rush_eve=std(history(1400:N_size,7));
 SD_gamma_off_eve=std(history(1400:N_size,8));
 SD_N_com=std(history(1405:N_size,9));
 SD_I=std(history(1405:N_size,10));



%% CI

CI_alpha_rush_mor_left=mu_alpha_rush_mor-1.96*(SD_alpha_rush_mor/sqrt(N_Sample_1));
CI_alpha_rush_mor_right=mu_alpha_rush_mor+1.96*(SD_alpha_rush_mor/sqrt(N_Sample_1));

CI_alpha_off_mor_left=mu_alpha_off_mor-1.96*(SD_alpha_off_mor/sqrt(N_Sample_2));
CI_alpha_off_mor_right=mu_alpha_off_mor+1.96*(SD_alpha_off_mor/sqrt(N_Sample_2));

CI_alpha_rush_eve_left=mu_alpha_rush_eve-1.96*(SD_alpha_rush_eve/sqrt(N_Sample_3));
CI_alpha_rush_eve_right=mu_alpha_rush_eve+1.96*(SD_alpha_rush_eve/sqrt(N_Sample_3));

CI_alpha_off_eve_left=mu_alpha_off_eve-1.96*(SD_alpha_off_eve/sqrt(N_Sample_4));
CI_alpha_off_eve_right=mu_alpha_off_eve+1.96*(SD_alpha_off_eve/sqrt(N_Sample_4));

CI_gamma_rush_mor_left=mu_gamma_rush_mor-1.96*(SD_gamma_rush_mor/sqrt(N_Sample_5));
CI_gamma_rush_mor_right=mu_gamma_rush_mor+1.96*(SD_gamma_rush_mor/sqrt(N_Sample_5));


CI_gamma_off_mor_left=mu_gamma_off_mor-1.96*(SD_gamma_off_mor/sqrt(N_Sample_6));
CI_gamma_off_mor_right=mu_gamma_off_mor+1.96*(SD_gamma_off_mor/sqrt(N_Sample_6));

CI_gamma_rush_eve_left=mu_gamma_rush_eve-1.96*(SD_gamma_rush_eve/sqrt(N_Sample_7));
CI_gamma_rush_eve_right=mu_gamma_rush_eve+1.96*(SD_gamma_rush_eve/sqrt(N_Sample_7));

CI_gamma_off_eve_left=mu_gamma_off_eve-1.96*(SD_gamma_off_eve/sqrt(N_Sample_8));
CI_gamma_off_eve_right=mu_gamma_off_eve+1.96*(SD_gamma_off_eve/sqrt(N_Sample_8));


CI_N_com_left=mu_N_com-1.96*(SD_N_com/sqrt(N_Sample_9));
CI_N_com_right=mu_N_com+1.96*(SD_N_com/sqrt(N_Sample_9));


CI_I_left=mu_I-1.96*(SD_I/sqrt(N_Sample_10))
CI_I_right=mu_I+1.96*(SD_I/sqrt(N_Sample_10))

%----------- Model Function----------------

    function [error_estimate  x_rush_mor x_off_mor x_rush_eve x_off_eve]=fmin_estimatettc_one(param)

alpha_rush_mor=param(1);    % inflow rate  from community to the hub in morning rush hours
alpha_off_mor=param(2);     % inflow rate  from community to the hub in morning off-peak hours
alpha_rush_eve=param(3);    % inflow rate  from community to the hub in evening rush hours
alpha_off_eve=param(4);     % inflow rate  from community to the hub in evening off-peak hours
gamma_rush_mor=param(5);    %  outflow rate  from hub to the community in morning rush hours
gamma_off_mor=param(6);     % outflow rate  from hub to the community in morning off-peak hours
gamma_rush_eve=param(7);    %  outflow rate  from hub to the community in evening rush hours
gamma_off_eve=param(8);     %  outflow rate  from hub to the community in evening off-peak hours
 N=param(9);               % Number of the total population of the community
 I_0c=param(10);           % Initial data of infectees in the community

%-------------------- Initial Data-----

 N_h_Finch=61323;    % Number of the indviduals visiting Finch station daily

E_0c=0;              % initial data for number of the exposed individuals in the community
R_0c=0;              % initial data for number of the recovered individuals in the community
S_0c=N-(I_0c+E_0c);   %initial data for number of the scuceptible individuals in the community
S_0h=389;              %initial data for number of the scuceptible individuals in the hub
E_0h=0;              % initial data for number of the exposed individuals in the hub
I_0h=0;              % initial data for number of the infected individuals in the hub
R_0h=0;               % initial data for number of the recovered individuals in the hub
IC_estimate=[S_0c,E_0c,I_0c,R_0c,S_0h,E_0h,I_0h,R_0h]; % Vector of the initial data

%---------------- Data Fitting -------

%% Rush-peak hour Morning: 6AM-9AM
t_rush_mor=linspace(0,3,4);
IC_estimate_rush_mor=IC_estimate;
p_rush_mor=[alpha_rush_mor,gamma_rush_mor];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_mor]=ode45(@(t,x_rush_mor)ttccase_mob(t,x_rush_mor,p_rush_mor),t_rush_mor,IC_estimate_rush_mor,op);


%% Off-peak hour midday: 9AM-3PM
t_off_mor=linspace(3,9,7);

IC_estimate_off_mor=[x_rush_mor(4,1),x_rush_mor(4,2),x_rush_mor(4,3),x_rush_mor(4,4),x_rush_mor(4,5),x_rush_mor(4,6),x_rush_mor(4,7),x_rush_mor(4,8)]; % Vector of the initial data in off-peak hours
p_off_mor=[alpha_off_mor,gamma_off_mor];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_mor]=ode45(@(t,x_off_mor)ttccase_mob(t,x_off_mor,p_off_mor),t_off_mor,IC_estimate_off_mor,op);

%% Rush hour evening: 3AM-7AM
t_rush_eve=linspace(9,13,5);

IC_estimate_rush_eve=[x_off_mor(7,1),x_off_mor(7,2),x_off_mor(7,3),x_off_mor(7,4),x_off_mor(7,5),x_off_mor(7,6),x_off_mor(7,7),x_off_mor(7,8)]% Vector of the initial data in rush hours
p_rush_eve=[alpha_rush_eve,gamma_rush_eve];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_rush_eve]=ode45(@(t,x_rush_eve)ttccase_mob(t,x_rush_eve,p_rush_eve),t_rush_eve,IC_estimate_rush_eve,op);

%% Off-peak hour evening: 7AM-2AM
t_off_eve=linspace(13,20,8);
IC_estimate_off_eve=[x_rush_eve(5,1),x_rush_eve(5,2),x_rush_eve(5,3),x_rush_eve(5,4),x_rush_eve(5,5),x_rush_eve(5,6),x_rush_eve(5,7),x_rush_eve(5,8)]; % Vector of the initial data in off-peak hours
p_off_eve=[alpha_off_eve,gamma_off_eve];
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[~,x_off_eve]=ode45(@(t,x_off_eve)ttccase_mob(t,x_off_eve,p_off_eve),t_off_eve,IC_estimate_off_eve,op);
%% Result

 
%% executing error estimation matrices 
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
%% total error estimate matrix
EE=[error_estimate_rush_mor'; error_estimate_off_mor';error_estimate_rush_eve';error_estimate_off_eve'];

% cumulative of EE matrix
EE_cum=cumsum(EE);



t_total=linspace(0,20,21);

x_rush_cum_total_mor_h=cumsum(x_rush_total_mor_h);
x_rush_cum_total_mor_c=cumsum(x_rush_total_mor_c);

x_off_cum_total_mor_h=cumsum(x_off_total_mor_h);
x_off_cum_total_mor_c=cumsum(x_off_total_mor_c);

x_rush_cum_total_eve_h=cumsum(x_rush_total_eve_h);
x_rush_cum_total_eve_c=cumsum(x_rush_total_eve_c);

x_off_cum_total_eve_h=cumsum(x_off_total_eve_h);
x_off_cum_total_eve_c=cumsum(x_off_total_eve_c);

% total number of individuals visiting hub over a single day
N_h_total=x_rush_cum_total_mor_h(end)+x_off_cum_total_mor_h(end)-x_off_cum_total_mor_h(1)+...
   x_rush_cum_total_eve_h(end)- x_rush_cum_total_eve_h(1)+x_off_cum_total_eve_h(end)-x_off_cum_total_eve_h(1)
error_estimate= EE_cum(end)+(N_h_total-N_h_Finch)^2;

%% plot the fitted curve
figure(3)
 plot(t_rush_mor,x_rush_total_mor_h,'-b', t_off_mor,x_off_total_mor_h,'-b', ...
     t_rush_eve,x_rush_total_eve_h,'-b',t_off_eve,x_off_total_eve_h,'-b', ...
     t_total, Data,'k .', 'markersize', 12,'LineWidth',2)
 xlabel('Time (Hour)')
ylabel(' Number of Individuals Visiting Finch Station in a Single Day')


    end

% the follwing outputs are used in CTMC method

Inf_cum_ode_rush_mor=[x_rush_mor(:,3)]
Inf_cum_ode_off_mor=[ x_off_mor(:,3)]
Inf_cum_ode_rush_eve=[x_rush_eve(:,3)]
Inf_cum_ode_off_eve=[ x_off_eve(:,3)]
Inf_hub_ode_rush_mor=[x_rush_mor(:,7)]
Inf_hub_ode_off_mor=[ x_off_mor(:,7)]
Inf_hub_ode_rush_eve=[x_rush_eve(:,7)]
Inf_hub_ode_off_eve=[ x_off_eve(:,7)]
x_rush_total_mor_h=[x_rush_total_mor_h(:)]
x_off_total_mor_h=[x_off_total_mor_h(:)]
x_rush_total_eve_h=[x_rush_total_eve_h(:)]
x_off_total_eve_h=[x_off_total_eve_h(:)]



end


 