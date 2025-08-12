function dIdt=ttccase_one(t, I, par)
beta_c=par(1);
beta_h=par(2);
alpha=par(3);
gamma=par(4);
k=par(5);         %saturation number for the force of infection of the community
k_1=par(6);
p=par(7);
epsilon_c=par(8);
epsilon_h=par(9);
delta_c=par(10);
delta_h=par(11);
%saturation number for the force of infection of the hub
    F_1=I(3)+I(7);
    F=F_1/(1+k*F_1);  %Force of infection for the community
   G=F_1/(1+k_1*F_1); %Force of infection for the hub
   
    dIdt=[-beta_c*F*I(1)-alpha*I(1)+(1-p)*gamma*I(5); %I(1)-S_c
        beta_c*F*I(1)+p*gamma*I(5)-alpha*I(2)+gamma*I(6)-epsilon_c*I(2); %I(2)-L_c
        epsilon_c*I(2)-alpha*I(3)+gamma*I(7)-delta_c*I(3);%I(3)-I_c
         delta_c*I(3)-alpha*I(4)+gamma*I(8); %I(4)-R_c
        -beta_h*G*I(5)-gamma*I(5)+(1-p)*alpha*I(1); %I(5)-S_h
        beta_h*G*I(5)-gamma*I(6)+alpha*I(2)+p*alpha*I(1)-epsilon_h*I(6); %I(6)-L_h
        epsilon_h*I(6)-delta_h*I(7)-gamma*I(7)+alpha*I(3);%I(7)-I_h
        delta_h*I(7)-gamma*I(8)+alpha*I(4); %I(8)-R_h
    ];