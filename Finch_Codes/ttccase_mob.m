function dIdt=ttccase_mob(t, I, par)

alpha=par(1);
gamma=par(2);


    dIdt=[-alpha*I(1)+gamma*I(5); %I(1)-S_c
       -alpha*I(2)+gamma*I(6); %I(2)-L_c
        -alpha*I(3)+gamma*I(7);%I(3)-I_c
        -alpha*I(4)+gamma*I(8); %I(4)-R_c
        -gamma*I(5)+alpha*I(1); %I(5)-S_h
        -gamma*I(6)+alpha*I(2); %I(6)-L_h
        -gamma*I(7)+alpha*I(3);%I(7)-I_h
        -gamma*I(8)+alpha*I(4); %I(8)-R_h
    ];