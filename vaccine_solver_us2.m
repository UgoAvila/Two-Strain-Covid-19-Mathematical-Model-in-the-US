function [t,S,E,I,A,V,EU,IU,AU,D,R]=vaccine_solver_us2(v,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0)
% This function computes the solutions of the model with vaccination rate
% v, initial values (S0,E0,...,R0), and total population size N,
% over the time range t


%% Solving the model
options=odeset('NonNegative',(1:10));
[t,y] = ode45(@(t,y) vaccine_odes_us2(t,y,v),t,[S0;E0;I0;A0;V0;EU0;IU0;AU0;D0;R0],options);


S=y(:,1);
E=y(:,2);
I=y(:,3);
A=y(:,4);
V=y(:,5);
EU=y(:,6);
IU=y(:,7);
AU=y(:,8);
D=y(:,9);
R=y(:,10);
end