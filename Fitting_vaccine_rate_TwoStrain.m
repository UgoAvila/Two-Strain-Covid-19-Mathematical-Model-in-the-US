% This code plots the reported number of vaccinated people (with one or two
% doses) and compares it with the model simulations

load('simul-before-vaccineTwoStrain3.mat')    % This is used to set the initial conditions

% Reading the reported number of vaccinations from the file VaccineUsa.csv
% only for Pfizer
vaccUSA=readtable('United States.csv','ReadRowNames',true);
vaccDates=datetime(vaccUSA{1:end,1},'Format','dd/MM/yy');

V_exp_us = vaccUSA{1:end,5}; % People who have received a dose from Pfizer Vaccine
V2_exp = vaccUSA{1:end,6}; % People who have received both doses
V1_exp = V_exp_us - V2_exp; % People who have received only the first dose




%% Seting initial conditions based on 'simul-before-vaccine.mat'

S0 = Sbv(end);
E0 = Ebv(end);
I0 = Ibv(end);
A0 = Abv(end);
V0 = Vbv(end);
EU0 = 456019;
IU0 = 456019;
AU0 = ((1-0.2)/0.2)*IU0;
R0 = Rbv(end);
RU0=RUbv(end);
D0 = Dbv(end);

newinitdate=357; 

numberofdays = 96;            %Number of days for the predictions
t=linspace(0,numberofdays,numberofdays+1);
newdata_long=dateshift(data_long(43),'start','day',0:t(end));

%% Plotting number of vaccinated people

epsilonl=1.; epsilona=0.; epsilonla=1.; epsilonlb=1.; alpha=0.; 
% We set efficacy rates to 1, so V1 and V2 represent the cumulative number
% of vaccinations

vacc = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*(X>=81);     % Vaccination rate (manually fitted to data)

[t,S,E,I,A,V,EU,IU,AU,D,R,RU]=vaccine_solver_us(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,S,E,I,A,V,EU,IU,AU,D,R,RU]
tStart=datetime('2021-01-14'); tEnd=datetime('2021-04-18');
% Initial and final dates for the graph

figure
subplot(1,2,1)
maxpopul=12e7;   %Max. population for the graph
perc=maxpopul*100/N;    %Percent of Mexican population
yyaxis left
title('Vaccinated')
plot(newdata_long,V,'-k','LineWidth',2)
hold on
plot(vaccDates,V2_exp,'ob')
%ylim([0 maxpopul])
ylabel('people')
legend({'Model solutions','Reported data'},'Location','northwest')
grid on
set(gca,'YColor','k')
yyaxis right
ylim([0 perc])
ylabel('% of USA population')
ytickformat('percentage')
set(gca,'YColor','k')
xlim([tStart tEnd])

%% Proyections of fully vaccinated with Pfizer

newinitdate=357;   % 12 Jan 2021

numberofdays = 735;            %Number of days for the predictions
t=linspace(0,numberofdays,numberofdays+1);
newdata_long=dateshift(data_long(43),'start','day',0:t(end));

%% Plotting number of vaccinated people

%%
epsilonl=1.; epsilona=0.; epsilonla=1.; epsilonlb=1.; alpha=0.;
vacc = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*(X>=81);    % Vaccination rate (manually fitted to data)
    % Vaccination rate (manually fitted to data)

[t,S,E,I,A,V,EU,IU,AU,D,R,RU]=vaccine_solver_us(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,S,E,I,A,V,EU,IU,AU,D,R,RU]


%% Plotting number of vaccinated people for 200%

epsilonl=1.; epsilona=0.; epsilonla=1.; epsilonlb=1.; alpha=0.;
vaccb = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*( (X>=81)&(X<109) ) ...
          + (0.070842).*(X>=109);    % 200% vaccination rate starting on 1st May

[t,S,E,I,A,V,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
%%
epsilonl=1.; epsilona=0.; epsilonla=1.; epsilonlb=1.; alpha=0.;
vaccd = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*( (X>=81)&(X<109) ) ...
          + (0.01521).*(X>=109);      % 50% vaccination rate starting on 1st May

[t,Sc,Ec,Ic,Ac,Vc,EUc,IUc,AUc,Dc,Rc,RUc]=vaccine_solver_us(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
%%
subplot(1,2,2)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-12-25');
maxpopul=3e8;   %Max. population for the graph
perc=maxpopul*100/N;    %Percent of Mexican population
yyaxis left
title('Vaccinated')
hold on
plot(newdata_long,V,'-k','LineWidth',2)
plot(newdata_long,Vb,'--b','LineWidth',2)
plot(newdata_long,Vc,'--r','LineWidth',2)
ylim([0 maxpopul])
ylabel('millions of people')
legend({'Baseline vaccination rate','200% vaccination rate', '50% vaccination rate'},'Location','northwest')
grid on
set(gca,'YColor','k')
yyaxis right
ylim([0 perc])
ylabel('% of US population')
ytickformat('percentage')
set(gca,'YColor','k')
xlim([tStart tEnd])



