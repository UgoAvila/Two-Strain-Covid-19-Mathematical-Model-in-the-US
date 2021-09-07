load('simul-before-vaccineTwoStrain3.mat')
%% Seting initial conditions for simulations:

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
newinitdate=357;   % 12 Ene 2021


numberofdays = 735;
t=linspace(0,numberofdays,numberofdays+1);
newdata_long=dateshift(data_long(43),'start','day',0:t(end));

%% Defining vaccination rates:
novac = @(X) 0.0;   % No vaccination
%%
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vacc = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*(X>=81);    % Vaccination rate (manually fitted to data)

%%
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vaccb = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*( (X>=81)&(X<109) ) ...
          + (0.070842).*(X>=109);     % 200% vaccination rate starting on 1st May%%
%%
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vaccd =@(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*( (X>=81)&(X<109) ) ...
          + (0.01521).*(X>=109);     % 50% vaccination rate starting on 1st May
      

%% LOW TRANSMISSION RATE
beta1 = 0.2*0.8;
beta2 = 0.033*0.8;
beta3=0.3*0.8;
beta4=0.0495*0.8;
[t,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(novac,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,Sc,Ec,Ic,Ac,Vc,EUc,IUc,AUc,Dc,Rc,RUc]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%%Figure 1 
figure
subplot(3,2,1)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-08-30');
plot(newdata_long,Io+IUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Ia+IUa,'-b','LineWidth',2)
plot(newdata_long,Ib+IUb,'-.k','LineWidth',2)
plot(newdata_long,Ic+IUc, '-.g', 'LineWidth',2)
title('Low Transmission Rate')
ylabel('I(t)+IU(t)')
xlim([tStart tEnd])
%ylim([0 1e6])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate', '50% vaccination rate'},'Location','northwest')
grid on

subplot(3,2,2)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-07-30');
plot(newdata_long,Ao+AUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Aa+AUa,'-b','LineWidth',2)
plot(newdata_long,Ab+AUb,'-.k','LineWidth',2)
plot(newdata_long,Ac+AUc, '-.g', 'LineWidth',2)
title('Low Transmission Rate')
ylabel('A(t)+AU(t)')
xlim([tStart tEnd])
%ylim([0 1e6])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate', '50% vaccination rate'},'Location','northwest')
grid on

%% Baseline transmission rate. 
beta1 = 0.2;
beta2 = 0.033;
beta3=0.3;
beta4=0.0495;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(novac,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,Sc,Ec,Ic,Ac,Vc,EUc,IUc,AUc,Dc,Rc,RUc]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

subplot(3,2,3)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-07-30');
plot(newdata_long,Io+IUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Ia+IUa,'-b','LineWidth',2)
plot(newdata_long,Ib+IUb,'-.k','LineWidth',2)
plot(newdata_long,Ic+IUc,'-.g','LineWidth',2)
title('Normal transmission rate')
ylabel('I(t)+IU(t)')
xlim([tStart tEnd])
%ylim([0 1e6])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate','50% vaccination rate' },'Location','northwest')
grid on

subplot(3,2,4)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-07-30');
plot(newdata_long,Ao+AUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Aa+AUa,'-b','LineWidth',2)
plot(newdata_long,Ab+AUb,'-.k','LineWidth',2)
plot(newdata_long,Ac+AUc,'-.g','LineWidth',2)
title('Normal Transmission Rate')
ylabel('A(t)+AU(t)')
xlim([tStart tEnd])
%ylim([0 1e6])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate','50% vaccination rate'},'Location','northwest')
grid on

%% High Transmission Rate
beta1 = 0.2*1.5;
beta2 = 0.033*1.5;
beta3=0.3*1.5;
beta4=0.0495*1.5;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(novac,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,Sc,Ec,Ic,Ac,Vc,EUc,IUc,AUc,Dc,Rc,RUc]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);


subplot(3,2,5)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-07-30');
plot(newdata_long,Io+IUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Ia+IUa,'-b','LineWidth',2)
plot(newdata_long,Ib+IUb,'-.k','LineWidth',2)
plot(newdata_long,Ic+IUc,'-.g','LineWidth',2)
title('High transmission rate')
ylabel('I(t)+IU(t)')
xlim([tStart tEnd])
%ylim([0 1e6])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate','50% vaccination rate' },'Location','northwest')
grid on

subplot(3,2,6)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-07-30');
plot(newdata_long,Ao+AUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Aa+AUa,'-b','LineWidth',2)
plot(newdata_long,Ab+AUb,'-.k','LineWidth',2)
plot(newdata_long,Ac+AUc,'-.g','LineWidth',2)
title('High transmission rate')
ylabel('A(t)+AU(t)')
xlim([tStart tEnd])
%ylim([0 1e6])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate','50% vaccination rate' },'Location','northwest')
grid on


%%Code for obtaining figure 2 of the supplementary material 
%% LOW TRANSMISSION RATE
beta1 = 0.2*0.8;
beta2 = 0.033*0.8;
beta3=0.3*0.8;
beta4=0.0495*0.8;
[t,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(novac,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,Sc,Ec,Ic,Ac,Vc,EUc,IUc,AUc,Dc,Rc,RUc]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

figure
subplot(3,2,1)
plot(newdata_long,Ro,'--r','LineWidth',2)
hold on
plot(newdata_long,Ra,'-b','LineWidth',2)
plot(newdata_long,Rb,'-.k','LineWidth',2)
plot(newdata_long,Rc,'-.g', 'LineWidth',2)
title('Low Transmission Rate')
ylabel('R(t)')
xlim([tStart tEnd])
%ylim([1e5 9e5])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate', '50% vaccination rate'},'Location','southeast')
grid on

subplot(3,2,2)
plot(newdata_long,Do,'--r','LineWidth',2)
hold on
plot(newdata_long,Da,'-b','LineWidth',2)
plot(newdata_long,Db,'-.k','LineWidth',2)
plot(newdata_long,Dc,'-.g', 'LineWidth',2)
title('Low Transmission Rate')
ylabel('D(t)')
xlim([tStart tEnd])
%ylim([1e5 9e5])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate', '50% vaccination rate'},'Location','southeast')
grid on
%% Baseline transmission rate. 
beta1 = 0.2;
beta2 = 0.033;
beta3=0.3;
beta4=0.0495;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(novac,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,Sc,Ec,Ic,Ac,Vc,EUc,IUc,AUc,Dc,Rc,RUc]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

subplot(3,2,3)
plot(newdata_long,Ro,'--r','LineWidth',2)
hold on
plot(newdata_long,Ra,'-b','LineWidth',2)
plot(newdata_long,Rb,'-.k','LineWidth',2)
plot(newdata_long,Rc,'-.g','LineWidth',2)
ylabel('R(t)')
title('Normal Transmission Rate')
xlim([tStart tEnd])
%ylim([1e5 9e5])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate','50% vaccination rate'},'Location','northwest')
grid on

subplot(3,2,4)
plot(newdata_long,Do,'--r','LineWidth',2)
hold on
plot(newdata_long,Da,'-b','LineWidth',2)
plot(newdata_long,Db,'-.k','LineWidth',2)
plot(newdata_long,Dc,'-.g', 'LineWidth',2)
title('Normal Transmission Rate')
ylabel('D(t)')
xlim([tStart tEnd])
%ylim([1e5 9e5])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate', '50% vaccination rate'},'Location','southeast')
grid on

%% High Transmission Rate
beta1 = 0.2*1.5;
beta2 = 0.033*1.5;
beta3=0.3*1.5;
beta4=0.0495*1.5;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(novac,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);
[t,Sc,Ec,Ic,Ac,Vc,EUc,IUc,AUc,Dc,Rc,RUc]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

subplot(3,2,5)
plot(newdata_long,Ro,'--r','LineWidth',2)
hold on
plot(newdata_long,Ra,'-b','LineWidth',2)
plot(newdata_long,Rb,'-.k','LineWidth',2)
plot(newdata_long,Rc,'-.g','LineWidth',2)
title('High transmission rate')
ylabel('R(t)')
xlim([tStart tEnd])
%ylim([1e5 9e5])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate','50% vaccination rate'},'Location','northwest')
grid on

subplot(3,2,6)
plot(newdata_long,Do,'--r','LineWidth',2)
hold on
plot(newdata_long,Da,'-b','LineWidth',2)
plot(newdata_long,Db,'-.k','LineWidth',2)
plot(newdata_long,Dc,'-.g', 'LineWidth',2)
title('High Transmission Rate')
ylabel('D(t)')
xlim([tStart tEnd])
%ylim([1e5 9e5])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate', '50% vaccination rate'},'Location','southeast')
grid on


% Plotting the solutions of the model using different values for the
% efficacy rates eta1 and eta2

load('simul-before-vaccineTwoStrain3.mat')
%% Seting initial conditions for simulations:

S0 = Sbv(end);
E0 = Ebv(end);
I0 = Ibv(end);
A0 = Abv(end);
V0 = Vbv(end);
EU0 = 44000;
IU0 = 44000;
AU0 = ((1-0.2)/0.2)*IU0;
R0 = Rbv(end);
RU0=RUbv(end);
D0 = Dbv(end);
newinitdate=357;   % 12 Ene 2021

numberofdays = 735;
t=linspace(0,numberofdays,numberofdays+1);
newdata_long=dateshift(data_long(43),'start','day',0:t(end));
%% Low Vaccination Rate
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vaccd = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*( (X>=81)&(X<109) ) ...
          + (0.01521).*(X>=109);      % 50% vaccination rate starting on 1st May
      
%% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% GRAPHS
%%Code for obtaininf figure S3 of the supplementary material
figure
subplot(3,2,1)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-07-30');
plot(newdata_long,Io+IUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Ia+IUa,'-b','LineWidth',2)
plot(newdata_long,Ib+IUb,'-.k','LineWidth',2)
title('Symptomatic Infected Individuals')
ylabel('I(t)+IU(t)')
xlim([tStart tEnd])
legend({'\epsilon_a=0.0067, \epsilon_La=0.72, \epsilon_Lb=0.87, \epsilon_L=0.89','\epsilon_a=3.545e-3, \epsilon_La=0.75, \epsilon_Lb=0.9, \epsilon_L=0.913', '\epsilon_a=3.77e-4, \epsilon_La=0.95, \epsilon_Lb=0.95, \epsilon_L=0.932'},'Location','northwest')
grid on

subplot(3,2,2)
plot(newdata_long,Ao+AUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Aa+AUa,'-b','LineWidth',2)
plot(newdata_long,Ab+AUb,'-.k','LineWidth',2)
title('Asymptomatic Infected Individuals')
ylabel('A(t)+AU(t)')
xlim([tStart tEnd])
legend({'\epsilon_a=0.0067, \epsilon_La=0.72, \epsilon_Lb=0.87, \epsilon_L=0.89','\epsilon_a=3.545e-3, \epsilon_La=0.75, \epsilon_Lb=0.9, \epsilon_L=0.913', '\epsilon_a=3.77e-4, \epsilon_La=0.95, \epsilon_Lb=0.95, \epsilon_L=0.932'},'Location','southwest')
grid on



%% (Normal Vaccination Rate)
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vacc = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*(X>=81);     % Vaccination rate (manually fitted to data)



%% Baseline efficacy (Normal Vaccination)
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%%Code for obtaining figure S3 of the supplementary material

subplot(3,2,3)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-07-30');
plot(newdata_long,Io+IUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Ia+IUa,'-b','LineWidth',2)
plot(newdata_long,Ib+IUb,'-.k','LineWidth',2)
title('Symptomatic Infected Individuals')
ylabel('I(t)+ IU(t)')
xlim([tStart tEnd])
legend({'\epsilon_a=0.0067, \epsilon_La=0.72, \epsilon_Lb=0.87, \epsilon_L=0.89','\epsilon_a=3.545e-3, \epsilon_La=0.75, \epsilon_Lb=0.9, \epsilon_L=0.913', '\epsilon_a=3.77e-4, \epsilon_La=0.95, \epsilon_Lb=0.95, \epsilon_L=0.932'},'Location','southwest')
grid on

subplot(3,2,4)

plot(newdata_long,Ao+AUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Aa+AUa,'-b','LineWidth',2)
plot(newdata_long,Ab+AUb,'-.k','LineWidth',2)
title('Asymptomatic Infected Individuals')
ylabel('A(t)+AU(t)')
xlim([tStart tEnd])
legend({'\epsilon_a=0.0067, \epsilon_La=0.72, \epsilon_Lb=0.87, \epsilon_L=0.89','\epsilon_a=3.545e-3, \epsilon_La=0.75, \epsilon_Lb=0.9, \epsilon_L=0.913', '\epsilon_a=3.77e-4, \epsilon_La=0.95, \epsilon_Lb=0.95, \epsilon_L=0.932'},'Location','southwest')
grid on

%%High Vaccination Rate
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vaccb = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*( (X>=81)&(X<109) ) ...
          + (0.070842).*(X>=109); % 200 vaccination rate 
      
 %% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% GRAPHS

subplot(3,2,5)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-07-30');
plot(newdata_long,Io+IUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Ia+IUa,'-b','LineWidth',2)
plot(newdata_long,Ib+IUb,'-.k','LineWidth',2)
title('Symptomatic Infected Individuals')
ylabel('I(t)+IU(t)')
xlim([tStart tEnd])
legend({'\epsilon_a=0.0067, \epsilon_La=0.72, \epsilon_Lb=0.87, \epsilon_L=0.89','\epsilon_a=3.545e-3, \epsilon_La=0.75, \epsilon_Lb=0.9, \epsilon_L=0.913', '\epsilon_a=3.77e-4, \epsilon_La=0.95, \epsilon_Lb=0.95, \epsilon_L=0.932'},'Location','southwest')
grid on

subplot(3,2,6)
plot(newdata_long,Ao+AUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Aa+AUa,'-b','LineWidth',2)
plot(newdata_long,Ab+AUb,'-.k','LineWidth',2)
title('Asymptomatic Infected Individuals')
ylabel('A(t)+AU(t)')
xlim([tStart tEnd])
legend({'\epsilon_a=0.0067, \epsilon_La=0.72, \epsilon_Lb=0.87, \epsilon_L=0.89','\epsilon_a=3.545e-3, \epsilon_La=0.75, \epsilon_Lb=0.9, \epsilon_L=0.913', '\epsilon_a=3.77e-4, \epsilon_La=0.95, \epsilon_Lb=0.95, \epsilon_L=0.932'},'Location','southwest')
grid on

%% Low Vaccination Rate
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vaccd = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*( (X>=81)&(X<109) ) ...
          + (0.01521).*(X>=109);      % 50% vaccination rate starting on 1st May
      
%% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% GRAPHS
subplot(3,3,1)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Io,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Low Efficiency (LVR)')
plot(newdata_long,IUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,4)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ia,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Baseline Efficiency (LVR)')
plot(newdata_long,IUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,7)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ib,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('High Efficiency (LVR) ')
plot(newdata_long,IUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%% (Normal Vaccination Rate)
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vacc = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*(X>=81);     % Vaccination rate (manually fitted to data)



%% Baseline efficacy (Normal Vaccination)
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% GRAPHS
subplot(3,3,2)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Io,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Low Efficiency (NVR)')
plot(newdata_long,IUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,5)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ia,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Baseline Efficiency (NVR)')
plot(newdata_long,IUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,8)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ib,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('High Efficiency (NVR)')
plot(newdata_long,IUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';


%% Defining vaccination rates:
%%High Vaccination Rate
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vaccb = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*( (X>=81)&(X<109) ) ...
          + (0.070842).*(X>=109); % 200 vaccination rate 
      
 %% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% GRAPHS
subplot(3,3,3)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Io,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Low Efficiency (HVR)')
plot(newdata_long,IUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,6)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ia,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Baseline Efficiency (HVR)')
plot(newdata_long,IUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,9)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ib,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('High Efficiency (HVR)')
plot(newdata_long,IUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';


%%Code for obtaining figure S4 
%% Low Vaccination Rate
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vaccd = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*( (X>=81)&(X<109) ) ...
          + (0.01521).*(X>=109);      % 50% vaccination rate starting on 1st May
      
%% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% GRAPHS
subplot(3,3,1)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ao,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Low Efficiency (LVR)')
plot(newdata_long,AUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,4)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Aa,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Baseline Efficiency (LVR)')
plot(newdata_long,AUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,7)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ab,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('High Efficiency (LVR) ')
plot(newdata_long,AUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%% (Normal Vaccination Rate)
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vacc = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*(X>=81);     % Vaccination rate (manually fitted to data)



%% Baseline efficacy (Normal Vaccination)
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% GRAPHS
subplot(3,3,2)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ao,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Low Efficiency (NVR)')
plot(newdata_long,IUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,5)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Aa,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Baseline Efficiency (NVR)')
plot(newdata_long,AUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,8)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ab,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('High Efficiency (NVR)')
plot(newdata_long,AUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%% Defining vaccination rates:
%%High Vaccination Rate
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vaccb = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*( (X>=81)&(X<109) ) ...
          + (0.070842).*(X>=109); % 200 vaccination rate 
      
 %% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% GRAPHS
subplot(3,3,3)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ao,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Low Efficiency (HVR)')
plot(newdata_long,AUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,6)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Aa,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Baseline Efficiency (HVR)')
plot(newdata_long,AUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,9)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ab,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('High Efficiency (HVR)')
plot(newdata_long,AUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%%Code for obtaining figure S5 
%% Low Vaccination Rate
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vaccd = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*( (X>=81)&(X<109) ) ...
          + (0.01521).*(X>=109);      % 50% vaccination rate starting on 1st May
      
%% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% GRAPHS
figure
subplot(3,3,1)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Do,':b','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Low Efficiency (LVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

subplot(3,3,4)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Da,':b','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Baseline Efficiency (LVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

subplot(3,3,7)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Db,':b','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('High Efficiency (LVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
%% (Normal Vaccination Rate)
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vacc = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*(X>=81);     % Vaccination rate (manually fitted to data)



%% Baseline efficacy (Normal Vaccination)
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% GRAPHS
subplot(3,3,2)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Do,':r','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Low Efficiency (NVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

subplot(3,3,5)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Da,':r','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Baseline Efficiency (NVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

subplot(3,3,8)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Db,':r','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('High Efficiency (NVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%% Defining vaccination rates:
%%High Vaccination Rate
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vaccb = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*( (X>=81)&(X<109) ) ...
          + (0.070842).*(X>=109); % 200 vaccination rate 
      
 %% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% GRAPHS
subplot(3,3,3)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Do,':g','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Low Efficiency (HVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

subplot(3,3,6)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Da,':g','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Baseline Efficiency (HVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

subplot(3,3,9)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Db,':g','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('High Efficiency (HVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%% Code for obtaining Figure S6 of the supplementary material 
%% Low Vaccination Rate
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vaccd = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*( (X>=81)&(X<109) ) ...
          + (0.01521).*(X>=109);      % 50% vaccination rate starting on 1st May
      
%% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% GRAPHS
subplot(3,3,1)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ro,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Low Efficiency (LVR)')
plot(newdata_long,RUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,4)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ra,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Baseline Efficiency (LVR)')
plot(newdata_long,RUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,7)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Rb,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('High Efficiency (LVR) ')
plot(newdata_long,AUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%% (Normal Vaccination Rate)
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vacc = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*(X>=81);     % Vaccination rate (manually fitted to data)



%% Baseline efficacy (Normal Vaccination)
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% GRAPHS
subplot(3,3,2)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ro,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Low Efficiency (NVR)')
plot(newdata_long,RUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,5)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ra,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Baseline Efficiency (NVR)')
plot(newdata_long,RUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,8)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Rb,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('High Efficiency (NVR)')
plot(newdata_long,RUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%% Defining vaccination rates:
%%High Vaccination Rate
epsilonl=0.913; epsilona=3.545e-3; epsilonla=0.75; epsilonlb=0.87;
vaccb = @(X)(0.001).*(X<19) ...
          + (0.004028).*( (X>=19)&(X<23) ) ...
          + (0.005070).*( (X>=23)&(X<53) ) ...
          + (0.011167).*( (X>=53)&(X<74) ) ...
          + (0.020129).*( (X>=74)&(X<81) ) ...
          + (0.030421).*( (X>=81)&(X<109) ) ...
          + (0.070842).*(X>=109); % 200 vaccination rate 
      
 %% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra,RUa]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.818;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro,RUo]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.907;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb,RUb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0,RU0);

%% GRAPHS
subplot(3,3,3)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ro,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Low Efficiency (HVR)')
plot(newdata_long,RUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,6)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Ra,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('Baseline Efficiency (HVR)')
plot(newdata_long,RUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,9)
tStart=datetime('2021-01-12'); tEnd=datetime('2021-10-30');
plot(newdata_long,Rb,':r','LineWidth',2)
hold on
ylabel('Original Variant')
hold on
yyaxis right
title('High Efficiency (HVR)')
plot(newdata_long,RUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Alpha Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';


%% Code for obtaining figure 4 of the article. Control Reproduction Number 
%Graph for the control reproduction number for strain a (original
%%strain)where we variate the vaccine efficiency on the reduction for
%%symptomatic individuals
delta1=1.8256e-4;
beta1=0.2;
beta2=0.0330;
gamrec1 = 0.0370;
p1=0.12;
epsilonl=0.913;
epsilonla=0.75;
epsilonlb=0.87;
epsilona=3.545e-3;
alpha=0.002739;

% x: final proportion of vaccinated people
% y: final proportion of recovered people
[X,Y]=meshgrid(0:0.01:0.7, 0:0.01:0.4);

beta1=0.8*0.2; beta2=0.8*0.0330;
epsilonl=0.89;
RC = @(x,y) ((1-p1).*[beta2.*(1-x-y)+(1-epsilonla).*beta2.*x])./gamrec1...
    + (p1.*[beta1.*(1-x-y)+(1-epsilonl).*beta1*x])./(delta1+gamrec1);
z = RC(X,Y);

figure
subplot(3,3,1)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Low transmission, \epsilon_L=0.89', 'FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y]=meshgrid(0:0.01:0.7, 0:0.01:0.4);
beta1=0.2; beta2=0.0330;
epsilonl=0.89;
RC = @(x,y) ((1-p1).*[beta2.*(1-x-y)+(1-epsilonla).*beta2.*x])./gamrec1...
    + (p1.*[beta1.*(1-x-y)+(1-epsilonl).*beta1*x])./(delta1+gamrec1);
z = RC(X,Y);
subplot(3,3,2)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Baseline transmission, \epsilon_L=0.89', 'FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1=1.2*0.2; beta2=1.2*0.0330;
epsilonl=0.89;
RC = @(x,y) ((1-p1).*[beta2.*(1-x-y)+(1-epsilonla).*beta2.*x])./gamrec1...
    + (p1.*[beta1.*(1-x-y)+(1-epsilonl).*beta1*x])./(delta1+gamrec1);
z = RC(X,Y);
subplot(3,3,3)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('High transmission, \epsilon_L=0.89','FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1=0.8*0.2; beta2=0.8*0.0330;
epsilonl=0.913;
RC = @(x,y) ((1-p1).*[beta2.*(1-x-y)+(1-epsilonla).*beta2.*x])./gamrec1...
    + (p1.*[beta1.*(1-x-y)+(1-epsilonl).*beta1*x])./(delta1+gamrec1);
z = RC(X,Y);
subplot(3,3,4)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Low transmission, \epsilon_L=0.913', 'FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1=0.2; beta2=0.0330;
epsilonl=0.913;
RC = @(x,y) ((1-p1).*[beta2.*(1-x-y)+(1-epsilonla).*beta2.*x])./gamrec1...
    + (p1.*[beta1.*(1-x-y)+(1-epsilonl).*beta1*x])./(delta1+gamrec1);
z = RC(X,Y);
subplot(3,3,5)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Baseline transmission, \epsilon_L=0.913','FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1=1.2*0.2; beta2=1.2*0.0330;
epsilonl=0.913;
RC = @(x,y) ((1-p1).*[beta2.*(1-x-y)+(1-epsilonla).*beta2.*x])./gamrec1...
    + (p1.*[beta1.*(1-x-y)+(1-epsilonl).*beta1*x])./(delta1+gamrec1);
z = RC(X,Y);
subplot(3,3,6)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',12,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('High transmission, \epsilon_L=0.913','FontSize',9)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1=0.8*0.2; beta2=0.8*0.0330;
epsilonl=0.932;
RC = @(x,y) ((1-p1).*[beta2.*(1-x-y)+(1-epsilonla).*beta2.*x])./gamrec1...
    + (p1.*[beta1.*(1-x-y)+(1-epsilonl).*beta1*x])./(delta1+gamrec1);
z = RC(X,Y);
subplot(3,3,7)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Low transmission, \epsilon_L=0.932','FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1=0.2; beta2=0.0330;
epsilonl=0.932;
RC = @(x,y) ((1-p1).*[beta2.*(1-x-y)+(1-epsilonla).*beta2.*x])./gamrec1...
    + (p1.*[beta1.*(1-x-y)+(1-epsilonl).*beta1*x])./(delta1+gamrec1);
z = RC(X,Y);
subplot(3,3,8)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Baseline transmission, \epsilon_L=0.932','FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1=1.2*0.2; beta2=1.2*0.0330;
epsilonl=0.932;
RC = @(x,y) ((1-p1).*[beta2.*(1-x-y)+(1-epsilonla).*beta2.*x])./gamrec1...
    + (p1.*[beta1.*(1-x-y)+(1-epsilonl).*beta1*x])./(delta1+gamrec1);
z = RC(X,Y);
subplot(3,3,9)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('High transmission, \epsilon_L=0.932','FontSize',10)

%% Code for figure 5 
%%Graph for the control reproduction number for strain b (UK
%%strain)where we variate the vaccine efficiency on the reduction for
%%symptomatic individuals
delta1=1.8256e-4;
beta1=0.3;
beta2=0.0495;
gamrec1 = 0.0370;
p1=0.12;
epsilonl=0.913;
epsilonla=0.75;
epsilonlb=0.87;
epsilona=3.545e-3;
alpha=0.002739;
q=0.2;
delta2=1.8256e-4;
gamrec2=0.0370;

% x: final proportion of vaccinated people
% y: final proportion of recovered people
[X,Y]=meshgrid(0:0.01:0.9, 0:0.01:0.8);

beta3=0.8*0.3; beta4=0.8*0.0495;
epsilonlb=0.818;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);


subplot(3,3,1)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Low transmission, \epsilon_L_B=0.818','FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=0.3; beta4=0.0495;
epsilonlb=0.818;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,2)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Baseline transmission, \epsilon_L_B=0.818','FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=1.2*0.3; beta4=1.2*0.0495;
epsilonlb=0.818;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,3)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('High transmission, \epsilon_L_B=0.818','FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=0.8*0.3; beta4=0.8*0.0495;
epsilonlb=0.87;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,4)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Low transmission, \epsilon_L_B=0.87','FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=0.3; beta4=0.0495;
epsilonlb=0.87;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,5)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Baseline transmission, \epsilon_L_B=0.87','FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=1.2*0.3; beta4=1.2*0.0495;
epsilonlb=0.87;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,6)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('High transmission, \epsilon_L_B=0.87','FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=0.8*0.3; beta4=0.8*0.0495;
epsilonlb=0.907;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,7)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Low transmission, \epsilon_L_B=0.907','FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=0.3; beta4=0.0495;
epsilonlb=0.907;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,8)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Baseline transmission rate, \epsilon_L_B=0.907','FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=1.2*0.3; beta4=1.2*0.0495;
epsilonlb=0.907;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,9)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('High transmission, \epsilon_L_B=0.907','FontSize',10)

%% Figura 6A Local Sensitivity Analysis for the original variant 
values=[0.29981 0.700184 -0.9985 -0.001472 0.2043 -1.54425 -1.5454 0.0009083 -0.2553 0.2553];
figure
subplot(2,1,1)
bh = bar(1:numel(values),diag(values),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10

title('Local Sensitivity Analysis of the Control Reproduction Number of the Wild-type variant')
set(gca,'XTick', 1:numel(values),'XTickLabel',{'\beta_1','\beta_2','\gamma_1','\delta_1','p','\epsilon_L_A','\epsilon_L','\epsilon_a','\rho','\alpha'});

%%Figure 6B Sensitivity Analysis alpha strain
values=[0.60123 0.39877 -0.9952 -0.0047 0.50153 -3.9534 -.001162 -0.3266 0.3266];
subplot(2,1,2)
bh = bar(1:numel(values),diag(values),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9


title('Local Sensitivity Analysis of the Control Reproduction Number of the Alpha variant')
set(gca,'XTick', 1:numel(values),'XTickLabel',{'\beta_3','\beta_4','\gamma_2','\delta_2','q','\epsilon_L_B','\epsilon_a','\rho','\alpha'});


%% Global Sensitivity Analysis of the wild-type vs alpha strain 
clear
clc
%   Intervalos
rangmin=[0.2 0.096 0.16 0.16 0.00264 0.24 0.0396 1.4605e-4 1.4605e-5 0.0296 0.024 0.7304 0.6 0.696 0.0028 0.0243 0.0022 0.0022];
rangmax=[0.3 0.144 0.24 0.24 0.0396 0.36 0.0594 2.1907e-4 2.1907e-5 0.0444 0.036 0.97 0.9 0.94 0.0043 0.0365 0.0033 0.033];
%   sampling
M=HCV3_lhs(1000,rangmin,rangmax);
PobS=zeros(1000,1);
PobE=zeros(1000,1);
PobI=zeros(1000,1);
PobA=zeros(1000,1);
PobV=zeros(1000,1);
PobEU=zeros(1000,1);
PobIU=zeros(1000,1);
PobAU=zeros(1000,1);
PobD=zeros(1000,1);
PobR=zeros(1000,1);
%   1000 runs
for i=1:1000;
    HCV3_odeSamplingTwoStrain(M(i,1),M(i,2),M(i,3),M(i,4),M(i,5),M(i,6),M(i,7),M(i,8),M(i,9),M(i,10),M(i,11),M(i,12),M(i,13),M(i,14),M(i,15),M(i,16),M(i,17),M(i,18))
    load('qfile.mat');
    PobS(i)=p1;
    PobE(i)=p2;
    PobI(i)=p3;
    PobA(i)=p4;
    PobV(i)=p5;
    PobEU(i)=p6;
    PobIU(i)=p7;
    PobAU(i)=p8;
    PobD(i)=p9;
    PobR(i)=p10;
    end
%%% Ahora DD1 es un vector cuyas entradas son la población de hepatocitos sanos en t=50
%%% y las entradas de DD2 son la población de hepatocitos infectados en t=50
%   Rank transform
for i=1:18; %%%asdasd
    M(:,i)=tiedrank(M(:,i));
end
   PobS=tiedrank(PobS);
    PobE=tiedrank(PobE);
    PobI=tiedrank(PobI);
    PobA=tiedrank(PobA);
    PobV=tiedrank(PobV);
    PobEU=tiedrank(PobEU);
    PobIU=tiedrank(PobIU);
    PobAU=tiedrank(PobAU);
    PobD=tiedrank(PobD);
    PobR=tiedrank(PobR);
%   Empezando las prcc
MS=zeros(1000,21); %%%% asdfsdf +2
ME=zeros(1000,21);
MI=zeros(1000,21);
MA=zeros(1000,21);
MP=zeros(1000,21);
MF=zeros(1000,21);
MEB=zeros(1000,21);
MIB=zeros(1000,21);
MAB=zeros(1000,21);
MR=zeros(1000,21);
MD=zeros(1000,21);
for i=1:18
    MS(:,i)=M(:,i);
    ME(:,i)=M(:,i);
    MI(:,i)=M(:,i);
    MA(:,i)=M(:,i);
    MV(:,i)=M(:,i);
    MEU(:,i)=M(:,i);
    MIU(:,i)=M(:,i);
    MAU(:,i)=M(:,i);
    MD(:,i)=M(:,i);
    MR(:,i)=M(:,i);
end
%   las variables dummy;
MS(:,19)=randperm(1000,1000)';
ME(:,19)=MS(:,19);
MI(:,19)=MS(:,19);
MA(:,19)=MS(:,19);
MV(:,19)=MS(:,19);
MEU(:,19)=MS(:,19);
MIU(:,19)=MS(:,19);
MAU(:,19)=MS(:,19);
MD(:,19)=MS(:,19);
MR(:,19)=MS(:,19);%%%% asdsac +1
%   terminando de crear la matriz para el prcc
MS(:,21)=PobS; %%%asdasdas+2
ME(:,21)=PobE;
MI(:,21)=PobI; %%%asdasdas+2
MA(:,21)=PobA;
MV(:,21)=PobV; %%%asdasdas+2
MEU(:,21)=PobEU; %%%asdasdas+2
MIU(:,21)=PobIU;
MAU(:,21)=PobAU; %%%asdasdas+2
MD(:,21)=PobD;
MR(:,21)=PobR;
%   PRCC
[PCS,PvS]=partialcorr(MS);
[PCE,PvE]=partialcorr(ME);
[PCI,PvI]=partialcorr(MI);
[PCA,PvA]=partialcorr(MA);
[PCV,PvV]=partialcorr(MV);
[PCEU,PvEU]=partialcorr(MEU);
[PCIU,PvIU]=partialcorr(MIU);
[PCAU,PvAU]=partialcorr(MAU);
[PCD,PvD]=partialcorr(MD);
[PCR,PvR]=partialcorr(MR);
%   bar graph
GS=zeros(19,1);
GE=zeros(19,1);
GI=zeros(19,1);
GA=zeros(19,1);
GV=zeros(19,1);
GEU=zeros(19,1);
GIU=zeros(19,1);
GAU=zeros(19,1);
GD=zeros(19,1);
GR=zeros(19,1);
GvS=zeros(19,1);
GvE=zeros(19,1);
GvI=zeros(19,1);
GvA=zeros(19,1);
GvV=zeros(19,1);
GvEU=zeros(19,1);
GvIU=zeros(19,1);
GvAU=zeros(19,1);
GvD=zeros(19,1);
GvR=zeros(19,1);
for i=1:20
    GS(i)=PCS(i,21);
    GE(i)=PCE(i,21);
    GI(i)=PCI(i,21);
    GA(i)=PCA(i,21);
    GV(i)=PCV(i,21);
    GEU(i)=PCEU(i,21);
    GIU(i)=PCIU(i,21);
    GAU(i)=PCAU(i,21);
    GD(i)=PCD(i,21);
    GR(i)=PCR(i,21);
    GvS(i)=PvS(i,21);
    GvE(i)=PvE(i,21);
    GvI(i)=PvI(i,21);
    GvA(i)=PvA(i,21);
    GvV(i)=PvV(i,21);
    GvEU(i)=PvEU(i,21);
    GvIU(i)=PvIU(i,21);
    GvAU(i)=PvAU(i,21);
    GvD(i)=PvD(i,21);
    GvR(i)=PvR(i,21);
end

%%Susceptible
figure
subplot(2,1,1)
bh= bar(1:numel(GS),diag(GS),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis Susceptible Individuals')
set(gca,'XTick', 1:numel(GS),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvS),diag(GvS),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
set(gca,'XTick', 1:numel(GvS),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
title('P values of the parameters')


%% Figure for infected individuals 
figure
subplot(2,1,1)
bh= bar(1:numel(GI),diag(GI),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Infected Symptomatic Individuals to the Wild-type variant')
set(gca,'XTick', 1:numel(GI),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvI),diag(GvI),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvI),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});

%% Figure for Asymptomatic individuals 
figure
subplot(2,1,1)
bh= bar(1:numel(GA),diag(GA),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Infected Asymptomatic Individuals to the Wild-type variant')
set(gca,'XTick', 1:numel(GA),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvA),diag(GvA),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvA),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});

%% Figure for vaccinated individuals 
figure
subplot(2,1,1)
bh= bar(1:numel(GV),diag(GV),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Vaccinated Individuals')
set(gca,'XTick', 1:numel(GV),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvV),diag(GvV),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvV),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});

%% Figure for exposed to the second variant 
figure
subplot(2,1,1)
bh= bar(1:numel(GEU),diag(GEU),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Exposed Individuals to the Alpha Variant')
set(gca,'XTick', 1:numel(GEU),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvEU),diag(GvEU),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvEU),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});

%% Figure for Infected of the second variant 
figure
subplot(2,1,1)
bh= bar(1:numel(GIU),diag(GIU),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Infected Symptomatic Individuals to the Alpha Variant')
set(gca,'XTick', 1:numel(GIU),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvIU),diag(GvIU),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvIU),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});

%% Figure for Asymptomatic infections with the second variant 
figure
subplot(2,1,1)
bh= bar(1:numel(GAU),diag(GAU),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Infected Asymptomatic Individuals to the Alpha Variant')
set(gca,'XTick', 1:numel(GAU),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvAU),diag(GvAU),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvAU),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});

%% Figure for death individuals 
figure
subplot(2,1,1)
bh= bar(1:numel(GD),diag(GD),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Diseased Individuals')
set(gca,'XTick', 1:numel(GD),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvD),diag(GvD),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvD),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});

%% Recovered individuals 
figure
subplot(2,1,1)
bh= bar(1:numel(GR),diag(GR),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Recovered Individuals')
set(gca,'XTick', 1:numel(GR),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvR),diag(GvR),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvR),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});



%% Graphs for the behaviour of the alpha variant once delta appeared in the US 
%% parameters 
global w p1 N beta1 beta2 delta1 gamrec1 beta3 beta4 gamrec2 q alpha epsilona epsilonl epsilonla epsilonlb delta2 eta
N=331449281;
beta1=0.3;
beta2=0.0495;
beta3=0.5;
beta4=0.08;
epsilona=0.0862;
alpha=0.005555;
w=0.25;
p1=0.12;
delta1=1.8256e-4;
gamrec1=0.037;
gamrec2=0.03;
epsilonl=0.87;
epsilonla=0.75;
epsilonlb=0.49;
q=0.2;
delta2=1.8256e-3;
eta=0.002739;

%% Seting initial conditions based on 'simul-before-vaccine.mat'
S0=331449281;
E0 = 25000000;
I0 = 25000000;
A0 = ((1-0.2)/0.2)*I0;
V0=140000000;
EU0 = 100000;
IU0 = 100000;
AU0 = ((1-0.2)/0.2)*IU0;
R0 = 1.1125e+06;
D0 = 621000;
str = ['2021-06-01'];
date = datetime(str,'InputFormat','yyyy-MM-dd');

%% Solving the model with the given parameters

finaldate = 700;            %Date until which we want to make the predictions
t=linspace(0,finaldate,finaldate+1);
newdata_long=dateshift(date,'start','day',0:t(end));

%% Defining vaccination rates:
novac = @(X) 0.0;   % No vaccination
%% (Normal Vaccination Rate)
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555; 
vacc = @(X)(0.030421).*(X>=1);     % Vaccination rate (manually fitted to data)
%% Defining vaccination rates:
%%
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555;
vaccb = @(X)(0.070842).*(X>=1); % 200 vaccination rate 
%% Low Vaccination Rate
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555;
vaccd = @(X)(0.01521).*(X>=1);      % 50% vaccination rate starting on 1st May

%% LOW TRANSMISSION RATE
beta1 = 0.3*0.8;
beta2 = 0.0495*0.8;
beta3=0.5*0.8;
beta4=0.08*0.8;
[t,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(novac,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sc,Ec,Ic,Ac,Vc,EUc,IUc,AUc,Dc,Rc]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%%Figure 1 
figure
subplot(3,2,1)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Io+IUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Ia+IUa,'-b','LineWidth',2)
plot(newdata_long,Ib+IUb,'-.k','LineWidth',2)
plot(newdata_long,Ic+IUc, '-.g', 'LineWidth',2)
title('Low Transmission Rate')
ylabel('I(t)+IU(t)')
xlim([tStart tEnd])
%ylim([0 1e6])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate', '50% vaccination rate'},'Location','northwest')
grid on

subplot(3,2,2)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ao+AUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Aa+AUa,'-b','LineWidth',2)
plot(newdata_long,Ab+AUb,'-.k','LineWidth',2)
plot(newdata_long,Ac+AUc, '-.g', 'LineWidth',2)
title('Low Transmission Rate')
ylabel('A(t)+AU(t)')
xlim([tStart tEnd])
%ylim([0 1e6])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate', '50% vaccination rate'},'Location','northwest')
grid on

%% Baseline transmission rate. 
global w p1 N beta1 beta2 delta1 gamrec1 beta3 beta4 gamrec2 q alpha epsilona epsilonl epsilonla epsilonlb delta2 eta
N=331449281;
beta1=0.3;
beta2=0.0495;
beta3=0.5;
beta4=0.08;
epsilona=0.0862;
alpha=0.005555;
w=0.25;
p1=0.12;
delta1=1.8256e-4;
gamrec1=0.037;
gamrec2=0.03;
epsilonl=0.87;
epsilonla=0.75;
epsilonlb=0.49;
q=0.2;
delta2=1.8256e-3;
eta=0.002739;

%% Seting initial conditions based on 'simul-before-vaccine.mat'
S0=331449281;
E0 = 25000000;
I0 = 25000000;
A0 = ((1-0.2)/0.2)*I0;
V0=140000000;
EU0 = 100000;
IU0 = 100000;
AU0 = ((1-0.2)/0.2)*IU0;
R0 = 1.1125e+06;
D0 = 621000;
str = ['2021-06-01'];
date = datetime(str,'InputFormat','yyyy-MM-dd');

%% Solving the model with the given parameters

finaldate = 700;            %Date until which we want to make the predictions
t=linspace(0,finaldate,finaldate+1);
newdata_long=dateshift(date,'start','day',0:t(end));

%% Defining vaccination rates:
novac = @(X) 0.0;   % No vaccination
%% (Normal Vaccination Rate)
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555; 
vacc = @(X)(0.030421).*(X>=1);     % Vaccination rate (manually fitted to data)
%% Defining vaccination rates:
%%
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555;
vaccb = @(X)(0.070842).*(X>=1); % 200 vaccination rate 
%% Low Vaccination Rate
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555;
vaccd = @(X)(0.01521).*(X>=1);      % 50% vaccination rate starting on 1st May
beta1 = 0.3;
beta2 = 0.0495;
beta3=0.5;
beta4=0.08;
[t,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(novac,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sc,Ec,Ic,Ac,Vc,EUc,IUc,AUc,Dc,Rc]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

subplot(3,2,3)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Io+IUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Ia+IUa,'-b','LineWidth',2)
plot(newdata_long,Ib+IUb,'-.k','LineWidth',2)
plot(newdata_long,Ic+IUc,'-.g','LineWidth',2)
title('Normal transmission rate')
ylabel('I(t)+IU(t)')
xlim([tStart tEnd])
%ylim([0 1e6])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate','50% vaccination rate' },'Location','northwest')
grid on

subplot(3,2,4)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ao+AUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Aa+AUa,'-b','LineWidth',2)
plot(newdata_long,Ab+AUb,'-.k','LineWidth',2)
plot(newdata_long,Ac+AUc,'-.g','LineWidth',2)
title('Normal Transmission Rate')
ylabel('A(t)+AU(t)')
xlim([tStart tEnd])
%ylim([0 1e6])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate','50% vaccination rate'},'Location','northwest')
grid on

%% High Transmission Rate
beta1 = 0.3*1.5;
beta2 = 0.0495*1.5;
beta3=0.5*1.5;
beta4=0.08*1.5;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(novac,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sc,Ec,Ic,Ac,Vc,EUc,IUc,AUc,Dc,Rc]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);


subplot(3,2,5)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Io+IUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Ia+IUa,'-b','LineWidth',2)
plot(newdata_long,Ib+IUb,'-.k','LineWidth',2)
plot(newdata_long,Ic+IUc,'-.g','LineWidth',2)
title('High transmission rate')
ylabel('I(t)+IU(t)')
xlim([tStart tEnd])
%ylim([0 1e6])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate','50% vaccination rate' },'Location','northwest')
grid on

subplot(3,2,6)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ao+AUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Aa+AUa,'-b','LineWidth',2)
plot(newdata_long,Ab+AUb,'-.k','LineWidth',2)
plot(newdata_long,Ac+AUc,'-.g','LineWidth',2)
title('High transmission rate')
ylabel('A(t)+AU(t)')
xlim([tStart tEnd])
%ylim([0 1e6])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate','50% vaccination rate' },'Location','northwest')
grid on


%%Code for obtaining figure 2 of the supplementary material 
%% LOW TRANSMISSION RATE
beta1 = 0.3*0.8;
beta2 = 0.0495*0.8;
beta3=0.5*0.8;
beta4=0.08*0.8;
[t,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(novac,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sc,Ec,Ic,Ac,Vc,EUc,IUc,AUc,Dc,Rc]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

figure
subplot(3,2,1)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ro,'--r','LineWidth',2)
hold on
plot(newdata_long,Ra,'-b','LineWidth',2)
plot(newdata_long,Rb,'-.k','LineWidth',2)
plot(newdata_long,Rc,'-.g', 'LineWidth',2)
title('Low Transmission Rate')
ylabel('R(t)')
xlim([tStart tEnd])
%ylim([1e5 9e5])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate', '50% vaccination rate'},'Location','southeast')
grid on

subplot(3,2,2)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Do,'--r','LineWidth',2)
hold on
plot(newdata_long,Da,'-b','LineWidth',2)
plot(newdata_long,Db,'-.k','LineWidth',2)
plot(newdata_long,Dc,'-.g', 'LineWidth',2)
title('Low Transmission Rate')
ylabel('D(t)')
xlim([tStart tEnd])
%ylim([1e5 9e5])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate', '50% vaccination rate'},'Location','southeast')
grid on
%% Baseline transmission rate. 
beta1 = 0.3;
beta2 = 0.0495;
beta3=0.5;
beta4=0.08;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(novac,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sc,Ec,Ic,Ac,Vc,EUc,IUc,AUc,Dc,Rc]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

subplot(3,2,3)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ro,'--r','LineWidth',2)
hold on
plot(newdata_long,Ra,'-b','LineWidth',2)
plot(newdata_long,Rb,'-.k','LineWidth',2)
plot(newdata_long,Rc,'-.g','LineWidth',2)
ylabel('R(t)')
title('Normal Transmission Rate')
xlim([tStart tEnd])
%ylim([1e5 9e5])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate','50% vaccination rate'},'Location','northwest')
grid on

subplot(3,2,4)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Do,'--r','LineWidth',2)
hold on
plot(newdata_long,Da,'-b','LineWidth',2)
plot(newdata_long,Db,'-.k','LineWidth',2)
plot(newdata_long,Dc,'-.g', 'LineWidth',2)
title('Normal Transmission Rate')
ylabel('D(t)')
xlim([tStart tEnd])
%ylim([1e5 9e5])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate', '50% vaccination rate'},'Location','southeast')
grid on

%% High Transmission Rate
beta1 = 0.3*1.5;
beta2 = 0.0495*1.5;
beta3=0.5*1.5;
beta4=0.08*1.5;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(novac,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);
[t,Sc,Ec,Ic,Ac,Vc,EUc,IUc,AUc,Dc,Rc]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

subplot(3,2,5)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ro,'--r','LineWidth',2)
hold on
plot(newdata_long,Ra,'-b','LineWidth',2)
plot(newdata_long,Rb,'-.k','LineWidth',2)
plot(newdata_long,Rc,'-.g','LineWidth',2)
title('High transmission rate')
ylabel('R(t)')
xlim([tStart tEnd])
%ylim([1e5 9e5])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate','50% vaccination rate'},'Location','northwest')
grid on

subplot(3,2,6)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Do,'--r','LineWidth',2)
hold on
plot(newdata_long,Da,'-b','LineWidth',2)
plot(newdata_long,Db,'-.k','LineWidth',2)
plot(newdata_long,Dc,'-.g', 'LineWidth',2)
title('High Transmission Rate')
ylabel('D(t)')
xlim([tStart tEnd])
%ylim([1e5 9e5])
legend({'No vaccination','Baseline vaccination rate','200% vaccination rate', '50% vaccination rate'},'Location','southeast')
grid on


% Plotting the solutions of the model using different values for the
% efficacy rates eta1 and eta2

global w p1 N beta1 beta2 delta1 gamrec1 beta3 beta4 gamrec2 q alpha epsilona epsilonl epsilonla epsilonlb delta2 eta
N=331449281;
beta1=0.3;
beta2=0.0495;
beta3=0.5;
beta4=0.08;
epsilona=0.0862;
alpha=0.005555;
w=0.25;
p1=0.12;
delta1=1.8256e-4;
gamrec1=0.037;
gamrec2=0.03;
epsilonl=0.87;
epsilonla=0.75;
epsilonlb=0.49;
q=0.2;
delta2=1.8256e-3;
eta=0.002739;

%% Seting initial conditions based on 'simul-before-vaccine.mat'
S0=331449281;
E0 = 25000000;
I0 = 25000000;
A0 = ((1-0.2)/0.2)*I0;
V0=140000000;
EU0 = 100000;
IU0 = 100000;
AU0 = ((1-0.2)/0.2)*IU0;
R0 = 1.1125e+06;
D0 = 621000;
str = ['2021-06-01'];
date = datetime(str,'InputFormat','yyyy-MM-dd');

%% Solving the model with the given parameters

finaldate = 700;            %Date until which we want to make the predictions
t=linspace(0,finaldate,finaldate+1);
newdata_long=dateshift(date,'start','day',0:t(end));
%% Low Vaccination Rate
%% Low Vaccination Rate
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555;
vaccd = @(X)(0.01521).*(X>=1);      % 50% vaccination rate starting on 1st May   % 50% vaccination rate starting on 1st May
      
%% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% GRAPHS
%%Code for obtaininf figure S3 of the supplementary material
figure
subplot(3,2,1)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Io+IUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Ia+IUa,'-b','LineWidth',2)
plot(newdata_long,Ib+IUb,'-.k','LineWidth',2)
title('Symptomatic Infected Individuals')
ylabel('I(t)+IU(t)')
xlim([tStart tEnd])
legend({'\epsilon_a=0.0067, \epsilon_La=0.72, \epsilon_Lb=0.22, \epsilon_L=0.89','\epsilon_a=3.545e-3, \epsilon_La=0.75, \epsilon_Lb=0.49, \epsilon_L=0.913', '\epsilon_a=3.77e-4, \epsilon_La=0.95, \epsilon_Lb=0.67, \epsilon_L=0.932'},'Location','northwest')
grid on

subplot(3,2,2)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ao+AUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Aa+AUa,'-b','LineWidth',2)
plot(newdata_long,Ab+AUb,'-.k','LineWidth',2)
title('Asymptomatic Infected Individuals')
ylabel('A(t)+AU(t)')
xlim([tStart tEnd])
legend({'\epsilon_a=0.0067, \epsilon_La=0.72, \epsilon_Lb=0.22, \epsilon_L=0.89','\epsilon_a=3.545e-3, \epsilon_La=0.75, \epsilon_Lb=0.49, \epsilon_L=0.913', '\epsilon_a=3.77e-4, \epsilon_La=0.95, \epsilon_Lb=0.67, \epsilon_L=0.932'},'Location','northwest')
grid on



%% (Normal Vaccination Rate)
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555; 
vacc = @(X)(0.030421).*(X>=1);  % Vaccination rate (manually fitted to data)



%% Baseline efficacy (Normal Vaccination)
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%%Code for obtaining figure S3 of the supplementary material

subplot(3,2,3)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Io+IUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Ia+IUa,'-b','LineWidth',2)
plot(newdata_long,Ib+IUb,'-.k','LineWidth',2)
title('Symptomatic Infected Individuals')
ylabel('I(t)+ IU(t)')
xlim([tStart tEnd])
legend({'\epsilon_a=0.0067, \epsilon_La=0.72, \epsilon_Lb=0.22, \epsilon_L=0.89','\epsilon_a=3.545e-3, \epsilon_La=0.75, \epsilon_Lb=0.49, \epsilon_L=0.913', '\epsilon_a=3.77e-4, \epsilon_La=0.95, \epsilon_Lb=0.67, \epsilon_L=0.932'},'Location','northwest')
grid on

subplot(3,2,4)
plot(newdata_long,Ao+AUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Aa+AUa,'-b','LineWidth',2)
plot(newdata_long,Ab+AUb,'-.k','LineWidth',2)
title('Asymptomatic Infected Individuals')
ylabel('A(t)+AU(t)')
xlim([tStart tEnd])
legend({'\epsilon_a=0.0067, \epsilon_La=0.72, \epsilon_Lb=0.22, \epsilon_L=0.89','\epsilon_a=3.545e-3, \epsilon_La=0.75, \epsilon_Lb=0.49, \epsilon_L=0.913', '\epsilon_a=3.77e-4, \epsilon_La=0.95, \epsilon_Lb=0.67, \epsilon_L=0.932'},'Location','northwest')
grid on

%% High Vaccination Rate
%%
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555;
vaccb = @(X)(0.070842).*(X>=1); %200 vaccination rate 
      
 %% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% GRAPHS

subplot(3,2,5)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Io+IUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Ia+IUa,'-b','LineWidth',2)
plot(newdata_long,Ib+IUb,'-.k','LineWidth',2)
title('Symptomatic Infected Individuals')
ylabel('I(t)+IU(t)')
xlim([tStart tEnd])
legend({'\epsilon_a=0.0067, \epsilon_La=0.72, \epsilon_Lb=0.22, \epsilon_L=0.89','\epsilon_a=3.545e-3, \epsilon_La=0.75, \epsilon_Lb=0.49, \epsilon_L=0.913', '\epsilon_a=3.77e-4, \epsilon_La=0.95, \epsilon_Lb=0.67, \epsilon_L=0.932'},'Location','northwest')
grid on

subplot(3,2,6)
plot(newdata_long,Ao+AUo,'--r','LineWidth',2)
hold on
plot(newdata_long,Aa+AUa,'-b','LineWidth',2)
plot(newdata_long,Ab+AUb,'-.k','LineWidth',2)
title('Asymptomatic Infected Individuals')
ylabel('A(t)+AU(t)')
xlim([tStart tEnd])
legend({'\epsilon_a=0.0067, \epsilon_La=0.72, \epsilon_Lb=0.22, \epsilon_L=0.89','\epsilon_a=3.545e-3, \epsilon_La=0.75, \epsilon_Lb=0.49, \epsilon_L=0.913', '\epsilon_a=3.77e-4, \epsilon_La=0.95, \epsilon_Lb=0.67, \epsilon_L=0.932'},'Location','northwest')
grid on

%% Low Vaccination Rate
%% Low Vaccination Rate
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555;
vaccd = @(X)(0.01521).*(X>=1);      % 50% vaccination rate starting on 1st May   % 50% vaccination rate starting on 1st May
      
%% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);


%% GRAPHS
subplot(3,3,1)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Io,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('Low Efficiency (LVR)')
plot(newdata_long,IUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,4)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ia,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('Baseline Efficiency (LVR)')
plot(newdata_long,IUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,7)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ib,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('High Efficiency (LVR) ')
plot(newdata_long,IUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%% (Normal Vaccination Rate)
%% (Normal Vaccination Rate)
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555; 
vacc = @(X)(0.030421).*(X>=1);  % Vaccination rate (manually fitted to data)



%% Baseline efficacy (Normal Vaccination)
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% GRAPHS
subplot(3,3,2)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Io,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('Low Efficiency (NVR)')
plot(newdata_long,IUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,5)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ia,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('Baseline Efficiency (NVR)')
plot(newdata_long,IUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,8)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ib,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('High Efficiency (NVR)')
plot(newdata_long,IUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';


%% High Vaccination Rate
%%
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555;
vaccb = @(X)(0.070842).*(X>=1); %200 vaccination rate 
      
 %% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% GRAPHS
subplot(3,3,3)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Io,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('Low Efficiency (HVR)')
plot(newdata_long,IUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,6)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ia,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('Baseline Efficiency (HVR)')
plot(newdata_long,IUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,9)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ib,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('High Efficiency (HVR)')
plot(newdata_long,IUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';


%%Code for obtaining figure S4 

%% Low Vaccination Rate
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555;
vaccd = @(X)(0.01521).*(X>=1);      % 50% vaccination rate starting on 1st May   % 50% vaccination rate starting on 1st May
      
%% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);


%% GRAPHS
subplot(3,3,1)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ao,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('Low Efficiency (LVR)')
plot(newdata_long,AUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,4)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Aa,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('Baseline Efficiency (LVR)')
plot(newdata_long,AUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,7)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ab,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('High Efficiency (LVR) ')
plot(newdata_long,AUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%% (Normal Vaccination Rate)
%% (Normal Vaccination Rate)
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555; 
vacc = @(X)(0.030421).*(X>=1);  % Vaccination rate (manually fitted to data)



%% Baseline efficacy (Normal Vaccination)
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% GRAPHS
subplot(3,3,2)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ao,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('Low Efficiency (NVR)')
plot(newdata_long,IUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,5)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Aa,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('Baseline Efficiency (NVR)')
plot(newdata_long,AUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,8)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ab,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('High Efficiency (NVR)')
plot(newdata_long,AUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';


%% High Vaccination Rate
%%
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555;
vaccb = @(X)(0.070842).*(X>=1); %200 vaccination rate 
      
 %% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% GRAPHS
subplot(3,3,3)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ao,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('Low Efficiency (HVR)')
plot(newdata_long,AUo,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Part b
subplot(3,3,6)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Aa,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('Baseline Efficiency (HVR)')
plot(newdata_long,AUa,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%Part c
subplot(3,3,9)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ab,':r','LineWidth',2)
hold on
ylabel('Alpha Variant')
hold on
yyaxis right
title('High Efficiency (HVR)')
plot(newdata_long,AUb,'-.g','LineWidth',2)
xlim([tStart tEnd])
%ylim([0 4e6])
ylabel('Delta Variant')
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';


%Part b%%Code for obtaining figure S5 
%% Low Vaccination Rate
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555;
vaccd = @(X)(0.01521).*(X>=1);      % 50% vaccination rate starting on 1st May   % 50% vaccination rate starting on 1st May
      
%% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);


%% GRAPHS
subplot(3,3,1)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Do,':g','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Low Efficiency (LVR)')
xlim([tStart tEnd])
%ylim([0 1e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
subplot(3,3,4)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Da,':g','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Baseline Efficiency (LVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%Part c
subplot(3,3,7)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Db,':g','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('High Efficiency (LVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%% (Normal Vaccination Rate)
%% (Normal Vaccination Rate)
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555; 
vacc = @(X)(0.030421).*(X>=1);  % Vaccination rate (manually fitted to data)



%% Baseline efficacy (Normal Vaccination)
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% GRAPHS
subplot(3,3,2)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Do,':b','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Baseline Efficiency (NVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%Part b
subplot(3,3,5)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Da,':b','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Baseline Efficiency (NVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%Part c
subplot(3,3,8)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Db,':b','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('High Efficiency (NVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%% High Vaccination Rate
%%
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555;
vaccb = @(X)(0.070842).*(X>=1); %200 vaccination rate 
      
 %% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% GRAPHS
subplot(3,3,3)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Do,':r','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Baseline Efficiency (HVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%Part b
subplot(3,3,6)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Da,':r','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Baseline Efficiency (HVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%Part c
subplot(3,3,9)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Db,':r','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('High Efficiency (HVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';


%% Code for obtaining Figure S6 of the supplementary material 
%% Low Vaccination Rate
%%Code for obtaining figure S5 
%% Low Vaccination Rate
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555;
vaccd = @(X)(0.01521).*(X>=1);      % 50% vaccination rate starting on 1st May   % 50% vaccination rate starting on 1st May
      
%% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccd,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);


%% GRAPHS
subplot(3,3,1)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ro,':g','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Low Efficiency (LVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
%Part b
subplot(3,3,4)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ra,':g','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Baseline Efficiency (LVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%Part c
subplot(3,3,7)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Rb,':g','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('High Efficiency (LVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%% (Normal Vaccination Rate)
%% (Normal Vaccination Rate)
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555; 
vacc = @(X)(0.030421).*(X>=1);  % Vaccination rate (manually fitted to data)



%% Baseline efficacy (Normal Vaccination)
[t,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vacc,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% GRAPHS
subplot(3,3,2)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ro,':b','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Baseline Efficiency (NVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%Part b
subplot(3,3,5)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ra,':b','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Baseline Efficiency (NVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%Part c
subplot(3,3,8)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Rb,':b','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('High Efficiency (NVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%% High Vaccination Rate
%%
epsilonl=0.87; epsilona=0.0862; epsilonla=0.75; epsilonlb=0.49; alpha=0.005555;
vaccb = @(X)(0.070842).*(X>=1); %200 vaccination rate 
      
 %% Baseline efficacy (High Vaccination)
[~,Sa,Ea,Ia,Aa,Va,EUa,IUa,AUa,Da,Ra]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% Low efficacy
epsilona = 0.0067;
epsilonla = 0.72;
epsilonlb=0.22;
epsilonl=0.89;
[~,So,Eo,Io,Ao,Vo,EUo,IUo,AUo,Do,Ro]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% High efficacy
epsilona = 3.77e-4;
epsilonla = 0.95;
epsilonlb=0.67;
epsilonl=0.932;
[t,Sb,Eb,Ib,Ab,Vb,EUb,IUb,AUb,Db,Rb]=vaccine_solver_us2(vaccb,t,S0,E0,I0,A0,V0,EU0,IU0,AU0,D0,R0);

%% GRAPHS
subplot(3,3,3)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ro,':r','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Baseline Efficiency (HVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%Part b
subplot(3,3,6)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Ra,':r','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('Baseline Efficiency (HVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';

%Part c
subplot(3,3,9)
tStart=datetime('2021-06-01'); tEnd=datetime('2021-11-30');
plot(newdata_long,Rb,':r','LineWidth',2)
hold on
ylabel('D(t)')
hold on
title('High Efficiency (HVR)')
xlim([tStart tEnd])
%ylim([0 4e6])
hold off
grid on
ax = gca;
ax.YAxis(1).Color = 'k';


%% Control Reproduction Number for delta variant 
%%Graph for the control reproduction number for strain c (delta
%%strain)where we variate the vaccine efficiency on the reduction for
%%symptomatic individuals
beta1=0.3;
beta2=0.0495;
beta3=0.5;
beta4=0.08;
epsilona=0.0862;
alpha=0.005555;
w=0.25;
p1=0.12;
delta1=1.8256e-4;
gamrec1=0.037;
gamrec2=0.03;
epsilonl=0.87;
epsilonla=0.75;
epsilonlb=0.49;
q=0.2;
delta2=1.8256e-3;
eta=0.002739;

% x: final proportion of vaccinated people
% y: final proportion of recovered people
[X,Y]=meshgrid(0:0.01:1.0, 0:0.01:0.8);

beta3=0.8*0.5; beta4=0.8*0.08;
epsilonlb=0.22;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);

figure
subplot(3,3,1)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Low transmission, \epsilon_L_B=0.22')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=0.5; beta4=0.08;
epsilonlb=0.22;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,2)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Baseline transmission, \epsilon_L_B=0.22')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=1.2*0.5; beta4=1.2*0.08;
epsilonlb=0.22;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,3)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('High transmission, \epsilon_L_B=0.22')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=0.8*0.5; beta4=0.8*0.08;
epsilonlb=0.49;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,4)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Low transmission, \epsilon_L_B=0.49')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=0.5; beta4=0.08;
epsilonlb=0.49;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,5)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Baseline transmission, \epsilon_L_B=0.49')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=1.2*0.5; beta4=1.2*0.08;
epsilonlb=0.49;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,6)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('High transmission, \epsilon_L_B=0.49')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=0.8*0.5; beta4=0.8*0.08;
epsilonlb=0.67;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,7)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Low transmission, \epsilon_L_B=0.67')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=0.5; beta4=0.08;
epsilonlb=0.67;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,8)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('Baseline transmission rate, \epsilon_L_B=0.67')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta3=1.2*0.5; beta4=1.2*0.08;
epsilonlb=0.67;
RC = @(x,y) ((1-q).*[beta4.*(1-x-y)+(1-epsilonlb)*beta4.*x + beta4*0.6*y])./gamrec2...
    + (q.*[beta3.*(1-x-y)+(1-epsilonlb)*beta3.*x+ beta3*y])./(delta2+gamrec2);
z = RC(X,Y);
subplot(3,3,9)
[C,h]=contourf(X,Y,z,[0.2:0.1:2.0]);
clabel(C,h,'FontSize',10,'Color','k')
xlabel('vaccinated')
ylabel('recovered')
title('High transmission, \epsilon_L_B=0.67')

%% Figura 14 Local Sensitivity Analysis for the delta variant 
values=[0.58878 0.41121 -0.9507 -0.04924 0.4859 -0.7053 0.0003912 -0.1099 0.1099];
figure
bh = bar(1:numel(values),diag(values),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9

title('Local Sensitivity Analysis of the Control Reproduction Number of the Delta variant')
set(gca,'XTick', 1:numel(values),'XTickLabel',{'\beta_5','\beta_6','\gamma_3','\delta_3','q','\epsilon_L_B','\epsilon_a','\rho','\alpha'});


%% Global Sensitivity Analysis of the alpha strain vs delta  
clear
clc
%   Intervalos
rangmin=[0.2 0.096 0.16 0.24 0.0396 0.4 0.064 1.4605e-4 1.4605e-3 0.024 0.02 0.696 0.6 0.392 0.0028 0.0243 0.0045 0.0022];
rangmax=[0.3 0.144 0.24 0.36 0.0594 0.6 0.096 2.1907e-4 2.1907e-3 0.036 0.03 0.94 0.9 0.588 0.0043 0.0365 0.0067 0.033];
%   sampling
M=HCV3_lhs(1000,rangmin,rangmax);
PobS=zeros(1000,1);
PobE=zeros(1000,1);
PobI=zeros(1000,1);
PobA=zeros(1000,1);
PobV=zeros(1000,1);
PobEU=zeros(1000,1);
PobIU=zeros(1000,1);
PobAU=zeros(1000,1);
PobD=zeros(1000,1);
PobR=zeros(1000,1);
%   1000 runs
for i=1:1000;
    HCV3_odeSamplingTwoStrain(M(i,1),M(i,2),M(i,3),M(i,4),M(i,5),M(i,6),M(i,7),M(i,8),M(i,9),M(i,10),M(i,11),M(i,12),M(i,13),M(i,14),M(i,15),M(i,16),M(i,17),M(i,18))
    load('qfile.mat');
    PobS(i)=p1;
    PobE(i)=p2;
    PobI(i)=p3;
    PobA(i)=p4;
    PobV(i)=p5;
    PobEU(i)=p6;
    PobIU(i)=p7;
    PobAU(i)=p8;
    PobD(i)=p9;
    PobR(i)=p10;
    end
%%% Ahora DD1 es un vector cuyas entradas son la población de hepatocitos sanos en t=50
%%% y las entradas de DD2 son la población de hepatocitos infectados en t=50
%   Rank transform
for i=1:18; %%%asdasd
    M(:,i)=tiedrank(M(:,i));
end
   PobS=tiedrank(PobS);
    PobE=tiedrank(PobE);
    PobI=tiedrank(PobI);
    PobA=tiedrank(PobA);
    PobV=tiedrank(PobV);
    PobEU=tiedrank(PobEU);
    PobIU=tiedrank(PobIU);
    PobAU=tiedrank(PobAU);
    PobD=tiedrank(PobD);
    PobR=tiedrank(PobR);
%   Empezando las prcc
MS=zeros(1000,21); %%%% asdfsdf +2
ME=zeros(1000,21);
MI=zeros(1000,21);
MA=zeros(1000,21);
MP=zeros(1000,21);
MF=zeros(1000,21);
MEB=zeros(1000,21);
MIB=zeros(1000,21);
MAB=zeros(1000,21);
MR=zeros(1000,21);
MD=zeros(1000,21);
for i=1:18
    MS(:,i)=M(:,i);
    ME(:,i)=M(:,i);
    MI(:,i)=M(:,i);
    MA(:,i)=M(:,i);
    MV(:,i)=M(:,i);
    MEU(:,i)=M(:,i);
    MIU(:,i)=M(:,i);
    MAU(:,i)=M(:,i);
    MD(:,i)=M(:,i);
    MR(:,i)=M(:,i);
end
%   las variables dummy;
MS(:,19)=randperm(1000,1000)';
ME(:,19)=MS(:,19);
MI(:,19)=MS(:,19);
MA(:,19)=MS(:,19);
MV(:,19)=MS(:,19);
MEU(:,19)=MS(:,19);
MIU(:,19)=MS(:,19);
MAU(:,19)=MS(:,19);
MD(:,19)=MS(:,19);
MR(:,19)=MS(:,19);%%%% asdsac +1
%   terminando de crear la matriz para el prcc
MS(:,21)=PobS; %%%asdasdas+2
ME(:,21)=PobE;
MI(:,21)=PobI; %%%asdasdas+2
MA(:,21)=PobA;
MV(:,21)=PobV; %%%asdasdas+2
MEU(:,21)=PobEU; %%%asdasdas+2
MIU(:,21)=PobIU;
MAU(:,21)=PobAU; %%%asdasdas+2
MD(:,21)=PobD;
MR(:,21)=PobR;
%   PRCC
[PCS,PvS]=partialcorr(MS);
[PCE,PvE]=partialcorr(ME);
[PCI,PvI]=partialcorr(MI);
[PCA,PvA]=partialcorr(MA);
[PCV,PvV]=partialcorr(MV);
[PCEU,PvEU]=partialcorr(MEU);
[PCIU,PvIU]=partialcorr(MIU);
[PCAU,PvAU]=partialcorr(MAU);
[PCD,PvD]=partialcorr(MD);
[PCR,PvR]=partialcorr(MR);
%   bar graph
GS=zeros(19,1);
GE=zeros(19,1);
GI=zeros(19,1);
GA=zeros(19,1);
GV=zeros(19,1);
GEU=zeros(19,1);
GIU=zeros(19,1);
GAU=zeros(19,1);
GD=zeros(19,1);
GR=zeros(19,1);
GvS=zeros(19,1);
GvE=zeros(19,1);
GvI=zeros(19,1);
GvA=zeros(19,1);
GvV=zeros(19,1);
GvEU=zeros(19,1);
GvIU=zeros(19,1);
GvAU=zeros(19,1);
GvD=zeros(19,1);
GvR=zeros(19,1);
for i=1:20
    GS(i)=PCS(i,21);
    GE(i)=PCE(i,21);
    GI(i)=PCI(i,21);
    GA(i)=PCA(i,21);
    GV(i)=PCV(i,21);
    GEU(i)=PCEU(i,21);
    GIU(i)=PCIU(i,21);
    GAU(i)=PCAU(i,21);
    GD(i)=PCD(i,21);
    GR(i)=PCR(i,21);
    GvS(i)=PvS(i,21);
    GvE(i)=PvE(i,21);
    GvI(i)=PvI(i,21);
    GvA(i)=PvA(i,21);
    GvV(i)=PvV(i,21);
    GvEU(i)=PvEU(i,21);
    GvIU(i)=PvIU(i,21);
    GvAU(i)=PvAU(i,21);
    GvD(i)=PvD(i,21);
    GvR(i)=PvR(i,21);
end

%%Susceptible
figure
subplot(2,1,1)
bh= bar(1:numel(GS),diag(GS),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis Susceptible Individuals')
set(gca,'XTick', 1:numel(GS),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvS),diag(GvS),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
set(gca,'XTick', 1:numel(GvS),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
title('P values of the parameters')


%% Figure for infected individuals 
figure
subplot(2,1,1)
bh= bar(1:numel(GI),diag(GI),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Infected Symptomatic Individuals to the Wild-type variant')
set(gca,'XTick', 1:numel(GI),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvI),diag(GvI),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvI),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});

%% Figure for Asymptomatic individuals 
figure
subplot(2,1,1)
bh= bar(1:numel(GA),diag(GA),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Infected Asymptomatic Individuals to the Wild-type variant')
set(gca,'XTick', 1:numel(GA),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvA),diag(GvA),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvA),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});

%% Figure for vaccinated individuals 
figure
subplot(2,1,1)
bh= bar(1:numel(GV),diag(GV),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Vaccinated Individuals')
set(gca,'XTick', 1:numel(GV),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvV),diag(GvV),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvV),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});

%% Figure for exposed to the second variant 
figure
subplot(2,1,1)
bh= bar(1:numel(GEU),diag(GEU),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Exposed Individuals to the Alpha Variant')
set(gca,'XTick', 1:numel(GEU),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvEU),diag(GvEU),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvEU),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});

%% Figure for Infected of the second variant 
figure
subplot(2,1,1)
bh= bar(1:numel(GIU),diag(GIU),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Infected Symptomatic Individuals to the Delta Variant')
set(gca,'XTick', 1:numel(GIU),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvIU),diag(GvIU),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvIU),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});

%% Figure for Asymptomatic infections with the second variant 
figure
subplot(2,1,1)
bh= bar(1:numel(GAU),diag(GAU),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Infected Asymptomatic Individuals to the Delta Variant')
set(gca,'XTick', 1:numel(GAU),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvAU),diag(GvAU),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvAU),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});

%% Figure for death individuals 
figure
subplot(2,1,1)
bh= bar(1:numel(GD),diag(GD),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Diseased Individuals')
set(gca,'XTick', 1:numel(GD),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvD),diag(GvD),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvD),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});

%% Recovered individuals 
figure
subplot(2,1,1)
bh= bar(1:numel(GR),diag(GR),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('Global Sensitivity Analysis of Recovered Individuals')
set(gca,'XTick', 1:numel(GR),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});
subplot(2,1,2)
bh= bar(1:numel(GvR),diag(GvR),'stacked','FaceColor', 'c');
bh(1).FaceColor = 'r'; %bar 1 
bh(2).FaceColor = 'b'; %bar 2....
bh(3).FaceColor='g'; %bar 3
bh(4).FaceColor='c'; %bar 4
bh(5).FaceColor='m'; %bar 5
bh(6).FaceColor='y'; %bar 6
bh(7).FaceColor = 'r'; %bar 7 
bh(8).FaceColor = 'b'; %bar 8....
bh(9).FaceColor='g'; %bar 9
bh(10).FaceColor='c'; %bar 10
bh(11).FaceColor='m'; %bar 11
bh(12).FaceColor='y'; %bar 12
bh(13).FaceColor = 'r'; %bar 13 
bh(14).FaceColor = 'b'; %bar 14....
bh(15).FaceColor='g'; %bar 15
bh(16).FaceColor='c'; %bar 16
bh(17).FaceColor='m'; %bar 17
bh(18).FaceColor='y'; %bar 18
title('P values of the parameters')
set(gca,'XTick', 1:numel(GvR),'XTickLabel',{'w','p','q','\beta_1','\beta_2','\beta_3','\beta_4','\delta_1','\delta_2','\gamma_1','\gamma_2','\epsilon_L','\epsilon_L_A','\epsilon_L_B','\epsilon_a','\rho','\alpha','\eta','dummy','dumy'});


