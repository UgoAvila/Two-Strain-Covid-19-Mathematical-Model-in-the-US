function HCV3_odeSampling(w, p,q, beta1, beta2,beta3,beta4, delta1,delta2, gamrec1,gamrec2,epsilonl,epsilonla,epsilonlb,epsilona,rho,alpha,eta)

t=0:1:100; %time scale

initial_S=331449281;
initial_E0 = 15000;
initial_I0 = 15000;
initial_A0 = ((1-0.2)/0.2)*initial_I0;
initial_V0=140000000;
initial_EU0 = 1000;
initial_IU0 = 1000;
initial_AU0 = ((1-0.2)/0.2)*initial_IU0;
initial_R0 = 1.1125e+06;
initial_D0 = 62100;

[t,x]=ode45(@rhs,t,[initial_S initial_E0 initial_I0 initial_A0 initial_V0 initial_EU0 initial_IU0 initial_AU0 initial_R0 initial_D0]);

D1=x(:,1);
D2=x(:,2);
D3=x(:,3);
D4=x(:,4);
D5=x(:,5);
D6=x(:,6);
D7=x(:,7);
D8=x(:,8);
D9=x(:,9);
D10=x(:,10);


D1=D1(end);
D2=D2(end);
D3=D3(end);
D4=D4(end);
D5=D5(end);
D6=D6(end);
D7=D7(end);
D8=D8(end);
D9=D9(end);
D10=D10(end);

p1=D1;
p2=D2;
p3=D3;
p4=D4;
p5=D5;
p6=D6;
p7=D7;
p8=D8;
p9=D9;
p10=D10;


save('qfile','p1','p2','p3','p4','p5','p6','p7','p8','p9','p10')
plot(t,x(:,1))
hold on
%   -x(1)*x(2)*delta2
    function dydt=rhs(t,y)
       % Susceptibles
dS_dt = - beta1*y(1).*y(3)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)) ...
        - beta2*y(1).*y(4)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)) -beta3*y(1).*y(7)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)) ...
        - beta4*y(1).*y(8)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)) + alpha*y(5) ...
    -(1-epsilona)*rho*y(1) + eta*y(10);

% Exposed
dE_dt = beta1*y(1).*y(3)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)) ...
        + beta2*y(1).*y(4)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10))... 
        + (1-epsilonl)*beta1*y(3)*y(5)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)) ... 
        + (1-epsilonla)*beta2*y(3)*y(5)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10))  -w*y(2);  

% Infected
dI_dt = p*w*y(2) - delta1*y(3) - gamrec1*y(3);

% Infected Asymptomatic
dA_dt = (1-p)*w*y(2) - gamrec1*y(4); 

% Immunity provided by two dose 
dV_dt = (1-epsilona)*rho*y(1)...
          - (1-epsilonl)*beta1*y(3)*y(5)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)) ... 
      - (1-epsilonla)*beta2*y(4)*y(5)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)) ...
    - alpha*y(5) ...
    -(1-epsilonlb)*beta3*y(5).*y(7)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)) ...
    -(1-epsilonlb)*beta4*y(5).*y(8)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)); 
%Exposed to the other variant 
dEU_dt = + beta3*y(1).*y(7)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)) ...
        + beta4*y(1).*y(8)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)) ... 
        +(1-epsilonlb)*beta3*y(5).*y(7)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)) ...
    +(1-epsilonlb)*beta4*y(5).*y(8)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10))...
    +beta3*y(10).*y(7)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10))...
    +beta4*y(10).*y(8)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10))- w*y(6);
% Infected with the second variant cases 
dIU_dt = q*w*y(6) - delta2*y(7) - gamrec2*y(7);
%Infected but asymptomatic breakthrough cases 
dAU_dt = (1-q)*w*y(6) - gamrec2*y(8); 
%death cases despite their nature
dD_dt = delta1*y(3) + delta2*y(7); 
%Death by breakthrough or not
dR_dt= gamrec1*(y(3)+y(4))+ gamrec2*(y(7)+y(8))-beta3*y(10).*y(7)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10))...
    -beta4*y(10).*y(8)./(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10))-eta*y(10);

dydt = [dS_dt; dE_dt; dI_dt; dA_dt; dV_dt; dEU_dt; dIU_dt; dAU_dt; dD_dt; dR_dt];
    end
end