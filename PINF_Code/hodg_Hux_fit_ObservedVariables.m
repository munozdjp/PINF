% Copyright (C) 2023 (King Abdullah University of Science and Technology)
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://github.com/munozdjp/IHCV/blob/main/LICENSE>.

clear all, clc, clear, close all

%===simulation time===
simulationTime = 100; %in milliseconds
deltaT=.01;
t=0:deltaT:simulationTime;


%===specify the external current I===
changeTimes = [0]; %in milliseconds
currentLevels = [50]; %Change this to see effect of different currents on voltage (Suggested values: 3, 20, 50, 1000)

%polinomial structure for the bifurcation parameter
p1 =         1.1;
p2 =      0.9285;
p3 =    -0.02771;
p4 =   0.0003581;

Initinterval = 250;
%neg_inter = 0;
%tinterval=[neg_inter Initinterval];
c=[p1,p2,p3,p4,Initinterval];
y2=c(1)+c(2)*(t)+c(3)*(t).^2+c(4)*(t).^3;
%Set externally applied current across time

I=y2;

%===constant parameters===%
%All of these can be found in Table 3
gbar_K=36; gbar_Na=120; g_L=.3;
E_K = -12; E_Na=115; E_L=10.6;
C=1;


%===set the initial states===%
V=0; %Baseline voltage
alpha_n = .01 * ( (10-V) / (exp((10-V)/10)-1) ); %Equation 12
beta_n = .125*exp(-V/80); %Equation 13
alpha_m = .1*( (25-V) / (exp((25-V)/10)-1) ); %Equation 20
beta_m = 4*exp(-V/18); %Equation 21
alpha_h = .07*exp(-V/20); %Equation 23
beta_h = 1/(exp((30-V)/10)+1); %Equation 24

n(1) = alpha_n/(alpha_n+beta_n); %Equation 9
m(1) = alpha_m/(alpha_m+beta_m); %Equation 18
h(1) = alpha_h/(alpha_h+beta_h); %Equation 18


for i=1:numel(t)-1 %Compute coefficients, currents, and derivates at each time step
   
    %---calculate the coefficients---%
    %Equations here are same as above, just calculating at each time step
    alpha_n(i) = .01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) );
    beta_n(i) = .125*exp(-V(i)/80);
    alpha_m(i) = .1*( (25-V(i)) / (exp((25-V(i))/10)-1) );
    beta_m(i) = 4*exp(-V(i)/18);
    alpha_h(i) = .07*exp(-V(i)/20);
    beta_h(i) = 1/(exp((30-V(i))/10)+1);
   
   
    %---calculate the currents---%
    I_Na = (m(i)^3) * gbar_Na * h(i) * (V(i)-E_Na); %Equations 3 and 14
    I_K = (n(i)^4) * gbar_K * (V(i)-E_K); %Equations 4 and 6
    I_L = g_L *(V(i)-E_L); %Equation 5
    %I_ion = I(i) - I_K - I_Na - I_L;
    I_ion = I(i) - I_K - I_Na - I_L;
   
   
    %---calculate the derivatives using Euler first order approximation---%
    V(i+1) = V(i) + deltaT*I_ion/C;
    n(i+1) = n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i)); %Equation 7
    m(i+1) = m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i)); %Equation 15
    h(i+1) = h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i)); %Equation 16

end

%%Construction of the hidden variable
% Hidden variable. 
mu_observed= (m.^3).* gbar_Na.* h.* (V-E_Na)...
    + (n.^4).* gbar_K.* (V-E_K)...
    + g_L.*(V-E_L);


%%Polinomial fit
%Define the function to fi: "Poly order 2"
F1=@(weightdx,time) weightdx(1)+weightdx(2)*time+time.^2*weightdx(3)+...
    time.^3*weightdx(4);
%Curve fitting
weights0=[0.5 0.5 0.5 0.5];
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, t, mu_observed)

%Data predicted from F1 evaluated in the polinomail fit
reproduced_data= F1(weightdx,t);

% 
 V = V-70; %Set resting potential to -70mv
% 
% 
%===plot Voltage===%
figure(1)
plot(t,V,'LineWidth',3)
hold on
l=legend({'voltage'})
l.FontSize = 14;
l.Location='northeast';

ylabel('Voltage (mv)')
xlabel('time (ms)')
set(gca,'FontSize',16);
title('Voltage over Time in Simulated Neuron')
hold off

figure(2)
plot(t,y2,'r-','LineWidth',1)% good real dynamic
yline(0,'k--','HandleVisibility','off');
hold on
plot(t,reproduced_data,'b--','LineWidth',1); %Beta from polinomial fit
plot(t,mu_observed,'m-')%,'LineWidth',1) %Beta from steady equations
xlabel('Time')
ylabel('Beta')
%axis([0 20 -100 5000])
set(gca,'FontSize',16);
l=legend('Real hidden variable:analitical polinomial ','Predicted hidden variable:With Polinomial weights','hidden observed from ODE data(x observations) ');
l.FontSize = 14;
l.Location='northeast';
title('Beta from: real dynamic, polinomial fit, Observed from steady equations')
hold off
%Other plots.


x = [t',y2',V'];

plotState_Beta_time_Pred(t,x,mu_observed,y2,reproduced_data,3,3)
%%
figure(nfig)
        subplot(2,2,1)
        hold on
        plot(t,x(:,3),'b-','linewidth',2)
        xlabel('time')
        ylabel('x')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('X VS Time');
        l.FontSize = 14;
        l.Location='northeast';
        title('Observation State')
        hold off
        
        subplot(2,3,2)
        plot(t,yzero,'r-','linewidth',2)% good True dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(t,reproduced_data,'b--','linewidth',2); %Alpha from polinomial fit
        xlabel('Time')
        ylabel('Alpha')
        %axis([0 20 -100 5000])
        set(gca,'FontSize',16);
        l=legend('True ','Inferred control parameter ','Intermediate step ');
        l.FontSize = 14;
        l.Location='northeast';
        title('True, Inferred ')

        hold off                
        
        subplot(2,3,3)
        plot(t,yzero,'r-','linewidth',2)% good True dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(t,reproduced_data,'b--','linewidth',2); %Alpha from polinomial fit
        plot(t,mu_observed,'m-','linewidth',1)%,'linewidth',2) %Alpha from steady equations
        xlabel('Time')
        ylabel('Alpha')
        %axis([0 20 -100 5000])
        set(gca,'FontSize',16);
        l=legend('True ','Inferred control parameter ','Intermediate step ');
        l.FontSize = 14;
        l.Location='northeast';
        title('True, Inferred, intermediate')
        hold off
        
        subplot(2,3,4)
        hold on
        plot(x(:,2),x(:,3),'b-','linewidth',2)
        xlabel('Alpha')
        ylabel('x')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend(' x VS Time');
        l.FontSize = 14;
        l.Location='northeast';
        title('State vs Alpha')
        hold off
