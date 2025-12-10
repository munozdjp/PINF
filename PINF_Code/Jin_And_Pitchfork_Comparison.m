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
% along with this program.  If not, see <https://github.com/munozdjp/HINNDy/blob/main/LICENSE>.

% Structure
% 1) It generates the pitchfork Normalized (to extract 
%   -timePitch scaling and 
%   -maxStatePitch )
%   -maxBetaPitch
% 2) It integrate the jingwang (To find the right scale for the
%   pitchfork)(it finds the maxState of Jin Wang) 
%   - timeJin
%   - maxStateJin       
%    -maxBetaJin
% 3) It normalize pitchfork using the scale of the jin wang and the pitchfork . 
%   - creates finalScale = maxStatepitch/maxStateJIn
% 4) It generates the prediction of jin using pitch equations.
% Equation  = mu = (x * finalscale )^2 / maxBeta

%note: everymodification in 4 needs to be reflected in 2 (Recall this)
%% 1) Generate data for the Pitchfork
clear all, close all, clc
figpath = '../figures/';
addpath('./utils');
addpath('./violin/');
addpath('./NoiseFunctions/');
addpath('./plottingfunctions/');
tic


%generate Data

n = 3; %number of variables
dt=0.01; %timestep

%Integrate interval

Endinterval = 1; %initial value of the pitchfor bifurcation
%Endinterval = 1.5; %initial value of the pitchfor bifurcation
intervalzero = -.001;

%real hidden variable interval
%tinterval=[0 Initinterval];
tinterval=[intervalzero Endinterval];

% Integrate
tspan=[.001:dt:tinterval(2)-tinterval(1)];
N = length(tspan);

Xpitch=[];
%Coeficients of polyorder 3 c1*t+c2*t^2+c3*t^3 & and final time value 
%c=[3,2,1,Endinterval];

cPitch=[1,0,0,Endinterval];
c = cPitch
%Analitic polinomial

%Integration for Pitchforksystem
reescale = 16
%xMax = [1,1,1];

yzNonreescaled=c(1)*(tspan+tinterval(1))+c(2)*(tspan+tinterval(1)).^2 ...
        +c(3).*(tspan+tinterval(1)).^3;

yzero = reescale*(c(1)*(tspan+tinterval(1))+c(2)*(tspan+tinterval(1)).^2 ...
        +c(3).*(tspan+tinterval(1)).^3);

maxBetapitch = max(yzero);
yzero = yzero/max(yzero);

figure(1)
plot(tspan, yzero,'r--')
xlabel('time')
ylabel('beta')
xline(0,'k--');
set(gca,'FontSize',16);
l=legend('Polynomial Pitchfork');
l.FontSize = 14;
l.Location='northeast';
title('polynomial Pitchfork');

xMax = ones((length(c)-1));
        mu0=[yzero(1)];  %initial condition of Beta  
        
        %sqrt(muo)%formula of coordinate x0 in function of beta(mu0)
        x0=real(sqrt(mu0))+0.01; %Initial condition of state variables
        [t,x] = ode45(@(t,x) PitchforkTimeReescale(t,x,c,xMax,reescale,maxBetapitch ),tspan,[tspan(1),mu0,x0]);%,options);
        Xpitch = [Xpitch;x];
        
        %plot variable x vs bifurcation parameter
        figure(2)
        hold on
        plot(x(:,2),x(:,3),'b-')  
        xlabel('beta')
        ylabel('State x3')
        xline(0,'k--');
        yline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('Trajectory with changin beta');
        l.FontSize = 14;
        l.Location='northeast';
        title('Pitchfork State vsparameter')
        hold off    

mufunc = @(x) ((x(:,3).^2))*(xMax(3))^2/maxBetapitch; %here max is one
%mufunc = @(x) reescale*(x(:,3)/xMax(3)).^2
mu_observed=mufunc(x); 

%%Polinomial fit
%In the fututure we can let unkonw the coefficient (tinterval(1)or tinterval(2))

F1=@(weightdx,xdata)  (weightdx(1)*(xdata+tinterval(1))+weightdx(2)*(xdata+tinterval(1)).^2 ...
        +weightdx(3).*(xdata+tinterval(1)).^3);

%curve fitting
weights0=[0.5 0.5 0.5];%initial guess for weights
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, tspan, mu_observed')

%data predicted from F1 evaluated with the predicted coefficients
reproduced_data= F1(weightdx,tspan);
        
plotState_Beta_time_Pred(tspan,Xpitch,mu_observed,yzero,reproduced_data,n,3)

xMaxPitch = max(abs(Xpitch));

%% 2) Obtaining maximum of Jin-Wang scale to reescale pitchfork (Doned)

% Reescale the time between zero and one.

% Do the transition 
%   It is necessary to construct the JingWang Scale Togueter

param = [0.3 0.3 1 1 0.5 4];
a1 = param(1);
a2 = param(2);
k1 = param(3);
k2 = param(4);
s = param(5);
nhill = param(6);

%Number of variables
n = 4;


dt=0.001; %timestep
% Integrate interval
Initinterval = 1;
%neg_inter = 0
neg_inter = -0;
tinterval=[neg_inter Initinterval]; 
tspan=[0:dt:tinterval(2)-tinterval(1)];
N = length(tspan);
%options = odeset('RelTol',1e-50,'AbsTol',1e-50*ones(1,n));
Xwang=[];

%Coeficients of polyorder 3 c1*t+c2*t^2+c3*t^3 & and final time value 
cJin=[ 0.05, 0, 0, Initinterval ]; %2.14,-5.05,3.9143,Initinterval];
c = cJin
%c=[0,0,0,Initinterval];

%Analitic polinomial
reescaleTimeJin = 50;

yzero = reescaleTimeJin*(c(1)*(tspan+tinterval(1))+c(2)*(tspan+tinterval(1)).^2 ...
        +c(3).*(tspan+tinterval(1)).^3);

maxBetaJin= max(yzero);
yzero = yzero/maxBetaJin;

%plot of inverted polinomial
figure(11)
plot(tspan,yzero,'r-')
xlabel('time')
ylabel('beta')
hold on
xline(0,'k--');
set(gca,'FontSize',16);
l=legend('Explicit polinomial for steady state generation','Used polynomial for ODE data generation');
l.FontSize = 14;
l.Location='northeast';
title('Inverted polinomial');

% Integration of hopf system

        mu0=0;   
            %x0=sqrt(yzero(1)/2) %Formula of radius of Initial Condition
        x3= 0.00212578945146621; %Formula of radius of Initial Condition
        x4= 0.00202578945146621;
        [t,x] = ode45(@(t,x) jinWangModel(t,x,c,a1,a2,k1,k2,s,nhill,reescaleTimeJin,maxBetaJin),tspan,[tspan(1),mu0,x3,x4]);%,options);
        Xwang = [Xwang;x];
        %Variable x3 vs bifurcation parameter
        figure(12)
        subplot(2,2,1)
        hold on
        plot(tspan,x(:,3),'b-')
        xlabel('time')
        ylabel('x_3')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('Trajectory using numerical integration');
        l.FontSize = 14;
        l.Location='northeast';
        title('x3 vs b')
        
        subplot(2,2,2)
        plot(x(:,2),x(:,4),'b-')
        xlabel('Beta')
        ylabel('x_4')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('Trajectory x4');
        l.FontSize = 14;
        l.Location='northeast';
        title('x4 vs beta')
        
        subplot(2,2,3)
        plot(x(:,2),x(:,3),'b-')
        xlabel('Beta')
        ylabel('x_3')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('Trajectory x4');
        l.FontSize = 14;
        l.Location='northeast';
        title('x3 vs time')
        
        subplot(2,2,4)
        plot(x(:,1),x(:,4),'b-')
        xlabel('time')
        ylabel('x_4')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('Trajectory x4');
        l.FontSize = 14;
        l.Location='northeast';
        title('x4 vs time')
        hold off    

mufunc = @(x) ((x(:,3).^2))/maxBetaJin;

mu_observed = mufunc(Xwang); 

%%Polinomial fit

%In the fututure we can let unkonw the coefficient (tinterval(1)or tinterval(2))
F1=@(weightdx,xdata)  weightdx(1)*(xdata+tinterval(1))+weightdx(2)*(xdata+tinterval(1)).^2 ...
        +weightdx(3).*(xdata+tinterval(1)).^3;


%curve fitting
weights0=[0.5 0.5 0.5];%initial guess for weights
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, tspan, mu_observed')

%data predicted from F1 evaluated with the predicted coefficients
reproduced_data= F1(weightdx,tspan);



plotState_Beta_time_Pred(tspan,Xwang,mu_observed,yzero,reproduced_data,n,13)

xMaxWang= max(abs(Xwang));
xScale = xMaxPitch./xMaxWang(1:3);
xMaxstateJinPitch = xScale;

%% 3) Reescaling pitchfork witht the maximum of JinWang and also reescaling beta parameter. 

%Verify where is comming from xMaxPitch
%integrate the system of pitchfork with Jing wang scaled: 

%Generate data for the Pitchfork

n = 3; %number of variables
dt=0.01; %timestep

%Integrate interval

Endinterval = 1; %initial value of the pitchfor bifurcation
%Endinterval = 1.5; %initial value of the pitchfor bifurcation
intervalzero = -.001;

%real hidden variable interval
%tinterval=[0 Initinterval];
tinterval=[intervalzero Endinterval];

% Integrate
tspan2=[tinterval(2):-dt:tinterval(1)]; %interval for solving the 
tspan=[.001:dt:tinterval(2)-tinterval(1)];
N = length(tspan);

%options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
Xpitch=[];

%Coeficients of polyorder 3 c1*t+c2*t^2+c3*t^3 & and final time value 
%c=[3,2,1,Endinterval];

c=cPitch;

reescaletime = 16;

yzNonreescaled=c(1)*(tspan+tinterval(1))+c(2)*(tspan+tinterval(1)).^2 ...
        +c(3).*(tspan+tinterval(1)).^3;

yzero = reescaletime*(c(1)*(tspan+tinterval(1))+c(2)*(tspan+tinterval(1)).^2 ...
        +c(3).*(tspan+tinterval(1)).^3);
maxBetapitch = max(abs(yzero));
yzero = yzero/maxBetapitch
figure(21)
plot(tspan, yzero,'r--')
xlabel('time')
ylabel('beta')
xline(0,'k--');
set(gca,'FontSize',16);
l=legend('Polynomial Pitchfork');
l.FontSize = 14;
l.Location='northeast';
title('polynomial Pitchfork');


        mu0=[yzero(1)];  %initial condition of Beta  
        
        %sqrt(muo)%formula of coordinate x0 in function of beta(mu0)
        x0=real(sqrt(mu0))+0.01; %Initial condition of state variables
        [t,x] = ode45(@(t,x) PitchforkTimeReescale(t,x,c,xMaxstateJinPitch,reescaletime,maxBetapitch),tspan,[tspan(1),mu0,x0]);%,options);
        Xpitch = [Xpitch;x];
        
        %plot variable x vs bifurcation parameter
        figure(22)
        hold on
        plot(x(:,2),x(:,3),'b-')  
        xlabel('beta')
        ylabel('State x3')
        xline(0,'k--');
        yline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('Trajectory with changin beta');
        l.FontSize = 14;
        l.Location='northeast';
        title('Pitchfork State vsparameter')
        hold off    

mufunc = @(x) ((x(:,3).^2))*xMaxstateJinPitch(3)^2/maxBetapitch
mu_observed=mufunc(x); 

%%Polinomial fit

%In the fututure we can let unkonw the coefficient (tinterval(1)or tinterval(2))
F1=@(weightdx,xdata)  weightdx(1)*(xdata+tinterval(1))+weightdx(2)*(xdata+tinterval(1)).^2 ...
        +weightdx(3).*(xdata+tinterval(1)).^3;

%curve fitting
weights0=[0.5 0.5 0.5];%initial guess for weights
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, tspan, mu_observed')

reproduced_data= F1(weightdx,tspan);

plotState_Beta_time_Pred(tspan, Xpitch, mu_observed, yzero, reproduced_data, n, 23)

%% 4) Prediction of dynamics Jin wang using pitchfork equations. 

param = [0.3 0.3 1 1 0.5 4];
a1 = param(1);
a2 = param(2);
k1 = param(3);
k2 = param(4);
s = param(5);
nhill = param(6);

%Number of variables
n = 4;


dt=0.001; %timestep
% Integrate interval
Initinterval = 1;
%neg_inter = 0
neg_inter = -0;
tinterval=[neg_inter Initinterval]; 
tspan=[0:dt:tinterval(2)-tinterval(1)];
N = length(tspan);
%options = odeset('RelTol',1e-50,'AbsTol',1e-50*ones(1,n));
Xwang=[];

%Coeficients of polyorder 3 c1*t+c2*t^2+c3*t^3 & and final time value 
%c=[1/50,0,0,Initinterval];
%c=[2.14,-5.05,3.9143,Initinterval];
c=cJin; %2.14,-5.05,3.9143,Initinterval];
%c=[0,0,0,Initinterval];

%Analitic polinomial
reescaleJin = 50;

yzero = reescaleJin*(c(1)*(tspan+tinterval(1))+c(2)*(tspan+tinterval(1)).^2 ...
        +c(3).*(tspan+tinterval(1)).^3);

maxBetaJin= max(yzero);
yzero = yzero/maxBetaJin;

%plot of inverted polinomial
figure(31)
plot(tspan,yzero,'r-')
xlabel('time')
ylabel('beta')
hold on
xline(0,'k--');
set(gca,'FontSize',16);
l=legend('Explicit polinomial for steady state generation','Used polynomial for ODE data generation');
l.FontSize = 14;
l.Location='northeast';
title('Inverted polinomial');

% Integration of hopf system

        mu0=0;   
        %for mu0=0.0;   
            %x0=sqrt(yzero(1)/2) %Formula of radius of Initial Condition
        x3= 0.00212578945146621; %Formula of radius of Initial Condition
        x4= 0.00202578945146621;
        %for x3= 0 % Formula of radius of Initial Condition
        %    x4= 0
        [t,x] = ode45(@(t,x) jinWangModel(t,x,c,a1,a2,k1,k2,s,nhill,reescaleJin,maxBetaJin),tspan,[tspan(1),mu0,x3,x4]);%,options);
        Xwang = [Xwang;x];
        %Variable x3 vs bifurcation parameter
        figure(32)
        subplot(2,2,1)
        hold on
        plot(tspan,x(:,3),'b-')
        xlabel('time')
        ylabel('x_3')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('Trajectory using numerical integration');
        l.FontSize = 14;
        l.Location='northeast';
        title('x3 vs b')
        
        subplot(2,2,2)
        plot(x(:,2),x(:,4),'b-')
        xlabel('Beta')
        ylabel('x_4')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('Trajectory x4');
        l.FontSize = 14;
        l.Location='northeast';
        title('x4 vs beta')
        
        subplot(2,2,3)
        plot(x(:,2),x(:,3),'b-')
        xlabel('Beta')
        ylabel('x_3')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('Trajectory x4');
        l.FontSize = 14;
        l.Location='northeast';
        title('x3 vs time')
        
        subplot(2,2,4)
        plot(x(:,1),x(:,4),'b-')
        xlabel('time')
        ylabel('x_4')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('Trajectory x4');
        l.FontSize = 14;
        l.Location='northeast';
        title('x4 vs time')
        hold off    

mufunc = @(x) ((x(:,3).^2))*xMaxstateJinPitch(3)^2/maxBetapitch;

mu_observed = mufunc(Xwang); 

%%Polinomial fit

%In the fututure we can let unkonw the coefficient (tinterval(1)or tinterval(2))
F1=@(weightdx,xdata)  weightdx(1)*(xdata+tinterval(1))+weightdx(2)*(xdata+tinterval(1)).^2 ...
        +weightdx(3).*(xdata+tinterval(1)).^3;


%curve fitting
weights0=[0.5 0.5 0.5];%initial guess for weights
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, tspan, mu_observed')

%data predicted from F1 evaluated with the predicted coefficients
reproduced_data= F1(weightdx,tspan);



plotState_Beta_time_Pred(tspan,Xwang,mu_observed,yzero,reproduced_data,n,33)
