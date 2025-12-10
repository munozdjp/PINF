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

clear all, close all, clc
figpath = '../figures/';
addpath('./utils');
addpath('./NoiseFunctions/');
addpath('./plottingfunctions/');

%coefficient of the right polinomial (polinomial of order 4 ...
%that perform the desired behavior (increasing))

p1 =         1.1;
p2 =      0.9285;
p3 =    -0.02771;
p4 =   0.0003581;

n = 4; % Number of variables
dt=0.01; % Timestep

% Integrate Interval
Initinterval = 250;
neg_inter = 0;
tinterval=[neg_inter Initinterval];
tspan=[.001:dt:tinterval(2)-tinterval(1)];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
X=[];

%Coeficients of polyorder 3; p1+p2*t+p3*t^2+p4*t^3 & and initial time value 
c=[p1,p2,p3,p4,Initinterval];
yzero=c(1)+c(2)*(tspan)+c(3)*(tspan).^2+c(4)*(tspan).^3;

%Plot of inverted polynomial
figure(1)
plot(tspan,yzero,'r-')
xlabel('time')
ylabel('beta parameter')
hold on
xline(0,'k--');
set(gca,'FontSize',16);
l=legend('Explicit polinomial for steady state generation');
l.FontSize = 14;
l.Location='northeast';
title('Inverted polinomial');
hold off

traslation=0; %traslation of the system in y coordinate
%integration of the lorentz system
for mu0=[c(1)]; %initial condition of Beta
    for x0=0.51; %Iniitial contions for x and y [x0,y0]        
        y0=0.51;
        [t,x] = ode45(@(t,x) vanderpol(t,x,c),tspan,...
            [tspan(1),mu0,x0,y0],options);
                X = [X;x];
            % Plot of state variable vs Beta parameter
            figure(2)
            hold on
            plot(x(:,2),x(:,3),'b-')
            xlabel('Beta')
            ylabel('x_3')
            xline(0,'k--');
            set(gca,'FontSize',16);
            l=legend('Trajectory using numerical integration');
            l.FontSize = 14;
            l.Location='northeast';
            title('State variable vs Biffurcation parameter')
            hold off    
    end
end

%Normalization
xMax = max(abs(x));
[x,yzero] = xNormalization(x,yzero);

%Mu observed from equations of steady state
mufunc = @(x) (x(:,3)*xMax(3))./((1-(x(:,3)*xMax(3)).^2).*(x(:,4)*xMax(4)));
mu_observed= mufunc(x);

%%Polinomial fit
%Define the function to fi: "Poly order 2"
F1=@(weightdx,time) weightdx(1)+weightdx(2)*time+time.^2*weightdx(3)+...
    time.^3*weightdx(4);
%Curve fitting
weights0=[0.5 0.5 0.5 0.5];
opts = optimoptions('lsqcurvefit','Display','off');
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, tspan, mu_observed',[],[],opts);

%Data predicted from F1 evaluated in the polinomail fit
reproduced_data= F1(weightdx,tspan);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot of State Variable vs Timespan
plotState_Beta_time_Pred(tspan,x,mu_observed,yzero,reproduced_data,n,4)

%% Mu with ratio State-Noise

noise_scale = 0.0001;

[xNoisy, Noise_normalized] = add_Noise_Max(x,noise_scale);

%next step is to use the noise for the figures. 
mu_obsNoisy = mufunc(xNoisy);

[weightdxNoise,resnormNoise,~,exitflagNoise,outputNoise] = lsqcurvefit(F1,weights0, tspan, mu_obsNoisy');

reproduced_dataNoisy= F1(weightdxNoise,tspan);

comparisonVector = [c(1:(length(c)-1));weightdx;weightdxNoise]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotNoisyStateFig_5(tspan,xNoisy,mu_obsNoisy,yzero,reproduced_dataNoisy,n)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotCleanNoiseStateFig_6(tspan,Noise_normalized,yzero,n)

%% Noise analysis 
% Creation of comparison vector for different noise variances: 
VectorOfNoise = [0:0.008:0.01];
name = mfilename;
name = name(1:10);
%noise_Scale_maximum(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc)
noise_5_Boxplots(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc)
 
