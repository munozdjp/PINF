clear all, close all, clc
figpath = '../figures/';
addpath('./utils');
addpath('./violin/');
addpath('./NoiseFunctions/');
addpath('./plottingfunctions/');
rng(1)
%%
% Parameter to specify the type of analysis
analysisType = 'pitchfork';  % Change this to 'pitchfork' as needed when running different scripts

% Base directory for results
scriptFullPath = mfilename('fullpath');
[scriptPath, ~, ~] = fileparts(scriptFullPath);
resultsBasePath = fullfile(scriptPath, 'Results figures');

% Dynamic directory based on the analysis type
analysisPath = fullfile(resultsBasePath, analysisType);

% Ensure the directory exists
if ~exist(analysisPath, 'dir')
    mkdir(analysisPath);
end

% Define the path for the PDF file within the analysis directory
pdfPath = fullfile(analysisPath, 'ResultsFigures.pdf');
tic

%generate Data

n = 3;%number of variables
dt=0.01;%timestep
%Integrate interval

Endinterval = 16; %initial value of the pitchfor bifurcation
%Endinterval = 1.5; %initial value of the pitchfor bifurcation
intervalzero = -.001;
%real hidden variable interval
%tinterval=[0 Initinterval];
tinterval=[intervalzero Endinterval];




% Integrate
tspan2=[tinterval(2):-dt:tinterval(1)]; %interval for solving the 
tspan=[.001:dt:tinterval(2)-tinterval(1)];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
X=[];
%Coeficients of polyorder 3 c1*t+c2*t^2+c3*t^3 & and final time value 
%c=[3,2,1,Endinterval];
c=[1,0,0,Endinterval];
%Analitic polinomial
    y=c(1)*(tspan2)+c(2)*(tspan2).^2+c(3)*(tspan2).^3;
    %inverted polinomial
    y2=c(1)*(-tspan+tinterval(2) )+c(2)*(-tspan...
        +tinterval(2)).^2+c(3)*(-tspan+tinterval(2)).^3;
    %plot(tspan,y,'r-',tspan2,y2,'b--')
    yzero = c(1)*(tspan+tinterval(1))+c(2)*(tspan+tinterval(1)).^2 ...
        +c(3).*(tspan+tinterval(1)).^3;
    %Derivative of inverted polinomilal is 

figure(1)
plot(tspan2,y,'r-',tspan,y2,'b--',tspan,yzero,'k--')
xlabel('Time', 'FontName', 'Times New Roman')  % Set xlabel font to Times New Roman

xlim([min(tspan) max(tspan)]);
ylabel('\alpha')  % Use \alpha for the Greek letter alpha
ylim([min(y) max(y)]);

% xline(0,'k--');
set(gca,'FontSize',16);
l=legend('Explicit polynomial steady-state generation', 'Implicit polynomial - ODE data generation', 'Polynomial moved - positive axis');
l.Position = [0.46, 0.83, 0.07, 0.071];

l.FontSize = 12;
% l.Location='northeast';
title('Inverted polynomial');
exportgraphics(gcf, pdfPath, 'ContentType', 'vector', 'Append', true);

%Integration for Pitchforksystem
xMax = [1,1,1];
for mu0=[yzero(1)];  %initial condition of Beta  
    %sqrt(muo)%formula of coordinate x0 in function of beta(mu0)
    for x0=real(sqrt(mu0))+0.01 %Initial condition of state variables
        [t,x] = ode45(@(t,x) pitchforkPolyorder3Left2Right(t,x,c,xMax),tspan,[tspan(1),mu0,x0],options);
         X = [X;x];
    end
end

xMax = max(abs(x))
% load("MaxJin.mat")
% jingWangScale = xMax./MaxJin(1:3)
% xMax = jingWangScale

for mu0=[yzero(1)];  %initial condition of Beta  
    %sqrt(muo)%formula of coordinate x0 in function of beta(mu0)
    for x0=real(sqrt(mu0))+0.01 %Initial condition of state variables
        [t,x] = ode45(@(t,x) pitchforkPolyorder3Left2Right(t,x,c,xMax),tspan,[tspan(1),mu0,x0],options);
         X = [X;x];
         %plot variable x vs bifurcation parameter
            figure(2)
            hold on
            plot(x(:,2),x(:,3),'b-')    
            xlabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font for the xlabel
            ylabel('State', 'FontName', 'Times New Roman')  % Times New Roman font for the ylabel
            xline(0,'k--');
            set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');  % Set the axes font to Times New Roman

            xline(0,'k--');
            yline(0,'k--');
            set(gca,'FontSize',16,'FontName','Times New Roman');
            % l=legend('Trajectory using numerical integration');
            % l.FontSize = 14;
            % l.Location='northeast';
            title('State vs \alpha (t)', 'FontName', 'Times New Roman')  % Use LaTeX markup for alpha and set font to Times New Roman
            hold off
            exportgraphics(gcf, pdfPath, 'ContentType', 'vector', 'Append', true);
    end
end

%Mu observed from equations of steady state
%[x,yzero] = xNormalization(x,yzero);

mufunc = @(x) ((x(:,3).^2))*xMax(3)^2
mu_observed=mufunc(x); 
%%Polinomial fit

%In the fututure we can let unkonw the coefficient (tinterval(1)or tinterval(2))
F1=@(weightdx,xdata)  weightdx(1)*(xdata+tinterval(1))+weightdx(2)*(xdata+tinterval(1)).^2 ...
        +weightdx(3).*(xdata+tinterval(1)).^3;


%curve fitting
weights0=[0.5 0.5 0.5];%initial guess for weights
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, tspan, mu_observed')

%data predicted from F1 evaluated with the predicted coefficients
reproduced_data= F1(weightdx,tspan);


plotState_Beta_time_Pred(tspan,x,mu_observed,yzero,reproduced_data,n,3,pdfPath)

%%Data with noise: 

noise_scale = 0.15;

[xNoisy, Noise_normalized] = add_Noise_Max(x,noise_scale);

%next step is to use the noise for the figures. 
mu_obsNoisy = mufunc(xNoisy);

[weightdxNoise,resnormNoise,~,exitflagNoise,outputNoise] = lsqcurvefit(F1,weights0, tspan, mu_obsNoisy');

reproduced_dataNoisy= F1(weightdxNoise,tspan);

comparisonVector = [c(1:(length(c)-1));weightdx;weightdxNoise]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotState_Beta_time_Pred(tspan,xNoisy,mu_obsNoisy,yzero,reproduced_dataNoisy,n,4,pdfPath)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotCleanNoiseStateFig_6(tspan,Noise_normalized,yzero,n)


%% Noise analysis 
% Creation of comparison vector for different noise variances: 
VectorOfNoise = [0:0.1:0.5];
name = mfilename
name = name(1:9)
%noise_Scale_maximum(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc)
figure
noise_5_Boxplots(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc,pdfPath)

