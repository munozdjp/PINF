clear all, close all, clc
figpath = '../figures/';
addpath('./utils');
addpath('./NoiseFunctions/');
addpath('./plottingfunctions/');

%%
rng(1)
% Parameter to specify the type of analysis
analysisType = 'Lorentz';  % Change this to 'pitchfork' as needed when running different scripts

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

%coefficient of the right polinomial (polinomial of order 4 ...
%that perform the desired behavior (increasing))
%%
p1 =         1.1;
p2 =      0.9285;
p3 =    -0.02771;
p4 =   0.0003581;

% p1 =         1.1;
% p2 =      1.9285;
% p3 =    2.02771;
% p4 =   0.43581;

omega= 10;
beta=2.66;
rho=1.1;
n = 5; % Number of variables
dt=0.01; % Timestep
%dt=1; % Timestep

% Integrate Interval
Initinterval = 250;
%Initinterval = 160;
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
xlabel('Time', 'FontName', 'Times New Roman')  % Set xlabel font to Times New Roman
ylabel('\alpha')  % Use \alpha for the Greek letter alpha
hold on
xline(0,'k--');
set(gca,'FontSize',16);
l = legend('Explicit polynomial steady-state generation')
l.FontSize = 12;
l.Location='northwest';
title('Inverted polinomial');
hold off
exportgraphics(gcf, pdfPath, 'ContentType', 'vector', 'Append', true);

traslation=0; %traslation of the system in y coordinate
% integration of the lorentz system
%     i=i+1
for mu0=[c(1)]; %initial condition of Beta
    for x0=0.51; %Iniitial contions for x and y [x0,x0,z0]        
        y0=0.51;
        z0=0.1;  %initial condition for z0
        [t,x] = ode45(@(t,x) lorentz2(t,x,c,omega,beta,traslation),tspan,[tspan(1),mu0,x0,y0,z0],options);
         X = [X;x];
        % Plot of state variable vs Beta parameter
        figure(2)
        hold on
        plot(x(:,2),x(:,3),'b-')
        xlabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font for the xlabel
        ylabel('State', 'FontName', 'Times New Roman')  % Times New Roman font for the ylabel
        xline(0,'k--');
        set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');  % Set the axes font to Times New Roman
        xlim([min(x(:,2)) max(x(:,2))]);
        title('State vs \alpha (t)', 'FontName', 'Times New Roman')  % Use LaTeX markup for alpha and set font to Times New Roman
        hold off
        exportgraphics(gcf, pdfPath, 'ContentType', 'vector', 'Append', true);
    end
end

%Normalization
xMax = max(abs(x));
[x,yzero] = xNormalization(x,yzero);

mufunc = @(x) ((x(:,3)*xMax(3)).*(x(:,5)*xMax(5))+(x(:,4)*xMax(4)))./(x(:,3)*xMax(3));
mu_observed = mufunc(x);


% Polinomial fit
%Define the function to fi: "Poly order 2"
F1=@(weightdx,time) weightdx(1)+weightdx(2)*time+time.^2*weightdx(3)+...
    time.^3*weightdx(4);
%Curve fitting
weights0=[0.5 0.5 0.5 0.5];
opts = optimoptions('lsqcurvefit','Display','off');
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, tspan, mu_observed',[],[],opts);

%Data predicted from F1 evaluated in the polinomail fit
reproduced_data= F1(weightdx,tspan);

plotState_Beta_time_Pred(tspan,x,mu_observed,yzero,reproduced_data,4,3,pdfPath)
%% Mu with ratio State-Noise

noise_scale = 0.15;

[xNoisy, Noise_normalized] = add_Noise_Max(x,noise_scale);

%next step is to use the noise for the figures. 
mu_obsNoisy = mufunc(xNoisy);

[weightdxNoise,resnormNoise,~,exitflagNoise,outputNoise] = lsqcurvefit(F1,weights0, tspan, mu_obsNoisy');

reproduced_dataNoisy= F1(weightdxNoise,tspan);

comparisonVector = [c(1:(length(c)-1));weightdx;weightdxNoise]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotNoisyStateFig_4(tspan,xNoisy,mu_obsNoisy,yzero,reproduced_dataNoisy,n)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% plotCleanNoiseStateFig_5(tspan,Noise_normalized,yzero,n)
%% Noise analysis 
% Creation of comparison vector for different noise variances: 
VectorOfNoise = [0:0.1:0.5];
name = mfilename;
name = name(1:10);
figure
%noise_Scale_maximum(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc)
noise_5_Boxplots(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc,pdfPath)

