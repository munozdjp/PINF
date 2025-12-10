clear all, close all, clc
figpath = '../figures/';
addpath('./utils');
addpath('./NoiseFunctions/');
addpath('./plottingfunctions/');

% 
load('variablesFitzHug.mat')
%%
rng(1)
% Parameter to specify the type of analysis
analysisType = 'FitzHugh';  % Change this to 'pitchfork' as needed when running different scripts

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

%% 
figure(1)
plot(tspan,yzero,'r-')
xlabel('Time', 'FontName', 'Times New Roman')  % Set xlabel font to Times New Roman
ylabel('\alpha')  % Use \alpha for the Greek letter alpha
hold on
xline(0,'k--');
set(gca,'FontSize',16);
l = legend('Explicit polynomial steady-state generation');
l.Position = [0.5, 0.84, 0.07, 0.071];
l.FontSize = 12;
% l.Location='northeast';
title('Inverted polynomial');
exportgraphics(gcf, pdfPath, 'ContentType', 'vector', 'Append', true);

% Plot of state variable vs Beta parameter
figure(2)
hold on
plot(x(:,2),x(:,3),'b-')
xlabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font for the xlabel
ylabel('State', 'FontName', 'Times New Roman')  % Times New Roman font for the ylabel
xline(0,'k--');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');  % Set the axes font to Times New Roman
% l = legend('state variable vs bifurcation parameter');
% l.FontSize = 14;
% l.FontName = 'Times New Roman';  % Set the legend font to Times New Roman
% l.Location = 'northeast';
title('State vs \alpha (t)', 'FontName', 'Times New Roman')  % Use LaTeX markup for alpha and set font to Times New Roman
hold off
exportgraphics(gcf, pdfPath, 'ContentType', 'vector', 'Append', true);

%%
xMax = max(abs(x));
[x,yzero] = xNormalization(x,yzero);

mufunc = @(x) ((x(:,4)*xMax(4))-epsilon*((x(:,3)*xMax(3)).*((x(:,3)*xMax(3))-lambda).*(1-(x(:,3)*xMax(3)))))/maxBeta;
mu_observed=mufunc(x);
            
%%Polinomial fit
%Define the function to fi: "Poly order 2"
F1=@(weightdx,time) (weightdx(1)+weightdx(2)*time+time.^2*weightdx(3)+...
    time.^3*weightdx(4))/maxBeta;
%Curve fitting
weights0=[0.5 0.5 0.5 0.5];
opts = optimoptions('lsqcurvefit','Display','off');
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, tspan, mu_observed',[],[],opts);

%Data predicted from F1 evaluated in the polinomail fit
reproduced_data= F1(weightdx,tspan);
%here we are working with all the variables normalized. 
polyorder1 = 3
plotState_Beta_time_Pred(tspan,x,mu_observed,yzero,reproduced_data,n,3,pdfPath)

%% Mu with Noise experimental single noise
noise_scale=0.5;   

[xNoisy, Noise_normalized] = add_Noise_Max(x,noise_scale);

%next step is to use the noise for the figures. 
mu_obsNoisy = mufunc(xNoisy);
opts = optimoptions('lsqcurvefit','Display','off');
[weightdxNoise,resnormNoise,~,exitflagNoise,outputNoise] = lsqcurvefit(F1,weights0, tspan, mu_obsNoisy',[],[],opts);

reproduced_dataNoisy= F1(weightdxNoise,tspan);

comparisonVector = [c(1:(length(c)-1));weightdx;weightdxNoise]



%% Noise analysis 
% Creation of comparison vector for different noise variances: 
VectorOfNoise = [0:0.1:.5];
name = mfilename;
name = name(1:10);
figure
noise_5_Boxplots(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc,pdfPath)
