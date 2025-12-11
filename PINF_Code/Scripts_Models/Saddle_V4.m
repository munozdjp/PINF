%The sadddle note bifurcation has problem when converging to the origin
%i.e. moving from left to right(it goes to infinity after reachin the biff
%point)

clear all, close all, clc
figpath = '../figures/';
addpath('../Functions/utils');
addpath('../Functions/violin/');
addpath('../Functions/NoiseFunctions/');
addpath('../Functions/plottingfunctions/');

%data directory
dataDir = fullfile(pwd, 'Data_from_scripts');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

%plot directory
fullScriptPath = mfilename('fullpath');

% split into folder, name and extension
[scriptFolder, scriptName, scriptExt] = fileparts(fullScriptPath);

% now build your plots_<scriptname> folder
plotsDir = fullfile(scriptFolder, [scriptName '_plots']);
if ~exist(plotsDir,'dir')
    mkdir(plotsDir);
end

% Number of variables
n = 3;

dt=0.1; %timestep
% Integrate interval
% Initinterval = 100;
%  Initinterval = 2;
% neg_inter = 0;
Initinterval = 16;
neg_inter = 0;
tinterval=[neg_inter Initinterval];
tspan2=[tinterval(2):-dt:tinterval(1)]; %interval for solving the 
tspan=[.01:dt:tinterval(2)-tinterval(1)];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
X=[];

%Coeficients of polyorder 3 c1*t+c2*t^2+c3*t^3 & and final time value 
c=[1,0,0,Initinterval];
%Analitic polynomial
%y=c(1)*(tspan2)+c(2)*(tspan2).^2+c(3)*(tspan2).^3;
y=c(1)*(tspan)+c(2)*(tspan).^2+c(3)*(tspan).^3;
%derivative of analitical polynomial is : c1+2*c(2)*(tspan2)+3*c(3)*(tspan2).^2;

%inverted polynomial
y2=c(1)*(-tspan+tinterval(2) )+c(2)*(-tspan...
    +tinterval(2)).^2+c(3)*(-tspan+tinterval(2)).^3;

yzero = c(1)*(tspan+tinterval(1))+c(2)*(tspan+tinterval(1)).^2 ...
        +c(3).*(tspan+tinterval(1)).^3;
maxBeta = max(yzero)
yzero = yzero/maxBeta
%The input of the system is the derivative of y2
% = -c(1)- 2*c(2)* (-y(1)+c(4)) - 3*c(3)*(-y(1)+c(4))^2;
%plot of inverted polynomial
figure(1)
%plot(tspan2,y,'r-',tspan,y2,'b--')
plot(tspan,y,'r-',tspan,y2,'b--')
xlabel('Time', 'Interpreter', 'latex')
ylabel('$\alpha$', 'Interpreter', 'latex')
% ylabel('beta parameter')
hold on
xline(0,'k--');
set(gca,'FontSize',16);
l=legend('Explicit','Implicit');
l.FontSize = 14;
l.Location='North';
title('Inverted polynomial');
hold off
%Integration of hopf system

%reescaling is zero for my origina variables. 
xMax = [1,1,1];

for mu0=[yzero(1)]    
    %x0=sqrt(y2(1)/2) %Formula of radius of Initial Condition
    %for x0=sqrt(y2(1))%Formula of radius of Initial Condition
    for x0=sqrt(yzero(1))%Formula of radius of Initial Condition
        [t,x] = ode45(@(t,x) saddlenodeNewPoli(t,x,c,xMax,1,maxBeta),tspan,[tspan(1),mu0,x0],options);
         X = [X;x];
         % Variable x3 vs bifurcation parameter
            figure(2)  
            subplot(1,3,1)
            hold on
            plot(x(:,2),x(:,3),'b-')
            xlabel('Beta')
            ylabel('x_3')
            xline(0,'k--');
            set(gca,'FontSize',16);
            l=legend('Trajectory');
            l.FontSize = 14;
            l.Location='northeast';
            title('State vs Beta for LV system')
            hold off               

            subplot(1,3,2)
            hold on
%             plot(tspanSaddle,x(:,2),'b-')
            plot(tspan,x(:,2),'b-')
            xlabel('time')
            ylabel('beta')
            xline(0,'k--');
            set(gca,'FontSize',16);
            l=legend('Trajectory');
            l.FontSize = 14;
            l.Location='northeast';
            title('Beta vs time')
            hold off

            subplot(1,3,3)
            hold on
            plot(tspan,x(:,3),'b-')
            xlabel('time')
            ylabel('x_3')
            xline(0,'k--');
            set(gca,'FontSize',16);
            l=legend('Trajectory');
            l.FontSize = 14;
            l.Location='northeast';
            title('Beta vs time')
            hold off
    end
end
% Section 2 
% the interval of mu goes to maximum 1. 


% Normalization variables: 

xMax = max(abs(x));
[x,yzero] = xNormalization(x,yzero);

mufunc = @(x) ((x(:,3).^2)*xMax(3)^2)/maxBeta; 

mu_observed= mufunc(x);

% Section 3 
% Polynomial fit
polyOrder = 2
%Define the function to fit: "Poly order 2"
F1=@(weightdx,xdata)  (weightdx(1)*(xdata+tinterval(1))+weightdx(2)*(xdata+tinterval(1)).^2 ...
        +weightdx(3).*(xdata+tinterval(1)).^3)/maxBeta;

%curve fitting
weights0=[0.5 0.5 0.5]; %initial guess for weights
opts = optimoptions('lsqcurvefit','Display','off');
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, tspan, mu_observed',[],[],opts);

%data predicted from F1 evaluated with the predicted coefficients
reproduced_data= F1(weightdx,tspan);

%plot of beta: real VS polynomial fit VS steady equations
plotState_Beta_time_Pred(tspan,x,mu_observed,yzero,reproduced_data,n,polyOrder,'SaddleReconstruct',3)

%% Mu with ratio State-Noise

%noise normal multiplied by 0.1
%noise_scale = 0.1
noise_scale = 0


%noise has the same variaance over all the trayectory. 
[xNoisy, Noise_normalized] = add_Noise_Max(x,noise_scale);

% next step is to use the noise for the figures. 
% (prediction of mu with noisy data)
mu_obsNoisy = mufunc(xNoisy);

[weightdxNoise,resnormNoise,~,exitflagNoise,outputNoise] = lsqcurvefit(F1,weights0, tspan, mu_obsNoisy');

reproduced_dataNoisy= F1(weightdxNoise,tspan);

comparisonVector = [c(1:(length(c)-1));weightdx;weightdxNoise]

% plotNoisyStateFig_4(tspan,xNoisy,mu_obsNoisy,yzero,reproduced_dataNoisy,n)

plotNoisyStateFig_5(tspan,xNoisy,mu_obsNoisy,yzero,reproduced_dataNoisy,n)

%% plotCleanNoiseStateFig_5(tspan,Noise_normalized,yzero,n)
% Here we change the c by the original mu_noisy
X = []

C_predicted = weightdxNoise
% Assuming yzero, C_predicted, xMax, tspan, and options are already defined

% Directly use the singular values for mu0 and x0
mu0 = yzero(1);
% x0 = sqrt(yzero(1)); % Formula of radius of Initial Condition
x0 = 0.025070610528195
reescaletimeSadd  = 1
% maxBetaSadd = 1
maxBetaSadd = maxBeta
%Here saddle is not normalizez in the second variable ReescaledTimeSadd.
% This code is to compare with the original one.
% Integrate the system[t,x] = ode45(@(t,x) saddlenodeNewPoli(t,x,c,xMax,1,maxBeta),tspan,[tspan(1),mu0,x0],options);
[t, xreconstruct] = ode45(@(t, x) saddlenodeNewPoli(t, x, C_predicted, xMax, reescaletimeSadd , maxBetaSadd), tspan, [tspan(1), mu0, x0], options);

% Plot Variable x3 vs bifurcation parameter
figure(6)
subplot(1,3,1)
hold on
plot(xreconstruct(:,2), xreconstruct(:,3), 'b-')
xlabel('$\alpha$', 'Interpreter', 'latex')
ylabel('State', 'Interpreter', 'latex')
xline(0, 'k--');
set(gca, 'FontSize', 16);
title(['$\mathbf{State \hspace{0.005\linewidth} vs \hspace{0.005\linewidth}\alpha(t)}$'], 'Interpreter', 'latex')
hold off

subplot(1,3,2)
            hold on
            plot(tspan,xreconstruct(:,2),'b-')
            xlabel('time')
            ylabel('beta')
            xline(0,'k--');
            set(gca,'FontSize',16);
            l=legend('Trajectory');
            l.FontSize = 14;
            l.Location='northeast';
            title('Beta vs time')
            hold off

            subplot(1,3,3)
            hold on
            plot(tspan,xreconstruct(:,3),'b-')
            xlabel('time')
            ylabel('x_3')
            xline(0,'k--');
            set(gca,'FontSize',16);
            l=legend('Trajectory');
            l.FontSize = 14;
            l.Location='northeast';
            title('Beta vs time')
            hold off
%% calculation of R^2 score sfor state variables
R2StateVariables = calculateAvgR2(xNoisy(:,2:end),xreconstruct(:,2:end))

% Calculation of r^2score of derivatives: 
params = {C_predicted, xMax, reescaletimeSadd , maxBetaSadd}; % Pack parameters into a cell array

%% Calculation of R^2 score for derivatives
% Example call to the modified function
% [R2Score,dxdt_noisy,dxdt_model] = function_CalcDerivatives(@saddlenodeNewPoli, params, xNoisy, tspan);

%% ============================================================
% Derivative comparison (inline version of function_CalcDerivatives)
% ============================================================

% Example inputs (make sure these exist in your workspace before running)
% model = @saddlenodeNewPoli;
% params = {...};         % your parameter list
% xNoisy = ...;           % trajectory data, e.g. 160x4
% tspan = ...;            % matching time vector

dxdt_model = zeros(size(xNoisy));

for idx = 1:length(tspan)
    dxdt_model(idx, :) = saddlenodeNewPoli(tspan(idx), xNoisy(idx, :), params{:});
end

% Finite difference derivatives
dxdt_fd = diff(xNoisy) ./ diff(tspan');
dxdt_fd = [dxdt_fd; dxdt_fd(end, :)]; % Extend last row to match dimensions

% --- Compute R² metrics ---
addpath('/Users/munozdjp/Library/CloudStorage/GoogleDrive-juan.munozdiaz@kaust.edu.sa/My Drive/PHD_Kaust_DriveFolder/Academic_KAUST/ACDC_project/projects1/HINNDy project 2022/NoiseFunctions');
[R2Score, r2perVariable] = calculateAvgR2(dxdt_fd(:,end:end), dxdt_model(:,end:end));

% --- Display results ---
fprintf('Global R² score: %.6f\n', R2Score);
for i = 1:length(r2perVariable)
    fprintf('R² for variable %d: %.6f\n', i, r2perVariable(i));
end

% --- Plot comparison ---
nVars = size(dxdt_model, 2);
colors  = ['b', 'r', 'g', 'm'];
markers = ['o', '*', '+', 'x'];

figure;
hold on;
for i = 1:nVars
    plot(tspan, dxdt_model(:,i), ['-', colors(i)], 'LineWidth', 1.5, ...
        'DisplayName', ['Model dx', num2str(i), '/dt']);
    plot(tspan, dxdt_fd(:,i), ['--', colors(i), markers(i)], 'LineWidth', 1, ...
        'MarkerSize', 4, 'DisplayName', ['FD dx', num2str(i), '/dt']);
end
title('Comparison of Model vs Finite Difference Derivatives');
xlabel('Time');
ylabel('Derivative Value');
legend('Location','best');
grid on;
hold off;

% Display the R^2 score
disp(['R^2 Score: ', num2str(R2Score)]);

%%
h = figure('Name','Reconstruction');

hold on
% plot(x(:,2),x(:,3),'b-')
plot(tspan,xNoisy(:,3),'b-')
plot(tspan,xreconstruct(:,3),'r-')
xlabel('$Time$','interpreter','latex')
ylabel('State','interpreter','latex')
xline(0,'k--');
legend('Original Noisy Data', 'Reconstructed Trajectory', 'Location', 'best') % Add a legend
set(gca,'FontSize',16);
% title(['$\mathbf{State \hspace{0.005\linewidth} vs \hspace{0.005\linewidth} Time, \hspace{0.005\linewidth} V = noise_scale}$'],'Interpreter','latex')
titleStr = sprintf('$\\mathbf{State \\hspace{0.005\\linewidth} vs \\hspace{0.005\\linewidth} Time, \\hspace{0.005\\linewidth} V = %.1f}$', noise_scale);
title(titleStr,'Interpreter','latex')
hold off    

        figName = get(h,'Name');

        if isempty(figName)
            % no Name property set → use figure number
            safeName = sprintf('Figure_%d', h.Number);
        else
            % sanitize the given Name
            safeName = regexprep(figName,'[^\w]','_');
        end

    set(h, 'Renderer', 'painters');                        % force vector renderer
    filename = fullfile(plotsDir, [safeName '.pdf']);     % build path
    exportgraphics(h, filename, ...
        'ContentType', 'vector', ...                       % vector PDF
        'BackgroundColor', 'none', ...                     % transparent background
        'Append', false);                                  % overwrite if exists

save(fullfile(dataDir, 'Saddle_Reconstruction.mat'), 'xNoisy', 'xreconstruct','tspan','noise_scale');

%% Creation of plot that add noise to the reconstructed signal. 
% Noise analysis 
% Creation of comparison vector for different noise variances: 

VectorOfNoise = [0:0.1:0.5];
name = mfilename;

%noise_Scale_maximum(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc)

% dataDir: file to save the scripts
noise_5_Boxplots(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc,mfilename,dataDir,plotsDir)