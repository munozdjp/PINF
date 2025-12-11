clear all, clc, clear

clear all, close all, clc
figpath = '../figures/';
addpath('../Functions/utils');
addpath('../Functions/violin/');
addpath('../Functions/NoiseFunctions/');
addpath('../Functions/plottingfunctions/');

%% data directory
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

%% Updated second code to match good code behavior
%% Corrected MATLAB Code for Lotka-Volterra System

% Define parameters
npar = 0.25;                % Nonlinear parameter
n = 4;                      % Number of variables
dt = 0.01;                  % Time step
Initinterval = 56;          % End of time interval
neg_inter = 0;              % Start of time interval
tinterval = [neg_inter Initinterval];
tspan = 0:dt:(tinterval(2) - tinterval(1)); % Time span for simulation
N = length(tspan);          % Number of time points
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12*ones(1,n)); % ODE solver options
XLV = [];                   % Array to store results

% Polynomial coefficients for beta parameter
polyorder = 1;
cLV = [0.3, -0.001, 0, 0, Initinterval];
c = cLV;

% Time scaling and beta polynomial
scaletimeLV = 1;
yzero = scaletimeLV * (c(1) + c(2)*tspan + c(3)*tspan.^2 + c(4)*tspan.^3);
maxBetaLV = 1

% Plot the beta parameter polynomial
figure(1)
plot(tspan, yzero, 'r-')
xlabel('time')
ylabel('beta parameter')
hold on
xline(0, 'k--');
set(gca, 'FontSize', 16);
l = legend('Explicit Poly');
l.FontSize = 14;
l.Location = 'northeast';
title('Beta Parameter Polynomial');
mu0 = yzero(1);             % Initial beta value


% Initial conditions (normalized)maxBetaLV = 1;              % Maximum beta scaling factor
% xMax = [1, 1, 50/16, 50/16];      % Normalization factors for variables
% x0 = 15.6/16 ;        % Normalized initial condition for x3
% y0 = 19.65/16;       % Normalized initial condition for x4

% xMax = [1, 1, 50/50, 50/50];      % Normalization factors for variables
% x0 = 15.6/50 ;        % Normalized initial condition for x3
% y0 = 19.65/50;       % Normalized initial condition for x4
xMax = [1, 1, 1, 1];      % Normalization factors for variables
x0 = (15.6/50) * xMax(3);        % Normalized initial condition for x3
y0 = (19.65/50)* xMax(4);       % Normalized initial condition for x4

        % Solve the ODE using the normalized Lotka-Volterra function
        [t, x] = ode45(@(t,x) lotkavolterrav1(t, x, c, npar, xMax, ...
            scaletimeLV, maxBetaLV), tspan, [tspan(1), mu0, x0, y0], options);
        XLV = [XLV; x];
        
        % Plotting the results
        h = figure('Name', 'Lotka-Volterra Phase Plots');
        subplot(2,3,1)
        hold on
        plot(tspan, x(:,3), 'b-')
        xlabel('time')
        ylabel('x_3 (normalized)')
        xline(0, 'k--');
        set(gca, 'FontSize', 16);
        l = legend('Trajectory x3');
        l.FontSize = 14;
        l.Location = 'northeast';
        title('x3 vs Time')
        
        subplot(2,3,2)
        plot(tspan, x(:,4), 'b-')
        xlabel('time')
        ylabel('x_4 (normalized)')
        xline(0, 'k--');
        set(gca, 'FontSize', 16);
        l = legend('Trajectory x4');
        l.FontSize = 14;
        l.Location = 'northeast';
        title('x4 vs Time')
        
        subplot(2,3,3)
        plot(tspan, x(:,2), 'b-')
        xlabel('time')
        ylabel('beta')
        xline(0, 'k--');
        set(gca, 'FontSize', 16);
        l = legend('Trajectory beta');
        l.FontSize = 14;
        l.Location = 'northeast';
        title('Beta vs Time')
        
        subplot(2,3,4)
        plot(x(:,2), x(:,3), 'b-')
        xlabel('beta')
        ylabel('x_3 (normalized)')
        xline(0, 'k--');
        set(gca, 'FontSize', 16);
        l = legend('Trajectory x3');
        l.FontSize = 14;
        l.Location = 'northeast';
        title('x3 vs Beta')
        
        subplot(2,3,5)
        plot(x(:,2), x(:,4), 'b-')
        xlabel('beta')
        ylabel('x_4 (normalized)')
        xline(0, 'k--');
        set(gca, 'FontSize', 16);
        l = legend('Trajectory x4');
        l.FontSize = 14;
        l.Location = 'northeast';
        title('x4 vs Beta')
        hold off

%normalization
xMax = max(abs(x));
[x,yzero] = xNormalization(x,yzero);
%%
%mufunc = @(x) ((x(:,4)+xtras)/(xscale) - betatras)/maxBetaLV;
mufunc = @(x) ((x(:,3))*xMax(3))/maxBetaLV;

mu_observed = mufunc(x); 

F1=@(weightdx,time) weightdx(1)+(time)*weightdx(2)+...
    (time).^2*weightdx(3)+(time).^3*weightdx(4);

% Curve fitting
weights0=[0.5 0.5 0.5 0.5]; %initial guess for weights
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, tspan, mu_observed');

%data predicted from F1 evaluated with the predicted coefficients
reproduced_data= F1(weightdx,tspan);

titleX = 'Lotka Volterra'

plotState_Beta_time_Pred(tspan,x,mu_observed,yzero,reproduced_data,n,polyorder,titleX,10)
% find all the figures that are currently open

        figure('Name','phaseportrait')
        plot(tspan,x(:,3),'b-','LineWidth',2)
        xlabel('Time','Interpreter','latex')
        ylabel('X','Interpreter','latex')
        xline(0,'k--');
        xlim([tspan(1) tspan(end)])
        set(gca,'FontSize',16);
%         l=legend('x VS T');
        l.FontSize = 14;
        l.Location='NorthWest';
        title(sprintf('Phase for %s', titleX));
        

        figure('Name','prediction')
        plot(tspan,yzero,'r-','linewidth',2)% good True dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_data,'b--','linewidth',2); %\alpha from polinomial fit
        % plot(tspan,x(:,2),'k--','linewidth',3); %\alpha from polinomial fit
        xlabel('Time','Interpreter','latex')
        ylabel('$\alpha$','Interpreter','latex')
        Ymin = min([yzero, reproduced_data]) - 0.005;
        Ymax = max([yzero, reproduced_data]) + 0.005;
        ylim([Ymin Ymax])
        xlim([tspan(1) tspan(end)])
        %axis([0 20 -100 5000])
        set(gca,'FontSize',16);
        l=legend('True ','Inferred','from system');
        l.FontSize = 14;
        l.Location='NorthWest';
        title(sprintf('Prediction, PolyOrder = %d', polyorder));
        



%% Mu with ratio State-Noise

%noise normal multiplied by 0.1
noise_scale = 0.1


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
%%

X = []
% xMax = ones(1,4)
C_predicted = weightdxNoise
x0 = (15.6/50) * 1/xMax(3);        % Normalized with the inverse because 
% the model is normalized the inversed way
y0 = (19.65/50)* 1/xMax(4);       % Normalized initial condition for x4
mu0 = yzero(1);             % Initial beta value

[t, xreconstruct] = ode45(@(t,x) lotkavolterrav1(t, x, c, npar, 1./xMax, ... %model normalized with the inverse since we are using inverse normalization. 
            scaletimeLV, maxBetaLV), tspan, [tspan(1), mu0, x0, y0], options);

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

save(fullfile(dataDir, 'Volterra_reconstruction.mat'), 'xNoisy', 'xreconstruct','tspan','noise_scale');

%% Noise analysis 
% Creation of comparison vector for different noise variances: 
VectorOfNoise = [0:0.1:0.5];
name = mfilename;

% noise_Scale_maximum(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc)
noise_5_Boxplots(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,...
    mufunc,mfilename,dataDir,plotsDir)

