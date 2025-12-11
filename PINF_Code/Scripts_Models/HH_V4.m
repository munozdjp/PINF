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

n = 5
%===simulation time===
simulationTime = 100; %in milliseconds
deltaT=.01;
t=0:deltaT:simulationTime;
tspan = t

%===specify the external current I===
changeTimes = [0]; %in milliseconds
currentLevels = [50]; %Change this to see effect of different currents on voltage (Suggested values: 3, 20, 50, 1000)

%polinomial structure for the bifurcation parameter
p1 =         1.1;
p2 =      0.9285;
p3 =    -0.02771;
p4 =   0.0003581;

Initinterval = 250;

c=[p1,p2,p3,p4,Initinterval];
polyfunc= @(c) c(1)+c(2)*(t)+c(3)*(t).^2+c(4)*(t).^3;
xMax = [1, 1, 1, 1, 1]; % [I (μA/cm²), n, m, h, V (mV)]



%% Definition of the model
% Input current (assumed provided as a constant or time series)
yzero = polyfunc(c)
I = polyfunc(c); % Input current in μA/cm²

% === Constant parameters ===
% All of these can be found in Table 3
gbar_K = 36;   % mS/cm²
gbar_Na = 120; % mS/cm²
g_L = 0.3;     % mS/cm²
E_K = -12;     % mV
E_Na = 115;    % mV
E_L = 10.6;    % mV
C = 1;         % μF/cm²

% === Normalization parameters ===
% Define maximum values for normalization

% xMax(1): Max input current (adjust based on your context)
% xMax(2:4): Max for n, m, h (gating variables, max = 1)
% xMax(5): Max for V (set to ~E_Na for max depolarization)
reescaletime = 1; % Time scaling factor (set to 1 if no rescaling needed)

% === Set the initial states ===
% Initialize normalized states
V_norm = 0 / xMax(5); % Normalize baseline voltage (V = 0 mV)
alpha_n = 0.01 * ((10 - V_norm * xMax(5)) / (exp((10 - V_norm * xMax(5)) / 10) - 1)); % Equation 12
beta_n = 0.125 * exp(-V_norm * xMax(5) / 80); % Equation 13
alpha_m = 0.1 * ((25 - V_norm * xMax(5)) / (exp((25 - V_norm * xMax(5)) / 10) - 1)); % Equation 20
beta_m = 4 * exp(-V_norm * xMax(5) / 18); % Equation 21
alpha_h = 0.07 * exp(-V_norm * xMax(5) / 20); % Equation 23
beta_h = 1 / (exp((30 - V_norm * xMax(5)) / 10) + 1); % Equation 24

n_norm(1) = alpha_n / (alpha_n + beta_n); % Equation 9
m_norm(1) = alpha_m / (alpha_m + beta_m); % Equation 18
h_norm(1) = alpha_h / (alpha_h + beta_h); % Equation 18

% Initialize arrays for normalized states
V_norm = zeros(1, numel(t));
n_norm = zeros(1, numel(t));
m_norm = zeros(1, numel(t));
h_norm = zeros(1, numel(t));
V_norm(1) = 0 / xMax(5); % Set initial normalized voltage

% If I is time-varying, normalize it
I_norm = I / xMax(1); % Normalize input current

for i = 1:numel(t)-1 % Compute coefficients, currents, and derivatives at each time step
    % --- Denormalize states for calculations ---
    V = V_norm(i) * xMax(5); % Denormalize voltage
    n = n_norm(i) * xMax(2); % Denormalize n
    m = m_norm(i) * xMax(3); % Denormalize m
    h = h_norm(i) * xMax(4); % Denormalize h

    % --- Calculate the coefficients ---
    alpha_n(i) = 0.01 * ((10 - V) / (exp((10 - V) / 10) - 1));
    beta_n(i) = 0.125 * exp(-V / 80);
    alpha_m(i) = 0.1 * ((25 - V) / (exp((25 - V) / 10) - 1));
    beta_m(i) = 4 * exp(-V / 18);
    alpha_h(i) = 0.07 * exp(-V / 20);
    beta_h(i) = 1 / (exp((30 - V) / 10) + 1);

    % --- Calculate the currents ---
    I_Na = (m^3) * gbar_Na * h * (V - E_Na); % Equations 3 and 14
    I_K = (n^4) * gbar_K * (V - E_K); % Equations 4 and 6
    I_L = g_L * (V - E_L); % Equation 5
    I_ion = I_norm(i) * xMax(1) - I_K - I_Na - I_L; % Use normalized I

    % --- Calculate the derivatives using Euler first order approximation ---
    % Normalize the derivatives by dividing by xMax and apply time scaling
    V_norm(i+1) = V_norm(i) + deltaT * (I_ion / C) * reescaletime / xMax(5);
    n_norm(i+1) = n_norm(i) + deltaT * (alpha_n(i) * (1 - n) - beta_n(i) * n) * reescaletime / xMax(2);
    m_norm(i+1) = m_norm(i) + deltaT * (alpha_m(i) * (1 - m) - beta_m(i) * m) * reescaletime / xMax(3);
    h_norm(i+1) = h_norm(i) + deltaT * (alpha_h(i) * (1 - h) - beta_h(i) * h) * reescaletime / xMax(4);
end

% --- Output normalized states ---
% Match the output format x = [I', n', m', h', V']
x = [I_norm', n_norm', m_norm', h_norm', V_norm'];
% define the max variables and normalize variables 
xMax = max(abs(x));
[x,yzero] = xNormalization(x,yzero);

mufunc = @(x) ...
    ((x(:,3)*xMax(3)).^3)   .* gbar_Na .* (x(:,4)*xMax(4))   .* ((x(:,5)*xMax(5)) - E_Na) ...
  + ((x(:,2)*xMax(2)).^4)   .* gbar_K  .* ((x(:,5)*xMax(5)) - E_K)  ...
  + g_L .* ((x(:,5)*xMax(5)) - E_L);


% construction of the hidden variable
mu_observed = mufunc(x)


F1=@(weightdx,time) weightdx(1)+weightdx(2)*time+time.^2*weightdx(3)+...
    time.^3*weightdx(4);

% curve fitting
weights0=[0.5 0.5 0.5 0.5];
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, t, mu_observed')

% data predicted from F1 evaluated in the polinomail fit
reproduced_data= F1(weightdx,t);


% x = [t',y2',V'];
plotState_Beta_time_Pred(t,x,mu_observed,yzero,reproduced_data,5,3,'Hodgkin-Huxley',10)

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

plotNoisyStateFig_5(tspan,xNoisy,mu_obsNoisy,yzero,reproduced_dataNoisy,5)

%% Prediction on the model itself. 

I = polyfunc(weightdxNoise); % Input current in μA/cm²

% === Constant parameters ===
% All of these can be found in Table 3
gbar_K = 36;   % mS/cm²
gbar_Na = 120; % mS/cm²
g_L = 0.3;     % mS/cm²
E_K = -12;     % mV
E_Na = 115;    % mV
E_L = 10.6;    % mV
C = 1;         % μF/cm²

% === Normalization parameters ===
% Define maximum values for normalization

% xMax(1): Max input current (adjust based on your context)
% xMax(2:4): Max for n, m, h (gating variables, max = 1)
% xMax(5): Max for V (set to ~E_Na for max depolarization)
reescaletime = 1; % Time scaling factor (set to 1 if no rescaling needed)

% === Set the initial states ===
% Initialize normalized states
V_norm = 0 / xMax(5); % Normalize baseline voltage (V = 0 mV)
alpha_n = 0.01 * ((10 - V_norm * xMax(5)) / (exp((10 - V_norm * xMax(5)) / 10) - 1)); % Equation 12
beta_n = 0.125 * exp(-V_norm * xMax(5) / 80); % Equation 13
alpha_m = 0.1 * ((25 - V_norm * xMax(5)) / (exp((25 - V_norm * xMax(5)) / 10) - 1)); % Equation 20
beta_m = 4 * exp(-V_norm * xMax(5) / 18); % Equation 21
alpha_h = 0.07 * exp(-V_norm * xMax(5) / 20); % Equation 23
beta_h = 1 / (exp((30 - V_norm * xMax(5)) / 10) + 1); % Equation 24

n_norm(1) = alpha_n / (alpha_n + beta_n); % Equation 9
m_norm(1) = alpha_m / (alpha_m + beta_m); % Equation 18
h_norm(1) = alpha_h / (alpha_h + beta_h); % Equation 18

% Initialize arrays for normalized states
V_norm = zeros(1, numel(t));
n_norm = zeros(1, numel(t));
m_norm = zeros(1, numel(t));
h_norm = zeros(1, numel(t));
V_norm(1) = 0 / xMax(5); % Set initial normalized voltage

% If I is time-varying, normalize it
I_norm = I / xMax(1); % Normalize input current

for i = 1:numel(t)-1 % Compute coefficients, currents, and derivatives at each time step
    % --- Denormalize states for calculations ---
    V = V_norm(i) * xMax(5); % Denormalize voltage
    n = n_norm(i) * xMax(2); % Denormalize n
    m = m_norm(i) * xMax(3); % Denormalize m
    h = h_norm(i) * xMax(4); % Denormalize h

    % --- Calculate the coefficients ---
    alpha_n(i) = 0.01 * ((10 - V) / (exp((10 - V) / 10) - 1));
    beta_n(i) = 0.125 * exp(-V / 80);
    alpha_m(i) = 0.1 * ((25 - V) / (exp((25 - V) / 10) - 1));
    beta_m(i) = 4 * exp(-V / 18);
    alpha_h(i) = 0.07 * exp(-V / 20);
    beta_h(i) = 1 / (exp((30 - V) / 10) + 1);

    % --- Calculate the currents ---
    I_Na = (m^3) * gbar_Na * h * (V - E_Na); % Equations 3 and 14
    I_K = (n^4) * gbar_K * (V - E_K); % Equations 4 and 6
    I_L = g_L * (V - E_L); % Equation 5
    I_ion = I_norm(i) * xMax(1) - I_K - I_Na - I_L; % Use normalized I

    % --- Calculate the derivatives using Euler first order approximation ---
    % Normalize the derivatives by dividing by xMax and apply time scaling
    V_norm(i+1) = V_norm(i) + deltaT * (I_ion / C) * reescaletime / xMax(5);
    n_norm(i+1) = n_norm(i) + deltaT * (alpha_n(i) * (1 - n) - beta_n(i) * n) * reescaletime / xMax(2);
    m_norm(i+1) = m_norm(i) + deltaT * (alpha_m(i) * (1 - m) - beta_m(i) * m) * reescaletime / xMax(3);
    h_norm(i+1) = h_norm(i) + deltaT * (alpha_h(i) * (1 - h) - beta_h(i) * h) * reescaletime / xMax(4);
end

% --- Output normalized states ---
% Match the output format x = [I', n', m', h', V']
xreconstruct = [I_norm', n_norm', m_norm', h_norm', V_norm'];
% define the max variables and normalize variables 

%%
h = figure('Name','Reconstruction');

hold on

plot(tspan,xNoisy(:,5),'b-')
plot(tspan,xreconstruct(:,5),'r-')
xlabel('$Time$','interpreter','latex')
ylabel('State','interpreter','latex')
xline(0,'k--');
legend('Original Noisy Data', 'Reconstructed Trajectory', 'Location', 'best') % Add a legend
set(gca,'FontSize',16);
% title(['$\mathbf{State \hspace{0.005\linewidth} vs \hspace{0.005\linewidth} Time, \hspace{0.005\linewidth} V = noise_scale}$'],'Interpreter','latex')
titleStr = sprintf('$\\mathbf{State \\hspace{0.005\\linewidth} vs \\hspace{0.005\\linewidth} Time, \\hspace{0.005\\linewidth} V = %.1f}$', noise_scale);
title(titleStr,'Interpreter','latex')
hold off    

%% 
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


%%

% save('LorentzOri_Reco.mat', 'xNoisy', 'xreconstruct','tspan','noise_scale');
save(fullfile(dataDir, 'HH_reconstruction.mat'), 'xNoisy', 'xreconstruct','tspan','noise_scale');
%% Noise analysis 
% Creation of comparison vector for different noise variances: 
VectorOfNoise = [0:0.1:0.5];
name = mfilename;
%noise_Scale_maximum(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc)
noise_5_Boxplots(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,...
    mufunc,mfilename,dataDir,plotsDir)
