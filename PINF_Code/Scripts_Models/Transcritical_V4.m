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

%Note: Every-Modification in 4 needs to be reflected in 2
%% 1) Generate data for the Pitchfork
clear all, close all, clc
figpath = '../figures/';
addpath('../Functions/utils');
addpath('../Functions/violin/');
addpath('../Functions/NoiseFunctions/');
addpath('../Functions/plottingfunctions/');

tic
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

%% Number of variables
%generate Data

n = 3; %number of variables
dt=0.01; %timestep

%Integrate interval

Endinterval = 5; %initial value of the pitchfor bifurcation
%Endinterval = 1.5; %initial value of the pitchfor bifurcation
%intervalzero = -.001;
intervalzero = 0;

%real hidden variable interval
%tinterval=[0 Initinterval];
tinterval=[intervalzero Endinterval];

% Integrate
tspan=[.001:dt:tinterval(2)-tinterval(1)];
N = length(tspan);

options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
Xpitch=[];
%Coeficients of polyorder 3 c1*t+c2*t^2+c3*t^3 & and final time value 
%c=[3,2,1,Endinterval];
polyorder= 1
cPitch=[1,1,0,Endinterval];
c = cPitch
%Analitic polinomial

%Integration for Pitchforksystem
reescale = 1
%xMax = [1,1,1];

yzNonreescaled=c(1)*(tspan+tinterval(1))+c(2)*(tspan+tinterval(1)).^2 ...
        +c(3).*(tspan+tinterval(1)).^3;
yzero = reescale*(c(1)*(tspan+tinterval(1))+c(2)*(tspan+tinterval(1)).^2 ...
        +c(3).*(tspan+tinterval(1)).^3);

% maxBetapitch = max(yzero);
% yzero = yzero/max(yzero);
maxBetapitch = 1;
%yzero = yzero/max(yzero);

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
        [t,x] = ode45(@(t,x) transcriticalReescaled(t,x,c,xMax,reescale,maxBetapitch ),tspan,[tspan(1),mu0,x0]);%,options);
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
%% Section 2 

% Normalization variables: 

xMax = max(abs(x));
[x,yzero] = xNormalization(x,yzero);
       
mufunc = @(x) ((x(:,3)))*(xMax(3))/maxBetapitch; %here max is one

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
        
titleX = 'Transcritical'
plotState_Beta_time_Pred(tspan,x,mu_observed,yzero,reproduced_data,n,polyorder,'Transcritical',3)



        h  = figure('Name','phaseportrait')
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

       h =  figure('Name','prediction')
        plot(tspan,yzero,'r-','linewidth',2)% good True dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_data,'b--','linewidth',2); %\alpha from polinomial fit
        xlabel('Time','Interpreter','latex')
        ylabel('$\alpha$','Interpreter','latex')
        Ymin = min([yzero, reproduced_data]) - 0.005;
        Ymax = max([yzero, reproduced_data]) + 0.005;
        ylim([Ymin Ymax])
        xlim([tspan(1) tspan(end)])
        %axis([0 20 -100 5000])
        set(gca,'FontSize',16);
        l=legend('True ','Inferred','Intermediate');
        l.FontSize = 14;
        l.Location='NorthWest';
        title(sprintf('Prediction, PolyOrder = %d', polyorder));

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

%% 

% Mu with ratio State-Noise

%noise normal multiplied by 0.1
noise_scale = 0.1;


%noise has the same variaance over all the trayectory. 
[xNoisy, Noise_normalized] = add_Noise_Max(x,noise_scale);

%next step is to use the noise for the figures. 
 %(prediction of mu with noisy data)
mu_obsNoisy = mufunc(xNoisy);

[weightdxNoise,resnormNoise,~,exitflagNoise,outputNoise] = lsqcurvefit(F1,weights0, tspan, mu_obsNoisy');

reproduced_dataNoisy= F1(weightdxNoise,tspan);

comparisonVector = [c(1:(length(c)-1));weightdx;weightdxNoise]
% Evaluate the 

% plotNoisyStateFig_4(tspan,xNoisy,mu_obsNoisy,yzero,reproduced_dataNoisy,n)

plotNoisyStateFig_5(tspan,xNoisy,mu_obsNoisy,yzero,reproduced_dataNoisy,n)
 % plotCleanNoiseStateFig_5(tspan,Noise_normalized,yzero,n)
%% Here we change the c by the original mu_noisy
X = []

C_predicted = weightdxNoise
% Assuming yzero, C_predicted, xMax, tspan, and options are already defined

% Directly use the singular values for mu0 and x0
mu0 = yzero(1);
x0 = sqrt(yzero(1)); % Formula of radius of Initial Condition
reescaletimeSadd  = 1
maxBetaSadd = 1
% Integrate the system
        % [t,x] = ode45(@(t,x) transcriticalReescaled(t,x,c,xMax,reescale,maxBetapitch ),tspan,[tspan(1),mu0,x0]);%,options);

[t, xreconstruct] = ode45(@(t, x) transcriticalReescaled(t, x, C_predicted, xMax, reescaletimeSadd , maxBetaSadd), tspan, [tspan(1), mu0, x0], options);

% Plot Variable x3 vs bifurcation parameter
h = figure('Name','Plain_Reconstruction2')
hold on
plot(xreconstruct(:,2), xreconstruct(:,3), 'b-')
xlabel('$\alpha$', 'Interpreter', 'latex')
ylabel('State', 'Interpreter', 'latex')
xline(0, 'k--');
set(gca, 'FontSize', 16);
title(['$\mathbf{State \hspace{0.005\linewidth} vs \hspace{0.005\linewidth}\alpha(t)}$'], 'Interpreter', 'latex')
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

R2StateVariables = calculateAvgR2(xNoisy(:,2:end),xreconstruct(:,2:end))

% Calculation of r^2score of derivatives: 
params = {C_predicted, xMax, reescaletimeSadd , maxBetaSadd}; % Pack parameters into a cell array

% Example call to the modified function
[R2Score,dxdt_noisy,dxdt_model] = function_CalcDerivatives(@saddlenodeNewPoli, params, xNoisy, tspan);

% Display the R^2 score
disp(['R^2 Score: ', num2str(R2Score)]);

h = figure('Name','Reconstruction')
hold on
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


% % save('SaddleOri_Reco.mat', 'xNoisy', 'xreconstruct','tspan','noise_scale');
save(fullfile(dataDir, 'Saddle_Reconstruction.mat'), 'xNoisy', 'xreconstruct','tspan','noise_scale');

%% Creation of plot that add noise to the reconstructed signal. 
% Noise analysis 
% Creation of comparison vector for different noise variances: 

VectorOfNoise = [0:0.1:0.5];
name = mfilename;

%noise_Scale_maximum(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc)

noise_5_Boxplots(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,...
    mufunc,mfilename,dataDir,plotsDir)
