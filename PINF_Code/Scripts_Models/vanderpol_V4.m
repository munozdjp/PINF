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

%% Number of variables

%coefficient of the right polinomial (polinomial of order 4 ...
%that perform the desired behavior (increasing))

p1 =         1.1;
p2 =      0.9285;
p3 =    -0.02771;
p4 =   0.0003581;

%%
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
xMax = ones(1,4)
traslation=0; %traslation of the system in y coordinate
%integration of the lorentz system
    mu0=[c(1)]; %initial condition of Beta
     x0=0.51; %Iniitial contions for x and y [x0,y0]        
        y0=0.51;
        [t,x] = ode45(@(t,x) vanderpol_normalized(t,x,c,xMax),tspan,...
            [tspan(1),mu0,x0,y0],options);
                X = [X;x]; 


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

%% plotting

plotState_Beta_time_Pred(tspan,x,mu_observed,yzero,reproduced_data,n,4,'VanDerPol',10)

%% reconstruction with noise

noise_scale = 0.0000001;

[xNoisy, Noise_normalized] = add_Noise_Max(x,noise_scale);

%next step is to use the noise for the figures. 
mu_obsNoisy = mufunc(xNoisy);

[weightdxNoise,resnormNoise,~,exitflagNoise,outputNoise] = lsqcurvefit(F1,weights0, tspan, mu_obsNoisy');

reproduced_dataNoisy= F1(weightdxNoise,tspan);

comparisonVector = [c(1:(length(c)-1));weightdx;weightdxNoise]


plotNoisyStateFig_5(tspan,xNoisy,mu_obsNoisy,yzero,reproduced_dataNoisy,n)

% reconstruction: 
% initial conditions
    C_predicted = c
    maxBetaSadd = 1
    mu0 = yzero(1);
    x0=0.51/xMax(3); 
    y0=0.51/xMax(4);
% insert Model 
[t, xreconstruct]  = ode45(@(t,x) vanderpol_normalized(t,x,C_predicted,xMax),tspan,...
            [tspan(1),mu0,x0,y0],options);


%insert plot: 

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
titleStr = sprintf('$\\mathbf{State \\hspace{0.005\\linewidth} vs \\hspace{0.005\\linewidth} Time, \\hspace{0.005\\linewidth} V =  1e-7}$');
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
save(fullfile(dataDir, 'vanderpol_reconstruction.mat'), 'xNoisy', 'xreconstruct','tspan','noise_scale');

%% Noise analysis 
% Creation of comparison vector for different noise variances: 
VectorOfNoise = [0:0.0001:0.0005];
name = mfilename;
name = 'vanderpol';
%noise_Scale_maximum(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc)
noise_5_Boxplots(F1, weights0, tspan, x,c,yzero, VectorOfNoise, name...
    , mufunc, mfilename, dataDir,plotsDir)
