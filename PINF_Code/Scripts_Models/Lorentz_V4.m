clear all, close all, clc
figpath = '../figures/';
addpath('../Functions/utils');
addpath('../Functions/violin/');
addpath('../Functions/NoiseFunctions/');
addpath('../Functions/plottingfunctions/');
dataDir = fullfile(pwd, 'Data_from_scripts');

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

% p1 =         1.1;
% p2 =      1.9285;
% p3 =    2.02771;
% p4 =   0.43581;

omega= 10;
beta=2.66;
rho=1.1;
n = 5; % Number of variables
% dt=0.01; % Timestep
dt=0.1; % Timestep
%dt=1; % Timestep

% Integrate Interval
Initinterval = 250;
%Initinterval = 160;
neg_inter = 0;
tinterval=[neg_inter Initinterval];
tspan=[.1:dt:tinterval(2)-tinterval(1)];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
X=[];

%Coeficients of polyorder 3; p1+p2*t+p3*t^2+p4*t^3 & and initial time value 
c=[p1,p2,p3,p4,Initinterval];
yzero=c(1)+c(2)*(tspan)+c(3)*(tspan).^2+c(4)*(tspan).^3;
maxBeta = max(yzero);
yzero = yzero/maxBeta;
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

xMax = ones(size(c))
xMax = [1,1,1,1,1];
reescaletime = 1;

traslation=0; %traslation of the system in y coordinate
% integration of the lorentz system
% i=i+1
mu0=[yzero(1)]; %initial condition of Beta
x0=0.51; %Iniitial contions for x and y [x0,x0,z0]        
y0=0.51;
z0=0.1;  %initial condition for z0
[t,x] = ode45(@(t,x) lorentz2(t,x,c,omega,beta,traslation,xMax,reescaletime,maxBeta),tspan,[tspan(1),mu0,x0,y0,z0],options);
 X = [X;x];
%             Plot of state variable vs Beta parameter
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


%Normalization
xMax = max(abs(x));
[x,yzero] = xNormalization(x,yzero);

mufunc = @(x) ((x(:,3)*xMax(3)).*(x(:,5)*xMax(5))+(x(:,4)*xMax(4)))./(x(:,3)*xMax(3)*maxBeta);
mu_observed = mufunc(x);

%Define the function to fi: "Poly order 2"
F1=@(weightdx,time) (weightdx(1)+weightdx(2)*time+time.^2*weightdx(3)+...
    time.^3*weightdx(4))/maxBeta;

%Curve fitting
weights0=[0.5 0.5 0.5 0.5];
opts = optimoptions('lsqcurvefit','Display','off');
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, tspan, mu_observed',[],[],opts);

%Data predicted from F1 evaluated in the polinomail fit
reproduced_data= F1(weightdx,tspan);

plotState_Beta_time_Pred(tspan,x,mu_observed,yzero,reproduced_data,n,3,'Lorentz',15)
%% Mu with ratio State-Noise

noise_scale = 0.5;

[xNoisy, Noise_normalized] = add_Noise_Max(x,noise_scale);

%next step is to use the noise for the figures. 
mu_obsNoisy = mufunc(xNoisy);

[weightdxNoise,resnormNoise,~,exitflagNoise,outputNoise] = lsqcurvefit(F1,weights0, tspan, mu_obsNoisy');

reproduced_dataNoisy= F1(weightdxNoise,tspan);

comparisonVector = [c(1:(length(c)-1));weightdx;weightdxNoise]

plotNoisyStateFig_5(tspan,xNoisy,mu_obsNoisy,yzero,reproduced_dataNoisy,n)

X = []
traslation=0; %traslation of the system in y coordinate
C_predicted = weightdxNoise
% integration of the lorentz system
%     i=i+1
mu0=[yzero(1)]; %initial condition of Beta
    %Iniitial contions for x and y [x0,x0,z0]
x0=0.51/xMax(3);         
y0=0.51/xMax(4);
z0=0.1/xMax(5);  %initial condition for z0
[t,xreconstruct] = ode45(@(t,x) lorentz2(t,x,C_predicted,omega,beta,traslation,xMax,reescaletime,maxBeta),tspan,[tspan(1),mu0,x0,y0,z0],options);
X = [X;xreconstruct];
% Plot of state variable vs Beta parameter
figure(2)
hold on
plot(xreconstruct(:,2),xreconstruct(:,3),'b-')
xlabel('Beta')
ylabel('x_3')
xline(0,'k--');
set(gca,'FontSize',16);
l=legend('Trajectory using numerical integration');
l.FontSize = 14;
l.Location='northeast';
title('State variable vs Biffurcation parameter')
hold off    

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
save(fullfile(dataDir, 'Lorentz_Reconstruction.mat'), 'xNoisy', 'xreconstruct','tspan','noise_scale');
%% Noise analysis 
% Creation of comparison vector for different noise variances: 
VectorOfNoise = [0:0.1:0.5];
name = mfilename;
%noise_Scale_maximum(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc)
resolution = 'normal';
noise_5_Boxplots(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc,mfilename,dataDir,plotsDir,resolution)

