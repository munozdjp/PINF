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
maxBeta = max(yzero);
yzero = yzero/maxBeta;

figure(1)
plot(tspan2,y,'r-',tspan,y2,'b--',tspan,yzero,'k--')

xline(0,'k--');
set(gca,'FontSize',16);
l=legend('Explicit polinomial for steady state generation','implicit polinmial for ODE data generation', 'polinomial moved in the positive axis');
l.FontSize = 14;
l.Location='northeast';
title('Inverted polinomial');
%Integration for Pitchforksystem
xMax = [1,1,1];
reescaletime = 1

for mu0=[yzero(1)];  %initial condition of Beta  
    %sqrt(muo)%formula of coordinate x0 in function of beta(mu0)
    for x0=real(sqrt(mu0))+0.01 %Initial condition of state variables
        [t,x] = ode45(@(t,x) pitchforkPolyorder3Left2Right(t,x,c,xMax,reescaletime,maxBeta),tspan,[tspan(1),mu0,x0],options);
         X = [X;x];
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
    end
end
%%
xMax = max(abs(x))

%Mu observed from equations of steady state
[x,yzero] = xNormalization(x,yzero);

mufunc = @(x) ((x(:,3).^2)*xMax(3)^2)/maxBeta;
mu_observed=mufunc(x); 
%%Polinomial fit

%In the fututure we can let unkonw the coefficient (tinterval(1)or tinterval(2))
F1=@(weightdx,xdata)  (weightdx(1)*(xdata+tinterval(1))+weightdx(2)*(xdata+tinterval(1)).^2 ...
        +weightdx(3).*(xdata+tinterval(1)).^3)/maxBeta;


%curve fitting
weights0=[0.5 0.5 0.5];%initial guess for weights
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, tspan, mu_observed')

%data predicted from F1 evaluated with the predicted coefficients
reproduced_data= F1(weightdx,tspan);


plotState_Beta_time_Pred(tspan,x,mu_observed,yzero,reproduced_data,n,4,'Pitchfork',10)

%% Data with noise: 

noise_scale = 0.1;
% Here we add noise that has the same variance over all the trayectory
[xNoisy, Noise_normalized] = add_Noise_Max(x,noise_scale);

%next step is to use the noise for the figures. 
mu_obsNoisy = mufunc(xNoisy);

[weightdxNoise,resnormNoise,~,exitflagNoise,outputNoise] = lsqcurvefit(F1,weights0, tspan, mu_obsNoisy');

reproduced_dataNoisy= F1(weightdxNoise,tspan);

comparisonVector = [c(1:(length(c)-1));weightdx;weightdxNoise]

X = []

C_predicted = weightdxNoise
for mu0=[yzero(1)];  %initial condition of Beta  
    %sqrt(muo)%formula of coordinate x0 in function of beta(mu0)
    for x0=real(sqrt(mu0))+0.01 %Initial condition of state variables
        [t,xreconstruct] = ode45(@(t,x) pitchforkPolyorder3Left2Right(t,x,C_predicted,xMax,reescaletime,maxBeta),tspan,[tspan(1),mu0,x0],options);
         X = [X;xreconstruct];
         %plot variable x vs bifurcation parameter
            figure(3)
            hold on
            plot(xreconstruct(:,2),xreconstruct(:,3),'b-')  
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
    end
end

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

% save('PitchforkOri_Reco.mat', 'xNoisy', 'xreconstruct','tspan','noise_scale');
save(fullfile(dataDir, 'Pitchfork_Reconstruction.mat'), 'xNoisy', 'xreconstruct','tspan','noise_scale');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plotCleanNoiseStateFig_6(tspan,Noise_normalized,yzero,n)


%% Noise analysis 
% Creation of comparison vector for different noise variances: 
VectorOfNoise = [0:0.1:0.5];
name = mfilename;

%noise_Scale_maximum(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc)
noise_5_Boxplots(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc,mfilename,dataDir,plotsDir)

