clear all, close all, clc
figpath = '../figures/';
addpath('../Functions/utils');
addpath('../Functions/violin/');
addpath('../Functions/NoiseFunctions/');
addpath('../Functions/plottingfunctions/');

% Directory for Data
dataDir = fullfile(pwd, 'Data_from_scripts');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

% Plot directory
fullScriptPath = mfilename('fullpath');

% split into folder, name and extension
[scriptFolder, scriptName, scriptExt] = fileparts(fullScriptPath);

% now build your plots_<scriptname> folder
plotsDir = fullfile(scriptFolder, [scriptName '_plots']);
if ~exist(plotsDir,'dir')
    mkdir(plotsDir);
end

omega = 1;
A = 1;
%Number of variables
n = 4;

dt=0.1; %timestep
% Integrate interval
Initinterval = 16;
%neg_inter = 0
neg_inter = -0;
tinterval=[neg_inter Initinterval];
tspan2=[tinterval(2):-dt:tinterval(1)]; %interval for solving the 
tspan=[.001:dt:tinterval(2)-tinterval(1)];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
X=[];
%Coeficients of polyorder 3 c1*t+c2*t^2+c3*t^3 & and final time value 
c=[1,1,1,Initinterval];
%Analitic polinomial
%y=c(1)*(tspan2)+c(2)*(tspan2).^2+c(3)*(tspan2).^3;
%inverted polinomial
yzero=c(1)*(-tspan+tinterval(2) )+c(2)*(-tspan...
    +tinterval(2)).^2+c(3)*(-tspan+tinterval(2)).^3;
maxBeta = max(yzero);
yzero = yzero/maxBeta;

%plot of inverted polinomial
figure(1)
%plot(tspan2,y,'r-',tspan,yzero,'b--')
plot(tspan,yzero,'b--')
xlabel('time')
ylabel('beta parameter')
hold on
xline(0,'k--');
set(gca,'FontSize',16);
l=legend('Explicit polinomial for steady state generation');
l.FontSize = 14;
l.Location='northeast';
title('Inverted polinomial');

% Integration of hopf system
xMax = ones(size(c))
i=0;
xMax = [1,1,1,1];
reescaletime = 1;


% for traslation=[0];
%     i=i+1;
% for mu0=[yzero(1)];
%     %x0=sqrt(yzero(1)/2) %Formula of radius of Initial Condition
%     for x0=sqrt(yzero(1)*maxBeta/2)/xMax(3)+traslation%Formula of radius of Initial Condition
%         [t,x] = ode45(@(t,x) hopfPolyOrder3(t,x,c,omega,A,traslation,xMax,reescaletime,maxBeta),tspan,[tspan(1),mu0,x0,x0],options);
%          X = [X;x];
%          %Variable x3 vs bifurcation parameter
%             figure(2)
%             hold on
%             plot(x(:,2),x(:,3),'b-')
%             xlabel('Beta')
%             ylabel('x_3')
%             xline(0,'k--');
%             set(gca,'FontSize',16);
%             l=legend('Trajectory using numerical integration');
%             l.FontSize = 14;
%             l.Location='northeast';
%             title('State variable vs Biffurcation parameter')
%             hold off    
%     end
% end
% xtrajectories(:,i)=[x(:,3)];
% end
%%
% Assuming hopfPolyOrder3 and other necessary variables are defined
X = []; % Initialize array to store state trajectories
D = []; % Initialize array to store derivatives

traslation = 0;
mu0 = yzero(1);
x0 = sqrt(yzero(1)*maxBeta/2)/xMax(3) + traslation;

[t, x] = ode45(@(t,x) hopfPolyOrder3(t, x, c, omega, A, traslation, xMax, reescaletime, maxBeta), tspan, [tspan(1), mu0, x0, x0], options);
X = [X; x]; % Store state trajectories

% Calculate derivatives for each time step
dxdt = zeros(size(x)); % Initialize matrix to store derivatives
for idx = 1:length(t)
    dxdt(idx, :) = hopfPolyOrder3(t(idx), x(idx, :), c, omega, A, traslation, xMax, reescaletime, maxBeta);
end
D = [D; dxdt]; % Store derivatives
            
            % Plot state variable vs bifurcation parameter
    figure(2)
    hold on
    plot(x(:,2:end), x(:,2:end), 'b-')
    xlabel('Beta')
    ylabel('x_3')
    xline(0, 'k--');
    set(gca, 'FontSize', 16);
    legend('Trajectory using numerical integration', 'FontSize', 14, 'Location', 'northeast');
    title('State variable vs Bifurcation parameter')
    hold off


% Plot derivatives
    figure(12); % New figure for derivatives
    plot(t, D);
    xlabel('Time');
    ylabel('Derivative values');
    title('Derivatives of the system over time');
    legend({'dx_1/dt', 'dx_2/dt', 'dx_3/dt', 'dx_4/dt'}, 'Location', 'best');
    set(gca, 'FontSize', 16);

%%
xMax = max(abs(x))
[x,yzero] = xNormalization(x,yzero);

% Derivative normalized: Take D and apply the same process. 


mufunc = @(x) ((x(:,4)*xMax(4))+...
    (x(:,3)*xMax(3)).*((x(:,3)*xMax(3)).^2+(x(:,4)*xMax(4)).^2))./(x(:,3)*xMax(3)*maxBeta);
mu_observed=mufunc(x);%mu observed from equations. (Beta values!!!!!!)


% polinomial fit
polyOrder = 2
%Define the function to fit: "Poly order 2"
F1=@(weightdx,time) ((-time+tinterval(2))...
    *weightdx(1)+(-time+tinterval(2)).^2*weightdx(2)+(-time+tinterval(2)).^3*weightdx(3))/maxBeta;
%curve fitting
weights0=[0.5 0.5 0.5]; %initial guess for weights
%Prediction of mu(t) depency on t of mu-observed using data integration 
[weightdx, resnorm,~,exitflag,output] = lsqcurvefit(F1,weights0, tspan, mu_observed') 

%data predicted from F1 evaluated with the predicted coefficients
reproduced_data= F1(weightdx,tspan); %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot of State Variable vs Timespan
plotState_Beta_time_Pred(tspan,x,mu_observed,yzero,reproduced_data,n,polyOrder,'Hopf_reconstruct',3)


%% Mu with Noise
noise_scale=0.1;


[xNoisy, Noise_normalized] = add_Noise_Max(x,noise_scale);

figure(6)
plot(tspan,xNoisy(:,3))

%mu observed from equations. (Beta values!!!!!!)
% mu_obsNoisy=(xNoisy(:,4)+xNoisy(:,3).*(xNoisy(:,3).^2+xNoisy(:,4).^2))./(xNoisy(:,3)*maxBeta);
mu_obsNoisy = mufunc(xNoisy);

[weightdxNoise,resnormNoise,~,exitflagNoise,outputNoise] = lsqcurvefit(F1,weights0, tspan, mu_obsNoisy')

reproduced_dataNoisy= F1(weightdxNoise,tspan);

comparisonVector = [c(1:3);weightdx;weightdxNoise]

plotNoisyStateFig_5(tspan,xNoisy,mu_obsNoisy,yzero,reproduced_dataNoisy,n)


X = []
%%

i=0;
C_predicted = weightdxNoise;
C_predicted(4)  = c(4) %we assume we know the initial point of the trajectory, 
%we only perform regression over the coefficients

% Assuming hopfPolyOrder3, C_predicted, and other necessary variables are defined
Xreconstruct = []; % Initialize array to store reconstructed state trajectories
Dreconstruct = []; % Initialize array to store reconstructed derivatives

traslation = 0;
mu0 = yzero(1);
x0 = sqrt(yzero(1)*maxBeta/2)/xMax(3) + traslation;

% [t, x] = ode45(@(t,x) hopfPolyOrder3(t, x, c, omega, A, traslation, xMax, reescaletime, maxBeta), tspan, [tspan(1), mu0, x0, x0], options);

[t, xreconstruct] = ode45(@(t,x) hopfPolyOrder3(t, x, C_predicted, omega, A, traslation, xMax, reescaletime, maxBeta), tspan, [tspan(1), mu0, x0, x0], options);
Xreconstruct = [Xreconstruct; xreconstruct]; % Store reconstructed state trajectories

% Plot reconstructed state variable vs bifurcation parameter
figure(4);
plot(xreconstruct(:,2), xreconstruct(:,3), 'b-');
xlabel('Beta');
ylabel('x_3');
xline(0, 'k--');
set(gca, 'FontSize', 16);
legend('Trajectory using numerical integration', 'FontSize', 14, 'Location', 'northeast');
title('State variable vs Bifurcation parameter');

%trayectory 1 = xnoisy
%model.trayectory  = xreconstruct
R2StateVariables = calculateAvgR2(xNoisy(:,2:end),xreconstruct(:,2:end))
%% Calculate derivative of model and prediction: 
% Trayectory 1 = xnoisy
% Model = same from the xpredicted hence: HopfPolyOrder3 with respective
% model.trayectoryDerivative = evaluate with xreconstruct  
% Example parameters for a specific model
params = {c, omega, A, traslation, xMax, reescaletime, maxBeta}; % Pack parameters into a cell array

% Example call to the modified function
[R2Score,dxdt_noisy,dxdt_model] = function_CalcDerivatives(@hopfPolyOrder3, params, xNoisy, tspan);

% Display the R^2 score
disp(['R^2 Score: ', num2str(R2Score)]);

%% figure(7)
%define figures; 
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
save(fullfile(dataDir, 'Hopf_Reconstruction.mat'), 'xNoisy', 'xreconstruct','tspan','noise_scale');


%% Noise analysis 
% Creation of comparison vector for different noise variances: 
VectorOfNoise = [0:0.1:0.5];
name = mfilename;

underscoreIndex = strfind(name, '_');
if ~isempty(underscoreIndex)
    name(underscoreIndex(1)) = ' ';
end
name = name(1:7);
%noise_Scale_maximum(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc)
noise_5_Boxplots(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,...
    mufunc,mfilename,dataDir,plotsDir)

