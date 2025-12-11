function plotState_Beta_time_Pred_without_hidden(tspan,x,mu_observed,yzero,reproduced_data,n,nfig)
figure(nfig)
subplot(2,4,1)
hold on
plot(tspan,x(:,2),'b-','linewidth',2)
xlabel('Time','Interpreter','latex')
ylabel('X','Interpreter','latex')
xline(0,'k--');
set(gca,'FontSize',16);
title('Phase portrait')
hold off
%%%

subplot(2,4,2)
hold on
plot(tspan,x(:,3),'b-','linewidth',2)
xlabel('Time','Interpreter','latex')
ylabel('X','Interpreter','latex')
xline(0,'k--');
set(gca,'FontSize',16);
title('Phase portrait')
hold off
%%%%%

subplot(2,4,3)
plot(tspan,yzero,'r-','linewidth',2)% good True dynamic
yline(0,'k--','HandleVisibility','off');
hold on
plot(tspan,reproduced_data,'b--','linewidth',2); %Alpha from polinomial fit
xlabel('Time','Interpreter','latex')
ylabel('$\alpha$','Interpreter','latex')
%axis([0 20 -100 5000])
set(gca,'FontSize',16);
l=legend('True ','Inferred');
l.FontSize = 14;
l.Location='NorthWest';
title('Prediction, V = 0')

hold off                

subplot(2,4,4)
plot(tspan,yzero,'r-','linewidth',2)% good True dynamic
yline(0,'k--','HandleVisibility','off');
hold on
plot(tspan,reproduced_data,'b--','linewidth',2); %Alpha from polinomial fit
plot(tspan,mu_observed,'m-','linewidth',1)%,'linewidth',2) %Alpha from steady equations
xlabel('Time','Interpreter','latex')
ylabel('$\alpha$','Interpreter','latex')
%axis([0 20 -100 5000])
set(gca,'FontSize',16);
l=legend('True ','Inferred','Intermediate');
l.FontSize = 14;
l.Location='NorthWest';
title('Prediction, V = 0')
hold off

subplot(2,4,5)
hold on
plot(yzero,x(:,3),'b-','linewidth',2)
xlabel('$\alpha$','Interpreter','latex')
ylabel('X','Interpreter','latex')
xline(0,'k--');
set(gca,'FontSize',16);
title('Observation State')
hold off