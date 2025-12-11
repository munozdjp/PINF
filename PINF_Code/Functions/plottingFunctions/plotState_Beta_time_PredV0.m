function plotState_Beta_time_PredV0(tspan,x,mu_observed,yzero,reproduced_data,...
    n,polyorder,titleX,nfig)
% 9 input arguments. 
    if n == 4
        % figure(nfig)
        figure
        subplot(2,4,1)
        hold on
        plot(tspan,x(:,3),'b-','LineWidth',2)
        xlabel('Time','Interpreter','latex')
        ylabel('X','Interpreter','latex')
        xline(0,'k--');
        set(gca,'FontSize',16);
%         l=legend('x VS T');
        l.FontSize = 14;
        l.Location='NorthWest';
        title(sprintf('Phase for %s', titleX));
        hold off
        
        subplot(2,4,2)
        hold on
        plot(tspan,x(:,4),'b-','LineWidth',2)
        xlabel('Time','Interpreter','latex')
        ylabel('y')
        xline(0,'k--');
        set(gca,'FontSize',16);
%         l=legend('state y');
        l.FontSize = 14;
        l.Location='NorthWest';
        title('State vs time')
        hold off
                    
    
        
        
        subplot(2,4,3)
        plot(tspan,yzero,'r-','linewidth',2)% good True dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_data,'b--','linewidth',2); %\alpha from polinomial fit
        xlabel('Time','Interpreter','latex')
        ylabel('$\alpha$','Interpreter','latex')
        %axis([0 20 -100 5000])
        set(gca,'FontSize',16);
        l=legend('True ','Inferred','Intermediate');
        l.FontSize = 14;
        l.Location='NorthWest';
        title(sprintf('Prediction, PolyOrder = %d', polyorder));
        hold off
        
        
        subplot(2,4,4)
        plot(tspan,yzero,'r-','linewidth',2)% good True dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_data,'b--','linewidth',2); %Alpha from polinomial fit
        plot(tspan,mu_observed,'m-')%,'linewidth',2) %Alpha from steady equations
        xlabel('Time','Interpreter','latex')
        ylabel('$\alpha$','Interpreter','latex')
        %axis([0 20 -100 5000])
        set(gca,'FontSize',16);
        l=legend('True ','Inferred','Intermediate');
        l.FontSize = 14;
        l.Location='NorthWest';
        title(sprintf('Prediction, PolyOrder = %d', polyorder));
        hold off
        
        subplot(2,4,5)
        hold on
        plot(yzero,x(:,3),'b-','LineWidth',2)
        xlabel('$\alpha$','Interpreter','latex')
        ylabel('X','Interpreter','latex')
        xline(0,'k--');
        set(gca,'FontSize',16);
%         l=legend('x VS T');
        l.FontSize = 14;
        l.Location='NorthWest';
        title('x3 vs Bif par/ter')
        hold off
        
        subplot(2,4,6)
        hold on
        plot(yzero,x(:,4),'b-','LineWidth',2)
        xlabel('$\alpha$','Interpreter','latex')
        ylabel('y')
        xline(0,'k--');
        set(gca,'FontSize',16);
%         l=legend('y x VS Time');
        l.FontSize = 14;
        l.Location='NorthWest';
        title('x VS biff/par')
        hold off
    elseif n == 3
        % figure(nfig)
        figure
        subplot(2,3,1)
        hold on
        plot(tspan,x(:,3),'b-','linewidth',2)
        xlabel('Time','Interpreter','latex')
        ylabel('X','Interpreter','latex')
        xline(0,'k--');
        set(gca,'FontSize',16);
        % title('Phase portrait')
        % string1 = 'sigmoid';
        title(sprintf('Phase portrait for %s', titleX));
        hold off
        
        subplot(2,3,2)
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
        title(sprintf('Prediction, PolyOrder = %d', polyorder));

        hold off                
        
        subplot(2,3,3)
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
        title(sprintf('Prediction, PolyOrder = %d', polyorder));
        hold off
        
        subplot(2,3,4)
        hold on
        plot(yzero,x(:,3),'b-','linewidth',2)
        xlabel('$\alpha$','Interpreter','latex')
        ylabel('X','Interpreter','latex')
        xline(0,'k--');
        set(gca,'FontSize',16);
        title('Observation State')
        hold off
    elseif n == 5
        figure(nfig)
        subplot(2,5,1)
        hold on
        plot(tspan,x(:,3),'b-','LineWidth',2)
        xlabel('Time','Interpreter','latex')
        ylabel('X','Interpreter','latex')
        xline(0,'k--');
        set(gca,'FontSize',16);
%         l=legend('X vs $t$','Interpreter','latex');
        l.FontSize = 16;
        l.Location='NorthWest';
        title(sprintf('Phase portrait for %s', titleX));
        hold off
        
        subplot(2,5,2)
        hold on
        plot(tspan,x(:,4),'b-','LineWidth',2)
        xlabel('Time','Interpreter','latex')
        ylabel('y','Interpreter','latex')
        xline(0,'k--');
        set(gca,'FontSize',16);
%         l=legend('y x vs Time');
        l.FontSize = 16;
        l.Location='NorthWest';
        title('y vs time')
        hold off
                    
        subplot(2,5,3)
        hold on
        plot(tspan,x(:,5),'b-')
        xlabel('Time','Interpreter','latex')
        ylabel('x_5','Interpreter','latex')
        xline(0,'k--');
        set(gca,'FontSize',16);
%         l=legend('X_5 vs $t$','interpreter','latex');
        l.FontSize = 16;
        l.Location='NorthWest';
        title('X5 vs $t$','interpreter','latex');
        hold off
        
        
        subplot(2,5,4)
        plot(tspan,yzero,'r-','linewidth',2)% good True dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_data,'b--','linewidth',2); %Alpha from polinomial fit
        xlabel('Time','Interpreter','latex')
        ylabel('$\alpha$','Interpreter','latex')
        %axis([0 20 -100 5000])
        set(gca,'FontSize',16);
        l=legend('True ','Inferred','Intermediate');
        l.FontSize = 14;
        l.Location='NorthWest';
        title('Prediction V = 0')
        hold off
        
        
        subplot(2,5,5)
        plot(tspan,yzero,'r-','linewidth',2)% good True dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_data,'b--','LineWidth',1); %Alpha from polinomial fit
        plot(tspan,mu_observed,'m-')%,'LineWidth',1) %Alpha from steady equations
        xlabel('Time','Interpreter','latex')
        ylabel('$\alpha$','Interpreter','latex')
        %axis([0 20 -100 5000])
        set(gca,'FontSize',16);
        l=legend('True ','Inferred','Intermediate');
        l.FontSize = 14;
        l.Location='NorthWest';
        title(sprintf('Prediction, PolyOrder = %d', polyorder));
        hold off
        
        subplot(2,5,6)
        hold on
        plot(yzero,x(:,3),'b-')
        xlabel('Time','Interpreter','latex')
        ylabel('X','Interpreter','latex')
        xline(0,'k--');
        set(gca,'FontSize',16);
%         l=legend('X3 VS $t$','interpreter','latex');
        l.FontSize = 14;
        l.Location='NorthWest';
        title('$x_3$ vs $\alpha$','Interpreter','latex')
        hold off
        
        subplot(2,5,7)
        hold on
        plot(yzero,x(:,4),'b-')
        xlabel('Time','Interpreter','latex')
        ylabel('y')
        xline(0,'k--');
        set(gca,'FontSize',16);
%         l=legend('y VS Time');
        l.FontSize = 14;
        l.Location='NorthWest';
        title('y vs $\alpha$','Interpreter','latex')
        hold off

        subplot(2,5,8)
        hold on
        plot(yzero,x(:,5),'b-')
        xlabel('Time','Interpreter','latex')
        ylabel('x_5')
        xline(0,'k--');
        set(gca,'FontSize',16);
%         l=legend('x_5  VS Time');
        l.FontSize = 14;
        l.Location='NorthWest';
        title('State variable vs Biffurcation parameter')
        hold off
    end
end
