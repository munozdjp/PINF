function plotNoisyStateFig_5(tspan,xNoisy,mu_obsNoisy,yzero,reproduced_dataNoisy,n)
    if n == 4
        figure %figure 2 in the preceding example 
        subplot(2,4,1)
        hold on
        plot(tspan,xNoisy(:,3),'b-')
        xlabel('time')
        ylabel('x_3')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('State variable x VS Time');
        l.FontSize = 14;
        l.Location='northeast';
        title('State x3 vs alpha')
        hold off
        
        subplot(2,4,2)
        hold on
        plot(tspan,xNoisy(:,4),'b-')
        xlabel('time')
        ylabel('x_3')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('State x_4 x VS Time');
        l.FontSize = 14;
        l.Location='northeast';
        title('Noisy State x_4 vs alpha')
        hold off
        
        subplot(2,4,3)%graph of the prediction vs real(Add it in Lorentz)
        plot(tspan,yzero,'r-','LineWidth',1)
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_dataNoisy,'b--','LineWidth',1); %Beta from polinomial fit
        xlabel('time')
        ylabel('alpha par/ter')
        set(gca,'FontSize',16);
        l=legend('Real Beta ','Predicted Beta');
        l.FontSize = 16;
        l.Location='northeast';
        title('Prediction ')
        hold off
        
        subplot(2,4,4) %Prediction, intermediate and real 
        plot(tspan,yzero,'r-','LineWidth',1)
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_dataNoisy,'b--','LineWidth',1); %Beta from polinomial fit
        plot(tspan,mu_obsNoisy,'m-')%,'LineWidth',1) %Beta from steady equations
        %axis([13 20 -60 80])
        set(gca,'FontSize',16);
        l=legend('Real Beta ','Predicted Beta','Inter/ate step in alg/thm');
        l.FontSize = 16;
        l.Location='northeast';
        title('Real, Predict, Int/ate')
        hold off
        
        subplot(2,4,5)
        hold on
        plot(yzero,xNoisy(:,3),'b-')    
        xline(0,'k--');
        yline(0,'k--');
        xlabel('Bifurcation parameter')
        ylabel('State variable x')
        set(gca,'FontSize',16);
        l=legend('Trajectory with noisy');
        l.FontSize = 16;
        l.Location='northeast';
        title('Noisy State x3 vs alpha')
        hold off 
        
        %create subplot to check comparison of the noise
        %then do calculation with several noises
        subplot(2,4,6)
        plot(yzero,xNoisy(:,4),'b-')    
        xline(0,'k--');
        yline(0,'k--');
        xlabel('alpha')
        ylabel('state x4 with noise')
        set(gca,'FontSize',16);
        l=legend('Noisy State');
        l.FontSize = 16;
        l.Location='northeast';
        title('Noisy State x4 vs alpha')

     elseif n == 3
        figure
        subplot(2,3,1)
        hold on
        plot(tspan,xNoisy(:,3),'b-')
        xlabel('time')
        ylabel('x_3')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('State x3 VS Time');
        l.FontSize = 14;
        l.Location='northeast';
        title('State x3 vs time')
        hold off
        
        subplot(2,3,2)
        plot(tspan,yzero,'r-','LineWidth',1)% good real dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_dataNoisy,'b--','LineWidth',1); %Beta from polinomial fit
        xlabel('Time')
        ylabel('Beta')
        %axis([0 20 -100 5000])
        set(gca,'FontSize',16);
        l=legend('Real ','Predicted ','Intermediate step ');
        l.FontSize = 14;
        l.Location='northeast';
        title('Real, Predict, Int/ate')
        hold off
        
        
        subplot(2,3,3)
        plot(tspan,yzero,'r-','LineWidth',1)% good real dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_dataNoisy,'b--','LineWidth',1); %Beta from polinomial fit
        plot(tspan,mu_obsNoisy,'m-')%,'LineWidth',1) %Beta from steady equations
        xlabel('Time')
        ylabel('Beta')
        %axis([0 20 -100 5000])
        set(gca,'FontSize',16);
        l=legend('Real ','Predicted ','Intermediate step ');
        l.FontSize = 14;
        l.Location='northeast';
        title('Real, Predict, Int/ate')
        hold off
        
        subplot(2,3,4)
        hold on
        plot(tspan,xNoisy(:,3),'b-')
        xlabel('time')
        ylabel('x_3')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('State variable x VS Time');
        l.FontSize = 14;
        l.Location='northeast';
        title('State x3 vs time')
        hold off
    elseif n == 5
        figure(5)
        subplot(2,5,1)
        hold on
        plot(tspan,xNoisy(:,3),'b-')
        xlabel('time')
        ylabel('x_3')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('State variable x VS Time');
        l.FontSize = 16;
        l.Location='northeast';
        title('State x3 vs time')
        hold off
        
        subplot(2,5,2)
        hold on
        plot(tspan,xNoisy(:,4),'b-')
        xlabel('time')
        ylabel('x_3')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('State x_4 x VS Time');
        l.FontSize = 16;
        l.Location='northeast';
        title('State variable vs time')
        hold off
                    
        subplot(2,5,3)
        hold on
        plot(tspan,xNoisy(:,5),'b-')
        xlabel('time')
        ylabel('x_3')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('State x_5 x VS Time');
        l.FontSize = 16;
        l.Location='northeast';
        title('State variable vs time')
        hold off
        
        
        subplot(2,5,4)
        plot(tspan,yzero,'r-','LineWidth',1)% good real dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_dataNoisy,'b--','LineWidth',1); %Beta from polinomial fit
        xlabel('Time')
        ylabel('Beta')
        %axis([0 20 -100 5000])
        set(gca,'FontSize',16);
        l=legend('Real ','Predicted ','Intermediate step ');
        l.FontSize = 14;
        l.Location='northeast';
        title('Real, Predict, Int/ate')
        hold off
        
        
        subplot(2,5,5)
        plot(tspan,yzero,'r-','LineWidth',1)% good real dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_dataNoisy,'b--','LineWidth',1); %Beta from polinomial fit
        plot(tspan,mu_obsNoisy,'m-')%,'LineWidth',1) %Beta from steady equations
        xlabel('Time')
        ylabel('Beta')
        %axis([0 20 -100 5000])
        set(gca,'FontSize',16);
        l=legend('Real ','Predicted ','Intermediate step ');
        l.FontSize = 14;
        l.Location='northeast';
        title('Real, Predict, Int/ate')
        hold off
        
        subplot(2,5,6)
        hold on
        plot(yzero,xNoisy(:,3),'b-')
        xlabel('alpha')
        ylabel('x_3')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('State variable x3 Vs alpha');
        l.FontSize = 14;
        l.Location='northeast';
        title('State x3 vs alpha')
        hold off
        
        subplot(2,5,7)
        hold on
        plot(yzero,xNoisy(:,4),'b-')
        xlabel('alpha')
        ylabel('x_3')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('State x_4 x VS Time');
        l.FontSize = 14;
        l.Location='northeast';
        title('State x_4 vs alpha')
        hold off

        subplot(2,5,8)
        hold on
        plot(yzero,xNoisy(:,5),'b-')
        xlabel('alpha')
        ylabel('x_3')
        xline(0,'k--');
        set(gca,'FontSize',16);
        l=legend('State x_5 VS Time');
        l.FontSize = 14;
        l.Location='northeast';
        title('State x5 vs alpha')
        hold off
    end
    


end