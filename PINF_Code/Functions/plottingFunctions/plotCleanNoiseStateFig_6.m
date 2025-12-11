function plotCleanNoiseStateFig_6(tspan, Noise_normalized,yzero,n)
    if n == 4 
        figure(6) %figure of the noise for each variable
        subplot(2,2,1)
        hold on
        plot(yzero,Noise_normalized(:,3),'b-')    
        xline(0,'k--');
        yline(0,'k--');
        xlabel('alpha par/ter')
        ylabel('Clean Noise: x3')
        set(gca,'FontSize',16);
        l=legend('Trajectory with noisy');
        l.FontSize = 16;
        l.Location='northeast';
        title('only noise: x3')
        hold off 
        %create subplot to check comparison of the noise
        %then do calculation with several noises
        subplot(2,2,2)
        plot(yzero,Noise_normalized(:,4),'r-')    
        xline(0,'k--');
        yline(0,'k--');
        xlabel('alpha par/ter')
        ylabel('Clean noise x4')
        set(gca,'FontSize',16);
        l=legend(' Noise added');
        l.FontSize = 16;
        l.Location='northeast';
        title('Clean Noise')
        
        subplot(2,2,3)
        hold on
        plot(tspan,Noise_normalized(:,3),'b-')    
        xline(0,'k--');
        yline(0,'k--');
        xlabel('time')
        ylabel('Clean Noise added t x3')
        set(gca,'FontSize',16);
        l=legend('Trajectory with noisy');
        l.FontSize = 16;
        l.Location='northeast';
        title('Clean Noise for x3')
        hold off 
        %create subplot to check comparison of the noise
        %then do calculation with several noises
        subplot(2,2,4)
        plot(tspan,Noise_normalized(:,4),'r-')    
        xline(0,'k--');
        yline(0,'k--');
        xlabel('time')
        ylabel('Clean noise for x4')
        set(gca,'FontSize',16);
        l=legend(' Noise added');
        l.FontSize = 16;
        l.Location='northeast';
        title('Clean Noise')
    elseif n==3
        figure(6) %figure of the noise for each variable
        subplot(2,1,1)
        hold on
        plot(yzero,Noise_normalized(:,3),'b-')    
        xline(0,'k--');
        yline(0,'k--');
        xlabel('alpha par/ter')
        ylabel('Clean Noise: x3')
        set(gca,'FontSize',16);
        l=legend('Trajectory with noisy');
        l.FontSize = 16;
        l.Location='northeast';
        title('only noise: x3')
        hold off 
        %create subplot to check comparison of the noise
        %then do calculation with several noises
        
        subplot(2,1,2)
        hold on
        plot(tspan,Noise_normalized(:,3),'b-')    
        xline(0,'k--');
        yline(0,'k--');
        xlabel('Bifurcation parameter')
        ylabel('Clean Noise added t x3')
        set(gca,'FontSize',16);
        l=legend('Trajectory with noisy');
        l.FontSize = 16;
        l.Location='northeast';
        title('Clean Noise for x3')
        hold off 
        %create subplot to check comparison of the noise
        %then do calculation with several noises        
    
    elseif n==5
        figure(6) %figure of the noise for each variable
        subplot(2,3,1)
        hold on
        plot(yzero,Noise_normalized(:,3),'b-')    
        xline(0,'k--');
        yline(0,'k--');
        xlabel('alpha par/ter')
        ylabel('Clean Noise: x3')
        set(gca,'FontSize',16);
        l=legend('Trajectory with noisy');
        l.FontSize = 16;
        l.Location='northeast';
        title('only noise: x3')
        hold off 
        %create subplot to check comparison of the noise
        %then do calculation with several noises
        subplot(2,3,2)
        plot(yzero,Noise_normalized(:,4),'r-')    
        xline(0,'k--');
        yline(0,'k--');
        xlabel('alpha par/ter')
        ylabel('Clean noise x4')
        set(gca,'FontSize',16);
        l=legend(' Noise added');
        l.FontSize = 16;
        l.Location='northeast';
        title('only noise: x4')
        
        subplot(2,3,3)
        hold on
        plot(yzero,Noise_normalized(:,5),'b-')    
        xline(0,'k--');
        yline(0,'k--');
        xlabel('Bifurcation parameter')
        ylabel('Clean Noise added t x3')
        set(gca,'FontSize',16);
        l=legend('Trajectory with noisy');
        l.FontSize = 16;
        l.Location='northeast';
        title('only noise: x5')
        hold off 
        %create subplot to check comparison of the noise
        %then do calculation with several noises
        subplot(2,3,4)
        plot(tspan,Noise_normalized(:,3),'r-')    
        xline(0,'k--');
        yline(0,'k--');
        xlabel('Biff Par/ter')
        ylabel('Clean noise for x3')
        set(gca,'FontSize',16);
        l=legend(' Noise added');
        l.FontSize = 16;
        l.Location='northeast';
        title('only noise: x3')

        subplot(2,3,5)
        plot(tspan,Noise_normalized(:,4),'r-')    
        xline(0,'k--');
        yline(0,'k--');
        xlabel('Biff Par/ter')
        ylabel('Clean noise for x4')
        set(gca,'FontSize',16);
        l=legend(' Noise added');
        l.FontSize = 16;
        l.Location='northeast';
        title('only noise: x4')

        subplot(2,3,6)
        plot(tspan,Noise_normalized(:,5),'r-')    
        xline(0,'k--');
        yline(0,'k--');
        xlabel('Biff Par/ter')
        ylabel('Clean noise for x5')
        set(gca,'FontSize',16);
        l=legend(' Noise added');
        l.FontSize = 16;
        l.Location='northeast';
        title('only noise: x5')
        

    end


end