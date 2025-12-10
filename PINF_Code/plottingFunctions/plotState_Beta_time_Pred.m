function plotState_Beta_time_Pred(tspan,x,mu_observed,yzero,reproduced_data,n,nfig,pdfPath)
    if n == 4
        figure(nfig)
        subplot(2,4,1)
        hold on
        plot(tspan,x(:,3),'b-','LineWidth',2)
        xlabel('Time', 'FontName', 'Times New Roman')  % Capital "Time" and Times New Roman font
        xlim([min(tspan) max(tspan)]);
        ylabel('X', 'FontName', 'Times New Roman')  % Capital "X" and Times New Roman font
        xline(0,'k--');
        
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        title('Observation State', 'FontName', 'Times New Roman')
        hold off
        
        subplot(2,4,2)
        hold on
        plot(tspan,x(:,4),'b-','LineWidth',2)
        xlabel('Time', 'FontName', 'Times New Roman')  % Capital "Time" and Times New Roman font
        ylabel('X', 'FontName', 'Times New Roman')  % Capital "X" and Times New Roman font
        xline(0,'k--');
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        title('Observation State', 'FontName', 'Times New Roman')
        hold off
                    
        subplot(2,4,3)
        plot(tspan,yzero,'r-','linewidth',2)% good True dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_data,'b--','linewidth',2); %Alpha from polinomial fit
        xlabel('Time', 'FontName', 'Times New Roman')  % Capital "Time" and Times New Roman font
        ylabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        if contains(pdfPath, 'Lorentz') 
            ylim([0 6000]);  % Set y-axis limits for 'Lorentz' analysis
        elseif contains(pdfPath, 'FitzHugh')
            ylim([0 1.4])
        end
        
        l1=legend('True ','Inferred','Intermediate');
        l1.FontSize = 14;
        if contains(pdfPath, 'FitzHugh')
            l1.Location='northwest';
        else
            l1.Location='northeast';
        end
        title('Prediction, V = 0', 'FontName', 'Times New Roman')
        hold off
        
        
        subplot(2,4,4)
        plot(tspan,yzero,'r-','linewidth',2)% good True dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_data,'b--','linewidth',2); %Alpha from polinomial fit
        plot(tspan,mu_observed,'m-')%,'linewidth',2) %Alpha from steady equations
        xlabel('Time', 'FontName', 'Times New Roman')  % Capital "Time" and Times New Roman font
        ylabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        if contains(pdfPath, 'Lorentz') || contains(pdfPath, 'Hopf')
            ylim([0 6000]);  % Set y-axis limits for 'Lorentz' or 'hopf' analysis
        elseif contains(pdfPath, 'FitzHugh')
            ylim([0 1.4]);  % Set y-axis limits for 'FitzHugh' analysis
        end
        l2=legend('True ','Inferred','Intermediate');
        l2.FontSize = 14;
        if contains(pdfPath, 'FitzHugh')
            l2.Location='northwest';
        else
            l2.Location='northeast';
        end
        title('Prediction, V = 0', 'FontName', 'Times New Roman')
        hold off
        
        subplot(2,4,5)
        hold on
        plot(x(:,2),x(:,3),'b-','LineWidth',2)
        xlabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font
        ylabel('X', 'FontName', 'Times New Roman')  % Capital "X" and Times New Roman font
        xline(0,'k--');
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        title('Observation State', 'FontName', 'Times New Roman')
        hold off
        
        subplot(2,4,6)
        hold on
        plot(x(:,2),x(:,4),'b-','LineWidth',2)
        xlabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font
        ylabel('X2', 'FontName', 'Times New Roman')  % Capital "X" and Times New Roman font
        xline(0,'k--');
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        title('Observation State', 'FontName', 'Times New Roman')
        hold off
        set(gcf, 'Position', [100, 100, 1000, 600]);
        exportgraphics(gcf, pdfPath, 'ContentType', 'vector', 'Append', true);

    elseif n == 3
        figure(nfig)

        % Subplot 1
        subplot(2,3,1)
        hold on
        plot(tspan,x(:,3),'b-','linewidth',2)
        xlabel('Time', 'FontName', 'Times New Roman')  % Capital "Time" and Times New Roman font
        xlim([min(tspan) max(tspan)]);
        ylabel('X', 'FontName', 'Times New Roman')  % Capital "X" and Times New Roman font
        xline(0,'k--');
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        title('Observation State', 'FontName', 'Times New Roman')
        hold off
        
        % Subplot 2
        subplot(2,3,2)
        plot(tspan,yzero,'r-','linewidth',2)  % True dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_data,'b--','linewidth',2);  % Alpha from polynomial fit
        xlabel('Time', 'FontName', 'Times New Roman')  % Capital "Time" and Times New Roman font
        xlim([min(tspan) max(tspan)]);
        ylabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        l1=legend('True ','Inferred','Intermediate');
        l1.FontSize = 14;
        l1.Location='northwest';
        title('Prediction, V = 0', 'FontName', 'Times New Roman')
        hold off
        
        % Subplot 3
        subplot(2,3,3)
        plot(tspan,yzero,'r-','linewidth',2)  % True dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_data,'b--','linewidth',2);  % Alpha from polynomial fit
        plot(tspan,mu_observed,'m-','linewidth',1)  % Alpha from steady equations
        xlabel('Time', 'FontName', 'Times New Roman')  % Capital "Time" and Times New Roman font
        ylabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font
        xlim([min(tspan) max(tspan)]);
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        l2=legend('True ','Inferred','Intermediate');
        l2.FontSize = 14;
        l2.Location='northwest';
        title('Prediction, V = 0', 'FontName', 'Times New Roman')
        hold off
        
        % Subplot 4
        subplot(2,3,4)
        hold on
        plot(x(:,2),x(:,3),'b-','linewidth',2)
        xlabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font
        xlim([min(x(:,2)) max(x(:,2))]);
        ylabel('X', 'FontName', 'Times New Roman')  % Capital "X" and Times New Roman font
        xline(0,'k--');
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        title('Observation State', 'FontName', 'Times New Roman')
        hold off
        set(gcf, 'Position', [100, 100, 1000, 600]);
        exportgraphics(gcf, pdfPath, 'ContentType', 'vector', 'Append', true);

    elseif n == 5
        figure(nfig)
        subplot(2,5,1)
        hold on
        plot(tspan,x(:,3),'b-','LineWidth',2)
        xlabel('Time', 'FontName', 'Times New Roman')  % Capital "Time" and Times New Roman font
        ylabel('X', 'FontName', 'Times New Roman')  % Capital "X" and Times New Roman font
        xline(0,'k--');
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        title('Observation State', 'FontName', 'Times New Roman')
        hold off
        
        subplot(2,5,2)
        hold on
        plot(tspan,x(:,4),'b-','LineWidth',2)
        xlabel('Time', 'FontName', 'Times New Roman')  % Capital "Time" and Times New Roman font
        ylabel('X_4', 'FontName', 'Times New Roman')  % Capital "X" and Times New Roman font
        xline(0,'k--');
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        title('Observation State', 'FontName', 'Times New Roman')
        hold off
                    
        subplot(2,5,3)
        hold on
        plot(tspan,x(:,5),'b-')
        xlabel('Time', 'FontName', 'Times New Roman')  % Capital "Time" and Times New Roman font
        ylabel('X_5', 'FontName', 'Times New Roman')  % Capital "X" and Times New Roman font
        xline(0,'k--');
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        title('Observation State', 'FontName', 'Times New Roman')
        hold off
        
        
        subplot(2,5,4)
        plot(tspan,yzero,'r-','linewidth',2)% good True dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_data,'b--','linewidth',2); %Alpha from polinomial fit
        xlabel('Time', 'FontName', 'Times New Roman')  % Capital "Time" and Times New Roman font
        ylabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        ylim([0 6000]); % Set y-axis limits
        l1=legend('True ','Inferred','Intermediate');
        l1.FontSize = 14;
        l1.Location='northeast';
        title('Prediction, V = 0', 'FontName', 'Times New Roman')
        hold off
        
        
        subplot(2,5,5)
        plot(tspan,yzero,'r-','linewidth',2)% good True dynamic
        yline(0,'k--','HandleVisibility','off');
        hold on
        plot(tspan,reproduced_data,'b--','LineWidth',1); %Alpha from polinomial fit
        plot(tspan,mu_observed,'m-')%,'LineWidth',1) %Alpha from steady equations
        xlabel('Time', 'FontName', 'Times New Roman')  % Capital "Time" and Times New Roman font
        ylabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        ylim([0 9000]); % Set y-axis limits
        l2=legend('True ','Inferred','Intermediate');
        l2.FontSize = 14;
        l2.Location='northeast';
        title('Prediction, V = 0', 'FontName', 'Times New Roman')
        hold off
        
        subplot(2,5,6)
        hold on
        plot(yzero,x(:,3),'b-')
        xlabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font
        ylabel('X_3', 'FontName', 'Times New Roman')  % Capital "X" and Times New Roman font
        xline(0,'k--');
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        title('Observation State', 'FontName', 'Times New Roman')
        hold off
        
        subplot(2,5,7)
        hold on
        plot(yzero,x(:,4),'b-')
        xlabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font
        ylabel('X_4', 'FontName', 'Times New Roman')  % Capital "X" and Times New Roman font
        xline(0,'k--');
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        title('Observation State', 'FontName', 'Times New Roman')
        hold off

        subplot(2,5,8)
        hold on
        plot(yzero,x(:,5),'b-')
        xlabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol and Times New Roman font
        ylabel('X', 'FontName', 'Times New Roman')  % Capital "X" and Times New Roman font
        xline(0,'k--');
        set(gca,'FontSize',16, 'FontName', 'Times New Roman');
        title('Observation State', 'FontName', 'Times New Roman')
        hold off
        set(gcf, 'Position', [100, 100, 1000, 600]);
        exportgraphics(gcf, pdfPath, 'ContentType', 'vector', 'Append', true);
    end
end
