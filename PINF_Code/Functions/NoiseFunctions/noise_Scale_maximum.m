function noise_Scale_maximum(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc)
    violinmatrix = [];
    iteration = 20;
    for i = [1:iteration]
        %Noise scale;
        comparisonVectorCum = [];
        distanceVector = [];
        forindex = 0;
        indexMiddlePlot = floor(length(VectorOfNoise)/2);
        indexToPlot = [1, indexMiddlePlot,length(VectorOfNoise)];
        columnnsFig = length(indexToPlot);
        rowsFig = 4;
        counterif = 0;
        for noise_scale = VectorOfNoise  
            forindex = forindex+1;
            %Create a for loop for different noise structures: 
            [xNoisy, Noise_normalized] = add_Noise_Max(x,noise_scale);
        %End of noise scaling        
            mu_obsNoisy = mufunc(xNoisy);
            opts = optimset('Display','off');
            %opts= optimset('Display','on');
            [weightdxNoise,resnormNoise,~,exitflagNoise,outputNoise] = lsqcurvefit(F1,weights0, tspan, mu_obsNoisy',[],[],opts);
            
            reproduced_dataNoisy= F1(weightdxNoise,tspan);
            
            comparisonVectorCum = [comparisonVectorCum;weightdxNoise];

            %This is the cost function I need to change. 
            %the real data is yzero
            % The approximation is reproducedatanoisy
            
            % Metric with trajectories:
            %check first why I am gettting normalized the real variable in 
            %fig 4 and fig 5
            distance = norm(reproduced_dataNoisy-yzero);
            distanceVector = [distanceVector;distance];            

            %Metric with coefficient
%             distance = norm(weightdxNoise-c(1:length(c)-1));
%             distanceVector = [distanceVector;distance];
            % Here I need to create the figures for each iteration
            if sum(forindex == indexToPlot)==1 && i==iteration
                counterif = counterif+1;
                figure(8)%create 
                subplot(rowsFig,columnnsFig,counterif)
                hold on
                plot(x(:,2),xNoisy(:,3),'b-')    
                xline(0,'k--');
                yline(0,'k--');
                xlabel('Bifurcation parameter')
                ylabel('State variable x')
                set(gca,'FontSize',16);
                l=legend('Trajectory with noisy');
                l.FontSize = 16;
                l.Location='northeast';
                title(['Noisy State vs Biff/n par/ter with Variance = ',num2str(noise_scale)])
                hold off 
                %create subplot to check comparison of the noise
                %then do calculation with several noises
                subplot(rowsFig, columnnsFig, counterif+columnnsFig)
                plot(x(:,2),Noise_normalized(:,3),'r-')    
                xline(0,'k--');
                yline(0,'k--');
                xlabel('Bifurcation Parameter')
                ylabel('Only Noise added over  x')
                set(gca,'FontSize',16);
                l=legend('Normal Noise added');
                l.FontSize = 16;
                l.Location='northeast';
                title('Normal Noise')
                
                subplot(rowsFig, columnnsFig ,counterif+columnnsFig*2)
                plot(tspan,yzero,'r-','LineWidth',1)
                yline(0,'k--','HandleVisibility','off');
                xlabel('time')
                ylabel('Biff/tion Par/ter')
                hold on
                plot(tspan,reproduced_dataNoisy,'b--','LineWidth',1); %Beta from polinomial fit
                set(gca,'FontSize',16);
                l=legend('Real Beta ','Predicted Beta');
                l.FontSize = 16;
                l.Location='northeast';
                title(['Prediction with noisy data ',num2str(noise_scale)])
                hold off
                
                subplot(rowsFig, columnnsFig ,counterif+columnnsFig*3)
                plot(tspan,yzero,'r-','LineWidth',1)
                yline(0,'k--','HandleVisibility','off');
                xlabel('time')
                ylabel('Biff/tion Par/ter')
                hold on
                plot(tspan,reproduced_dataNoisy,'b--','LineWidth',1); %Beta from polinomial fit
                plot(tspan,mu_obsNoisy,'m-')%,'LineWidth',1) %Beta from steady equations
                %axis([13 20 -60 80])
                set(gca,'FontSize',16);
                l=legend('Real Beta ','Predicted Beta','Intermediate step in alg/thm');
                l.FontSize = 16;
                l.Location='northeast';
                title(['Prediction with noisy data',num2str(noise_scale)])
                hold off
                sgtitle(name) 
            end
        end
    
        violinmatrix = [violinmatrix,distanceVector];
        if i == iteration
            figure(9)
            %subplot(rowsFig,columnnsFig,[rowsFig*columnnsFig-lenght(rowsFig):1:rowsFig*columnnsFig])
            plot(VectorOfNoise,distanceVector)
            xlabel('Noise Variance')
            ylabel('L2 Error of polynomial coefficiente')
            title(['error vs Noise , iteration = ',num2str(iteration)])
            set(gca,'FontSize',16);
            l=legend('error');
            l.FontSize = 14;
            l.Location='northeast';
        end
    end
    print('Predicted weights for each variance')
    VectorAnd_distances = [[c(1:(length(c)-1));comparisonVectorCum],[0;distanceVector]]

    %% [h,L,MX,MED]=violin(violinmatrix')%,'xlabel',[0:0.2:1.8])
    figure(10)
    %subplot(rowsFig,columnnsFig,[rowsFig*columnnsFig-length(rowsFig):1:rowsFig*columnnsFig])
    arraystring = arrayfun(@num2str,VectorOfNoise,'Uni',0);
    h = boxplot(violinmatrix','Labels',arraystring); 
    set(h,{'linew'},{2})
    ylabel('Error HINNDy','FontSize',16)
    xlabel('Noise variance ','FontSize',16)
    title([name,' Bifurcation-Normalized Error '],'FontSize',16)
    set(gca,'FontSize',16);

end