function noise_5_Boxplots(F1,weights0, tspan, x,c,yzero,VectorOfNoise,name,mufunc,pdfPath)
    violinmatrix = [];
    iteration = 20;
    for i = [1:iteration]
        %Noise scale;
        comparisonVectorCum = [];
        distanceVector = [];
        forindex = 0;
        indexMiddlePlot = floor(length(VectorOfNoise)/2);
        %indexToPlot = [1, indexMiddlePlot,length(VectorOfNoise)];
        indexToPlot = [1:1:length(VectorOfNoise)];
        columnnsFig = length(indexToPlot);
        rowsFig = 4;
        counterif = 0;
        % figure
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
                %for pitchfork is the plot (5)
                %for saddle is the plot (4)
                
                figure(gcf)

                % Subplot 1
                subplot(rowsFig, columnnsFig, counterif)
                hold on
                plot(x(:,2),xNoisy(:,3),'b-')    
                xline(0,'k--');
                xlim([min(x(:,2)) max(x(:,2))]);
                yline(0,'k--');
                xlabel('\alpha', 'FontName', 'Times New Roman')  % X label with Times New Roman font
                ylabel('X', 'FontName', 'Times New Roman')  % Alpha symbol for the ylabel with Times New Roman font
                set(gca,'FontSize',16, 'FontName', 'Times New Roman');
                % legend('Trajectory with noisy', 'FontName', 'Times New Roman', 'FontSize', 16, 'Location', 'northeast');
                title(['State X, V = ', num2str(noise_scale)], 'FontName', 'Times New Roman')
                hold off 
                
                % Subplot 2
                subplot(rowsFig, columnnsFig, counterif + columnnsFig)
                plot(x(:,2),Noise_normalized(:,3),'r-')    
                xline(0,'k--');
                xlim([min(x(:,2)) max(x(:,2))]);
                yline(0,'k--');
                xlabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol for the xlabel with Times New Roman font
                ylabel('Amplitude', 'FontName', 'Times New Roman')  % Y label with Times New Roman font
                set(gca,'FontSize',16, 'FontName', 'Times New Roman');
                % legend('Noise', 'FontName', 'Times New Roman', 'FontSize', 16, 'Location', 'northeast');
                title(['N(0, ', num2str(noise_scale), ')'], 'FontName', 'Times New Roman', 'Interpreter', 'latex')  % Use LaTeX interpreter for math format
                
                % Subplot 3
                subplot(rowsFig, columnnsFig, counterif + columnnsFig * 2)
                plot(tspan,yzero,'r-', 'LineWidth', 1)
                yline(0,'k--','HandleVisibility','off');
                xlabel('Time', 'FontName', 'Times New Roman')  % Time label with Times New Roman font
                ylabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol for the ylabel with Times New Roman font
                hold on
                plot(tspan,reproduced_dataNoisy,'b--','LineWidth',1); % Beta from polynomial fit
                set(gca,'FontSize',16, 'FontName', 'Times New Roman');
                xlim([min(tspan) max(tspan)]);
                l1 = legend('Real Beta', 'Prediction', 'FontName', 'Times New Roman', 'FontSize', 12, 'Location', 'northeast');
                % legend('Real Beta', 'Prediction', 'FontName', 'Times New Roman', 'FontSize', 16, 'Location', 'northeast');
                l1.Position = [0.92, 0.4, 0.06, 0.052]; %(Near the top left)    
                title(['Prediction = ', num2str(noise_scale)], 'FontName', 'Times New Roman')
                hold off
                
                % Subplot 4
                subplot(rowsFig, columnnsFig, counterif + columnnsFig * 3)
                plot(tspan,yzero,'r-', 'LineWidth', 1)
                yline(0,'k--','HandleVisibility','off');
                xlabel('Time', 'FontName', 'Times New Roman')  % Time label with Times New Roman font
                ylabel('\alpha', 'FontName', 'Times New Roman')  % Alpha symbol for the ylabel with Times New Roman font
                hold on
                plot(tspan,reproduced_dataNoisy,'b--','LineWidth', 1); % Beta from polynomial fit
                plot(tspan,mu_obsNoisy,'m-')%,'LineWidth',1) % Beta from steady equations
                xlim([min(tspan) max(tspan)]);
                if contains(pdfPath, 'Hopf')
                    ylim([0 9000]);  % Set y-axis limits for 'Lorentz' or 'hopf' analysis
                % elseif contains(pdfPath, 'FitzHugh')
                %     ylim([0 1.4]);  % Set y-axis limits for 'FitzHugh' analysis
                end
                set(gca,'FontSize',16, 'FontName', 'Times New Roman');
                l2 = legend('Real Beta', 'Prediction', 'Intermediate step', 'FontName', 'Times New Roman', 'FontSize', 12, 'Location', 'northeast'); % Reduced font size
                % legend('Real Beta', 'Prediction', 'Intermediate step', 'FontName', 'Times New Roman', 'FontSize', 12, 'Location', 'northeast'); % Reduced font size
                % Manually set the legend position
                % l.Position = [x y width height];
                % Example: l.Position = [0.7, 0.7, 0.2, 0.1]; Adjust the values accordingly
                l2.Position = [0.92, 0.19, 0.06, 0.052]; %(Near the top left)
                title(['Noise = ', num2str(noise_scale)], 'FontName', 'Times New Roman')
                hold off
                % Global title
                sgtitle(name, 'FontName', 'Times New Roman')

            end
        end
        
        violinmatrix = [violinmatrix, distanceVector];
        if i == iteration
            set(gcf, 'Position', [100, 100, 1600, 1200]);
            exportgraphics(gcf, pdfPath, 'ContentType', 'vector', 'Append', true);
            figure;
            % If you need to use subplot, uncomment and adjust the next line accordingly
            % subplot(rowsFig, columnnsFig, [rowsFig * columnnsFig - length(rowsFig):1:rowsFig * columnnsFig])
            plot(VectorOfNoise, distanceVector);
            xlabel('Noise Variance', 'FontName', 'Times New Roman', 'FontSize', 16);  % Set X label with Times New Roman font
            ylabel('L2 Error', 'FontName', 'Times New Roman', 'FontSize', 16);  % Corrected Y label text and set font
            title('Iteration = 20', 'FontName', 'Times New Roman', 'FontSize', 16);  % Set a fixed title with Times New Roman font
            
            set(gca, 'FontSize', 16);  % Set the font size of the axes
            
            % l = legend('Error');  % Define the legend to label the plot data
            % l.FontSize = 14;      % Set font size of the legend
            % l.Location = 'northeast';  % Position the legend in the top right corner
            exportgraphics(gcf, pdfPath, 'ContentType', 'vector', 'Append', true);
        end
    end
    % print('Predicted weights for each variance')
    VectorAnd_distances = [[c(1:(length(c)-1)); comparisonVectorCum],[0; distanceVector]]

    %% [h,L,MX,MED]=violin(violinmatrix')%,'xlabel',[0:0.2:1.8])
    figure;
    % If subplot or other configuration is needed, include it here
    arraystring = arrayfun(@num2str, VectorOfNoise, 'UniformOutput', false);
    h = boxplot(violinmatrix', 'Labels', arraystring);
    % h = boxplot(violinmatrix', 'Labels', arraystring, 'Whisker', Inf, 'Symbol', ''); % Inf whisker length and no outlier plot
    % h = boxplot(violinmatrix', 'Labels', arraystring, 'Symbol', ''); % Inf whisker length and no outlier plot
    
    % Change line width and color
    set(h, {'LineWidth', 'Color'}, {2, [0 0 0]}); % Set all lines to width 2 and color black
    
    % Set properties for axes
    set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
    
    % Label and title configurations
    ylabel('Error IHCV', 'FontSize', 16, 'FontName', 'Times New Roman');
    xlabel('Noise variance ', 'FontSize', 16, 'FontName', 'Times New Roman');
    title([name, ' Bifurcation-Normalized Error '], 'FontSize', 16, 'FontName', 'Times New Roman');
    
    % To ensure all parts of the box are black, additional handle adjustments are needed:
    boxColors = findobj(gca, 'Tag', 'Box'); % Find the boxes
    set(boxColors, 'Color', [0 0 0]); % Set color of boxes to black
    
    medLines = findobj(gca, 'Tag', 'Median'); % Find the median lines
    set(medLines, 'Color', [0 0 0]); % Set color of median lines to black
    
    whiskers = findobj(gca, 'Tag', 'Whisker'); % Find the whiskers
    set(whiskers, 'Color', [0 0 0]); % Set color of whiskers to black
    
    outliers = findobj(gca, 'Tag', 'Outliers'); % Find the outliers
    set(outliers, 'MarkerEdgeColor', [0 0 0]); % Set color of outliers to black

    exportgraphics(gcf, pdfPath, 'ContentType', 'vector', 'Append', true);
end