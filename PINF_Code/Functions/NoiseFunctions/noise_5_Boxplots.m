% The lines: 
% xline(0,'k--'); yline(0,'k--');
% has been commented since 

function noise_5_Boxplots(...
    F1, weights0, tspan, x, c, yzero, VectorOfNoise, ...
    name, mufunc, calling_script, dataDir, plotdir, resolution)

    arguments
        F1, weights0, tspan, x, c, yzero, VectorOfNoise
        name, mufunc, calling_script, dataDir, plotdir
        resolution {mustBeMember(resolution,["normal","low"])} = "normal"
    end

    % Create output folder
    plotsDir = plotdir;
    if ~exist(plotsDir,'dir')
        mkdir(plotsDir);
    end

    % Prepare title
    name = name(1:end-2);
    name = strrep(name,'_',' ');
    firstIteration = true;

    violinmatrix = [];
    iteration    = 20;

    for i = 1:iteration
        distanceVector = [];

        % On last iteration, set up Figure 8
        if i==iteration
            hFig = figure(8);
            clf(hFig);
            set(hFig, ...
                'Name','Noise_Analysis', ...
                'NumberTitle','off', ...
                'Renderer','painters', ...
                'Units','pixels', ...
                'Position',[100 100 1600 900]);

            % compute uniform y-limits from max noise
            [~, Noise_max] = add_Noise_Max(x,VectorOfNoise(end));
            globalAmp = max(abs(Noise_max(:,3)));
            yL = [-globalAmp, globalAmp];
        end

        for noise_scale = VectorOfNoise
            % 1) add noise & fit
            [xNoisy, Noise_norm] = add_Noise_Max(x,noise_scale);
            muN = mufunc(xNoisy);
            opts = optimset('Display','off');
            wN   = lsqcurvefit(F1,weights0,tspan,muN',[],[],opts);
            pred = F1(wN,tspan);

            % 2) compute L2 error
            err = sqrt(mean((pred - yzero).^2));
            distanceVector(end+1,1) = err;

            % 3) on final iteration, draw subplots
            if i==iteration
                idx     = find(VectorOfNoise==noise_scale);
                rowsFig = 4;
                colsFig = numel(VectorOfNoise);

                % subplot 1: noisy state X
                subplot(rowsFig,colsFig,idx)
                plot(x(:,2),xNoisy(:,3),'b-','LineWidth',1);
                % xline(0,'k--'); yline(0,'k--');
                xlabel('$\alpha$','Interpreter','latex')
                ylabel('$X$','Interpreter','latex')
                set(gca,'FontSize',16)

                % subplot 2: noise (fixed y-scale)
                subplot(rowsFig,colsFig,idx+colsFig)
                plot(x(:,2),Noise_norm(:,3),'r-','LineWidth',1);
                % xline(0,'k--'); yline(0,'k--');
                xlabel('$\alpha$','Interpreter','latex')
                ylabel('$Amplitude$','Interpreter','latex')
                ylim(yL)
                set(gca,'FontSize',16)
                title(sprintf('$\\mathbf{N(0,%.1f)}$',noise_scale),'Interpreter','latex')

                % subplot 3: real vs pred
                subplot(rowsFig,colsFig,idx+2*colsFig)
                plot(tspan,yzero,'r-','LineWidth',1); hold on
                % yline(0,'k--','HandleVisibility','off')
                plot(tspan,pred,'b--','LineWidth',1); hold off
                xlabel('Time','Interpreter','latex')
                ylabel('$\alpha$','Interpreter','latex')
                legend('Real \beta','Prediction','Location','NorthWest','FontSize',16)
                set(gca,'FontSize',16)

                % subplot 4: plus intermediate
                subplot(rowsFig,colsFig,idx+3*colsFig)
                plot(tspan,yzero,'r-','LineWidth',1); hold on
                % yline(0,'k--','HandleVisibility','off')
                plot(tspan,pred,'b--','LineWidth',1)
                plot(tspan,muN,'m-','LineWidth',1); hold off
                xlabel('Time','Interpreter','latex')
                ylabel('$\alpha$','Interpreter','latex')
                legend('Real \beta','Prediction','Intermediate','Location','NorthWest','FontSize',14)
                set(gca,'FontSize',16)

                if firstIteration
                    sgtitle(['$\mathbf{' name '}$'],'Interpreter','latex','FontSize',20)
                    firstIteration = false;
                end
            end
        end

        % on last iteration, save Figure 8
        if i==iteration
            figure(hFig);
            figName = get(hFig,'Name');
            safeName = regexprep(figName,'[^\w]','_');
            outFile  = fullfile(plotsDir,[safeName '.pdf']);

            if resolution=="low"
                % rasterize entire figure at 150 dpi inside PDF
                exportgraphics(hFig,outFile, ...
                    'ContentType','image', ...
                    'Resolution',30, ...
                    'BackgroundColor','none');
                % shrink physical size to 6×3.5 inches
                set(hFig, ...
                    'Units','inches', ...
                    'Position',[1 1 6 3.5], ...
                    'PaperUnits','inches', ...
                    'PaperPosition',[0 0 6 3.5]);
                % now rasterize at 40dpi
                exportgraphics(hFig,outFile, ...
                    'ContentType','image', ...
                    'Resolution',40, ...
                    'BackgroundColor','none');

            else
                % full vector output
                exportgraphics(hFig,outFile, ...
                    'ContentType','vector', ...
                    'BackgroundColor','none', ...
                    'Append',false);
            end

            % also plot and save Figure 9
            figure(9);
            plot(VectorOfNoise, distanceVector,'-o','LineWidth',1.5)
            xlabel('Noise variance','Interpreter','latex')
            ylabel('L2 Error','Interpreter','latex')
            title(sprintf('Iteration = %d',iteration))
            set(gca,'FontSize',16)
        end

        violinmatrix = [violinmatrix, distanceVector];
    end

    % save the violinmatrix
    save(fullfile(dataDir,[calling_script '_violinMatrix.mat']),'violinmatrix');

    % summary boxplot in Figure 10
    figure(10);
    labels = arrayfun(@num2str,VectorOfNoise,'Uni',false);
    h = boxplot(violinmatrix','Labels',labels,'Colors','k');
    set(h,{'linew'},{1});
    xlabel('Noise variance','FontSize',16,'Interpreter','latex');
    ylabel('Error IHCV','FontSize',16,'Interpreter','latex');
    title([name ' Bifurcation-Normalized Error'],'FontSize',16);
    set(gca,'FontSize',16);
end

