function noise_5_Boxplots(F1, weights0, tspan, x, c, yzero, VectorOfNoise, name, mufunc, calling_script, dataDir, plotdir)
    % Create the output folder if it doesn't exist
    plotsDir = plotdir;
    if ~exist(plotsDir, 'dir')
        mkdir(plotsDir);
    end

    % Preprocess the name for the suptitle
    name = name(1:end-2);
    name = strrep(name, '_', ' ');
    firstIteration = true;

    % We'll collect error metrics here
    violinmatrix = [];

    % How many iterations?
    iteration = 20;

    % Loop over all iterations
    for i = 1:iteration
        distanceVector = [];
        forindex = 0;

        % If this is the last iteration, initialize figure 8
        if i == iteration
            hFig = figure(8);
            clf(hFig);
            set(hFig, ...
                'Name',        'Noise_Analysis', ...
                'NumberTitle', 'off', ...
                'Renderer',    'painters');
                    % ← insert the pixel-based resize here:
            set(hFig, 'Units','pixels', 'Position',[100 100 1600 900]);


        end

        % For each noise level
        for noise_scale = VectorOfNoise
            forindex = forindex + 1;

            % 1) add noise, compute mu, fit, reproduce
            [xNoisy, Noise_normalized] = add_Noise_Max(x, noise_scale);
            mu_obsNoisy = mufunc(xNoisy);
            opts = optimset('Display', 'off');
            [w_noisy, ~, ~, ~, ~] = lsqcurvefit(F1, weights0, tspan, mu_obsNoisy', [], [], opts);
            reproduced = F1(w_noisy, tspan);

            % 2) compute L2 error
            n = numel(yzero);
            err = sqrt(sum((reproduced - yzero).^2) / n);
            distanceVector(end+1,1) = err;

            % 3) on the last iteration, build subplots in fig 8
            if i == iteration
                rowsFig    = 4;
                colsFig    = numel(VectorOfNoise);
                idx        = forindex;

                % ---- subplot 1: noisy state X ----
                subplot(rowsFig, colsFig, idx)
                plot(x(:,2), xNoisy(:,3), 'b-');
                xline(0,'k--'); yline(0,'k--');
                xlabel('$\alpha$','Interpreter','latex')
                ylabel('$X$','Interpreter','latex')
                set(gca,'FontSize',16)

                % ---- subplot 2: normalized noise ----
                subplot(rowsFig, colsFig, idx + colsFig)
                plot(x(:,2), Noise_normalized(:,3), 'r-');
                xline(0,'k--'); yline(0,'k--');
                xlabel('$\alpha$','Interpreter','latex')
                ylabel('$Amplitude$','Interpreter','latex')
                set(gca,'FontSize',16)
                title(sprintf('$\\mathbf{N~(0,%.1f)}$',noise_scale),'Interpreter','latex')

                % ---- subplot 3: real vs predicted beta ----
                subplot(rowsFig, colsFig, idx + 2*colsFig)
                plot(tspan, yzero, 'r-','LineWidth',1); hold on
                yline(0,'k--','HandleVisibility','off')
                plot(tspan, reproduced, 'b--','LineWidth',1); hold off
                xlabel('Time','Interpreter','latex')
                ylabel('$\alpha$','Interpreter','latex')
                legend('Real \beta','Prediction','Location','NorthWest','FontSize',16)
                set(gca,'FontSize',16)

                % ---- subplot 4: plus intermediate step ----
                subplot(rowsFig, colsFig, idx + 3*colsFig)
                plot(tspan, yzero, 'r-','LineWidth',1); hold on
                yline(0,'k--','HandleVisibility','off')
                plot(tspan, reproduced, 'b--','LineWidth',1)
                plot(tspan, mu_obsNoisy, 'm-'); hold off
                xlabel('Time','Interpreter','latex')
                ylabel('$\alpha$','Interpreter','latex')
                legend('Real \beta','Prediction','Intermediate','Location','NorthWest','FontSize',14)
                set(gca,'FontSize',16)

                % Add a super-title once
                if firstIteration
                    sgtitle(['$\mathbf{' name '}$'], 'Interpreter','latex','FontSize',20)
                    firstIteration = false;
                end
            end
        end

        % Save the fully assembled figure 8 after the last iteration
        if i == iteration
            figure(hFig);
            figName = get(hFig,'Name');
            if isempty(figName)
                safeName = sprintf('Figure_%d', 8);
            else
                safeName = regexprep(figName, '[^\w]', '_');
            end
            outFile = fullfile(plotsDir, [safeName '.pdf']);
            exportgraphics(hFig, outFile, ...
                'ContentType',       'vector', ...
                'BackgroundColor',   'none', ...
                'Append',            false);
        end

        % Collect this iteration's errors
        violinmatrix = [violinmatrix, distanceVector];

        % Plot L2 error vs noise for the last iteration in figure 9
        if i == iteration
            figure(9);
            plot(VectorOfNoise, distanceVector, '-o','LineWidth',1.5)
            xlabel('Noise variance','Interpreter','latex')
            ylabel('L2 Error','Interpreter','latex')
            title(sprintf('Iteration = %d', iteration))
            set(gca,'FontSize',16)
        end
    end

    % Save the violinmatrix to disk
    save(fullfile(dataDir, [calling_script '_violinMatrix.mat']), 'violinmatrix');

    % Finally, create the summary boxplot in figure 10
    figure(10);
    arraystr = arrayfun(@num2str, VectorOfNoise, 'UniformOutput', false);
    h = boxplot(violinmatrix', 'Labels', arraystr, 'Colors', 'k');
    set(h, {'linew'}, {1});
    xlabel('Noise variance','FontSize',16,'Interpreter','latex');
    ylabel('Error IHCV','FontSize',16,'Interpreter','latex');
    title([name ' Bifurcation-Normalized Error'],'FontSize',16);
    set(gca,'FontSize',16);
end
