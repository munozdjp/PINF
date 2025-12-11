function [R2Score, r2perVariable, dxdt_fd, dxdt_model] = function_CalcDerivatives(model, params, xModel, tspan)
    % modelFunc: Function handle for the model, e.g., @hopfPolyOrder3
    % params: Structure containing all necessary parameters for the model
    % xNoisy: Noisy trajectory data (160x4 matrix)
    % tspan: Time span for the integration
    % Unpack parameters from the cell array
%     unpackedParams = num2cell(params{:});
     % Calculate model-derived derivatives for each time step
    dxdt_model = zeros(size(xModel));
    for idx = 1:length(tspan)
        dxdt_model(idx, :) = model(tspan(idx), xModel(idx, :), params{:});
    end

    % Calculate finite difference derivatives
    dxdt_fd = diff(xModel) ./ diff(tspan');
    dxdt_fd = [dxdt_fd; dxdt_fd(end, :)]; % Append the last row to match dimensions

    % Calculate R^2 score between model-based and finite difference derivatives
    addpath('/Users/munozdjp/Library/CloudStorage/GoogleDrive-juan.munozdiaz@kaust.edu.sa/My Drive/PHD_Kaust_DriveFolder/Academic_KAUST/ACDC_project/projects1/HINNDy project 2022/NoiseFunctions'); % Ensure calculateAvgR2.m is accessible
    [R2Score, r2perVariable]= calculateAvgR2(dxdt_fd(:,2:end), dxdt_model(:,2:end));
    %display mean square error
    figure;
    nVars = size(dxdt_model, 2); % Assuming dxdt_model and dxdt_fd have the same number of columns

    % Colors and markers for differentiation
    colors = ['b', 'r', 'g', 'm']; % Example colors for up to 4 variables
    markers = ['o', '*', '+', 'x']; % Example markers for up to 4 variables

    for i = 1:nVars
        % Plot model-based derivatives
        plot(tspan, dxdt_model(:,i), ['-', colors(i)], 'LineWidth', 1.5, 'DisplayName', ['Model dx', num2str(i), '/dt']);
        hold on; % Keep the plot active to overlay the next set of derivatives

        % Plot finite difference derivatives
        plot(tspan, dxdt_fd(:,i), ['--', colors(i), markers(i)], 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', ['FD dx', num2str(i), '/dt']);
    end

    title('Comparison of Derivatives');
    xlabel('Time');
    ylabel('Derivative Value');
    legend('Location', 'best');
    hold off; % Release the plot for further commands

%     figure;
% 
%     % Plot model-based derivatives with continuous lines
%     plot(tspan, dxdt_model,'--', 'LineWidth', 1);
%     hold on; % Keep the plot active to overlay the next set of derivatives
% 
%     % Plot finite difference derivatives with star markers for distinction
%     plot(tspan, dxdt_fd, '*--', 'LineWidth', 1, 'MarkerSize', 2);
% 
%     title('Comparison of Derivatives');
%     xlabel('Time');
%     ylabel('Derivative Value');
% 
%     % Adjust the legend to clearly differentiate between model-based and finite difference derivatives
%     legend({'Model dx1/dt', 'Model dx2/dt', 'Model dx3/dt', 'Model dx4/dt', ...
%             'FD dx1/dt', 'FD dx2/dt', 'FD dx3/dt', 'FD dx4/dt'}, ...
%             'Location', 'best');
% 
%     hold off; % Release the plot for further commands


    

    % Plotting
%     figure;
%     subplot(2,1,1);
%     plot(tspan, dxdt_model, 'LineWidth', 2);
%     title('Model-based Derivatives');
%     xlabel('Time');
%     ylabel('Derivative Value');
%     legend('dx1/dt', 'dx2/dt', 'dx3/dt', 'dx4/dt');
% 
%     subplot(2,1,2);
%     plot(tspan, dxdt_fd, 'LineWidth', 2);
%     title('Finite Difference Derivatives');
%     xlabel('Time');
%     ylabel('Derivative Value');
%     legend('dx1/dt', 'dx2/dt', 'dx3/dt', 'dx4/dt');

    % The function now returns R2Score, dxdt_fd, and dxdt_model
end
