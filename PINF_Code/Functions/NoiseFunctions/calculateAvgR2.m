function [avg_r2,r2] = calculateAvgR2(original, predicted)
    % Calculate the mean of the original values
    y_mean = mean(original, 1);
    
    % Calculate the total sum of squares
    ss_tot = sum((original - y_mean).^2, 1);
    
    % Calculate the residual sum of squares
    ss_res = sum((original - predicted).^2, 1);
    
    % Calculate R^2 score for each column
    r2 = 1 - (ss_res ./ ss_tot);
    
    % Calculate the average R^2 score across all columns
    avg_r2 = mean(r2);
end
