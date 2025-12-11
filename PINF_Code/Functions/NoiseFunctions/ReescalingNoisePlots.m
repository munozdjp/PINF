% Get the minimum and maximum values from your matrix
min_val = min(violinmatrix(:));
max_val = max(violinmatrix(:));

% Subtract the minimum value from all values in the matrix
violinmatrix = violinmatrix - min_val;

% Divide all values in the matrix by the range
violinmatrix = violinmatrix / (max_val - min_val);
