function n = nrows(x)
% NROWS Number of rows of X.

% Copyright 2016 Alexander Schaefer

%% Validate input and output.
nargoutchk(0,1)
narginchk(1,1)

%% Compute number of rows.
n = size(x,1);

end
