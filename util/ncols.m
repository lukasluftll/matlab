function n = ncols(x)
% NCOLS Number of columns of x.

% Copyright 2016 Alexander Schaefer

%% Validate input and output.
nargoutchk(0,1)
narginchk(1,1)

%% Compute number of columns.
n = size(x,2);
