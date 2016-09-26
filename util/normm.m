function y = normm(x)
% NORMM Normalize matrix.
%    Y = NORMM(X) normalizes the matrix X, so that all elements of Y stay 
%    within [0; 1.0].
%
%    Example:
%       normm(magic(5))
%
%   See also NORMR, NORMC.

% Copyright 2016 Alexander Schaefer

%% Validate input and output.
nargoutchk(0, 1)
narginchk(1, 1)

if iscell(x)
    error('X must not be of type CELL.')
end

%% Normalize input.
y = (x - min(x(:))) / range(x(:));

end
