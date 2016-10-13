function y = sq(x)
% SQ Square value.
%   Y = SQ(X) squares X and returns Y.
%
%   Example:
%      sq(34.1)
%
%   See also SQRT, POW.

% Copyright 2016 Alexander Schaefer

%% Validate input.
narginchk(1, 1)
nargoutchk(1, 1)

%% Square input.
y = pow(x, 2);

end
