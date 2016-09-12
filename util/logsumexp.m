function lse = logsumexp(x, varargin)
% LOGSUMEXP Logarithm of sum of exponentials.
%   LOGSUMEXP(X) computes LOG(SUM(EXP(X))) avoiding numerical underflow and 
%   overflow.
%
%   LOGSUMEXP(X,DIM) computes LOG(SUM(EXP(X),DIM)) along the specified
%   dimension of the sum.
%
%   Example:
%      logsumexp([1e3,1e4])
%
%   See also LOG, SUM, EXP.

% Copyright 2016 Alexander Schaefer
% The original algorithm was developed by Tom Minka, then changed by 
% Iain Murray, and finally commented and documented MATLAB-style by 
% Alexander Schaefer.

%% Validate input.
narginchk(1, 2)

% Handle empty input data.
if isempty(x)
    lse = log(sum(exp(x), varargin{:}));
    return;
end

%% Compute output.
% Determine maximum value of input data.
if isempty(varargin)
    xmax = max(x);
else
    xmax = max(x, [], varargin{:});
end

% Subtract maximum from input data.
xshift = bsxfun(@minus, x, xmax);

% Compute logarithm of sum of exponentials and add maximum value again.
lse = bsxfun(@plus, log(sum(exp(xshift), varargin{:})), xmax);

% If the input maximum values contain infinite or NaN values, propagate 
% them to the output.
lse(isinf(xmax)) = xmax(isinf(xmax));
lse(any(isnan(x), varargin{:})) = NaN;

end
