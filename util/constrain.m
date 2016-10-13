function y = constrain(x, lim)
% CONSTRAIN Fit value into interval.
%   Y = CONSTRAIN(X, LIM) fits all values of matrix X into the interval
%   defined by the ordered 2-element vector LIM.
%
%   Example:
%      constrain(1:9, [3,6])
%
%   See also MIN, MAX.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check the number of input arguments.
narginchk(2, 2)

% Make sure the interval is an ordered 2-element vector of numeric values.
if ~isnumeric(lim)
    error('LIM must contain numeric values.')
end
if numel(lim) ~= 2
    error('LIM must have exactly 2 elements.')
end
if diff(lim) < 0
    error('LIM must be ordered.')
end

%% Fit input into interval.
x(x < lim(1)) = lim(1);
x(x > lim(2)) = lim(2);
y = x;

end
