function hline(w, c)
% HLINE Horizontal line in command window.
%   HLINE(W) prints a horizontal line of width W in the command window.
%
%   HLINE(W, C) allows to specify a string C that will be repeated as often
%   as needed to make the line W characters long.
%
%   Example:
%      hline(50, '#')
%
%   See also DISPLAY.

% Copyright 2016 Alexander Schaefer

%% Validate input arguments.
% Check the number of input arguments.
narginchk(1, 2)

% Set the default line character.
c(nargin<2) = '=';

% Check the types and values of the input arguments.
if ~isnumeric(w) || rem(w, 1) ~= 0 || w <= 0 
    error('W must be a positive integer.')
end
if ~ischar(c)
    error('C must be a character.')
end

%% Display line.
% Repeat the given character sequence as often as needed to reach the
% specified width.
lineStr = repmat(c, 1, ceil(w/numel(c)));

% Display the line.
display(lineStr(1:w))

end
