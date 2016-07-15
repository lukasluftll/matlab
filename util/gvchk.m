function gvchk(varargin)
% GVCHK Checks validity of grid vectors.
%
%   GVCHK(GVX, GVY, ...) checks whether each of the given grid vectors
%   contains at least two elements and whether it increases monotonically.
%   If not, it displays an error message and stops the calling function.

% Copyright 2016 Alexander Schaefer

% Check all input arguments.
for i = 1 : numel(varargin)
    % Check if the grid vector has at least 2 elements.
    if numel(varargin{i}) < 2
        error('Every grid vector must contain at least 2 elements.')
    end
    
    % Check whether the grid vector is ordered.
    if any(diff(varargin{i}(:)) <= 0)
        error('Grid vectors must monotonically increase.')
    end
end

end
