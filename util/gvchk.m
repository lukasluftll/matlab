function gvchk(varargin)
% GVCHK Checks validity of grid vectors.
%
%   GVCHK(GVX, GVY, ...) checks whether each of the given grid vectors
%   contains at least two elements, whether their values are finite and 
%   whether they increase monotonically.
%   If not, it displays an error message and stops the calling function.

% Copyright 2016 Alexander Schaefer

% Check all input arguments.
for i = 1 : numel(varargin)
    % Check if the grid vector has at least 2 elements.
    if numel(varargin{i}) < 2
        error('Every grid vector must contain at least 2 elements.')
    end
    
    % Check if the grid vector values are finite.
    if any(~isfinite(varargin{i}(:)))
        error('Grid vector elements must not be NaN or Inf.')
    end
    
    % Check whether the grid vector is ordered.
    if any(diff(varargin{i}(:)) <= 0)
        error('Grid vectors must monotonically increase.')
    end
end

end
