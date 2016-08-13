function plotht(ht, l, varargin)
% PLOTPOSE Plot 3D homogeneous transformation.
%   PLOTHT(HT) plots a Cartesian coordinate system that represents the
%   rotation and translation contained in the 4x4 homogeneous
%   transformation matrix HT. Given HT is a transformation from system A 
%   to system B, PLOTHT(HT) plots the x, y, and z unit vectors of B in
%   system A as red, green, and blue arrows.
%
%   PLOTHT(HT, L) plots the axis vectors with length L.
%
%   PLOTHT(HT, VARARGIN) additionally specifies the plot linestyle for the
%   axis vectors. Any marker in LINESPEC is drawn at the base instead of an
%   arrow on the tip. Use a marker of '.' to specify no marker at all. 
%   See PLOT for available linestyles.
%
%   Example:
%      plotht(eye(4))
%
%   See also PLOT, QUIVER3.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(1, +Inf)

% Check HT is a homogeneous transformation matrix.
if ~ishrt(ht)
    error('HT must be a 4x4 homogeneous translation-rotation matrix.')
end

%% Plot unit vectors.
% Define the origin of B in A.
o = ht(1:3,4);

% Compute the unit vectors of B in A.
e = ht(1:3,1:3) * l * eye(3);

% Plot the unit vectors of B in A.
quiver3(o(1), o(2), o(3), e(1,1), e(2,1), e(3,1), 0, 'r', varargin{:})
hold on
quiver3(o(1), o(2), o(3), e(1,2), e(2,2), e(3,2), 0, 'g', varargin{:})
quiver3(o(1), o(2), o(3), e(1,3), e(2,3), e(3,3), 0, 'b', varargin{:})
hold off

end
