function plotht(ht, varargin)
% PLOTPOSE Plot 3D homogeneous transformation.
%   PLOTHT(HT) plots a Cartesian coordinate system that represents the
%   rotation and translation contained in the 4x4 homogeneous
%   transformation matrix HT. Given HT is a transformation from system A 
%   to system B, PLOTHT(HT) plots the x, y, and z unit vectors of B in
%   system A as red, green, and blue arrows.
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
% Define the origin and the unit vectors of system B in B in homogeneous 
% coordinates.
ob = [0; 0; 0; 1];
xb = [1; 0; 0; 1];
yb = [0; 1; 0; 1];
zb = [0; 0; 1; 1];

% Transform the origin and the unit vectors of B into system A.
oba = ht * ob;
xba = ht * xb;
yba = ht * yb;
zba = ht * zb;

% Compute the unit vector directions in A.
dx = xba - oba;
dy = yba - oba;
dz = zba - oba;

% Plot the unit vectors of B in A.
quiver3(oba(1), oba(2), oba(3), dx(1), dx(2), dx(3), 0, 'r', varargin{:})
hold on
quiver3(oba(1), oba(2), oba(3), dy(1), dy(2), dy(3), 0, 'g', varargin{:})
quiver3(oba(1), oba(2), oba(3), dz(1), dz(2), dz(3), 0, 'b', varargin{:})
hold off

end
