function cuboid(limits, varargin)
% CUBOID Plot a 3D axis-aligned cuboid.
%   CUBOID(LIMITS) plots an axis-aligned cuboid with the limits specified 
%   by the 6-element vector or matrix LIMITS.
%
%   If LIMITS is a vector, it is composed of the minimum limits followed by
%   the maximum limits: [xmin, ymin, zmin, xmax, ymax, zmax].
%
%   If LIMITS is a matrix, its columns contain the minimum and maximum 
%   limits, respectively. The rows correspond to the coordinates:
%   [xmin, xmax; ymin, ymax; zmin, zmax].
%
%   CUBOID(LIMITS, VARARGIN) plots a cuboid with the properties 
%   indicated by the name-value pair arguments VARARGIN. 
%   For possible name-value pairs, see the documentation of PATCH.
%
%   Example:
%      cuboid([0, 2; 3, 4; 15, 11], 'FaceColor', 'red', 'FaceAlpha', 0.5);
%
%   See also PATCH, CAT.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check if the user provided the required input argument.
narginchk(1, inf);

% Make sure the limits matrix has the right size.
if numel(limits) ~= 6
    error('LIMITS must have exactly 6 elements.');
end

%% Calculate the cube vertices and faces.
% Compute the vertices.
vertices = [limits(1), limits(2), limits(3);
    limits(4), limits(2), limits(3);
    limits(4), limits(5), limits(3);
    limits(1), limits(5), limits(3);
    limits(1), limits(2), limits(6);
    limits(4), limits(2), limits(6);
    limits(4), limits(5), limits(6);
    limits(1), limits(5), limits(6)];

% Define all combinations of the corners to form all 6 faces.
faces = [1, 2, 3, 4; ...
    1, 2, 6, 5; ...
    2, 3, 7, 6; ...
    3, 4, 8, 7; ...
    4, 1, 5, 8; ...
    5, 6, 7, 8];

%% Plot all faces.
patch('Faces', faces, 'Vertices', vertices, varargin{:});

end
