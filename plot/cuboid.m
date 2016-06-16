function cuboid(vol, varargin)
% CUBOID Plot a 3D axis-aligned cuboid.
%   CUBOID(VOL) plots an axis-aligned cuboid with the limits specified 
%   by the 6-element vector or matrix VOL.
%
%   VOL is composed of the minimum limits followed by the maximum limits: 
%   VOL(:) = [xmin; ymin; zmin; xmax; ymax; zmax].
%
%   CUBOID(VOL, VARARGIN) plots a cuboid with the properties indicated by 
%   the name-value pair arguments VARARGIN. 
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
if numel(vol) ~= 6
    error('LIMITS must have exactly 6 elements.');
end

%% Calculate the cube vertices and faces.
% Compute the vertices.
vertices = [vol(1), vol(2), vol(3);
    vol(4), vol(2), vol(3);
    vol(4), vol(5), vol(3);
    vol(1), vol(5), vol(3);
    vol(1), vol(2), vol(6);
    vol(4), vol(2), vol(6);
    vol(4), vol(5), vol(6);
    vol(1), vol(5), vol(6)];

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
