function cuboid(limits, varargin)
% CUBOID Plot a 3D cuboid.
%   CUBOID(LIMITS) plots N axis-aligned cuboids with the limits specified 
%   by the 2x3xN matrix LIMITS.
%
%   CUBOID(LIMITS, VARARGIN) plots N cuboids with the properties 
%   indicated by the name-value pair arguments VARARGIN. 
%   For possible name-value pairs, see the documentation of PATCH.
%
%   Example:
%      smallCuboid = [0, 0, 0; 1, 1, 1];
%      largeCuboid = [0, 2, 3; 4, 15, 11];
%      limits = cat(3, smallCuboid, largeCuboid);
%      cuboid(limits, 'FaceColor', 'red', 'FaceAlpha', 0.5);
%
%   See also PATCH, CAT.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check if the user provided enough input arguments.
narginchk(2, inf);

% Make sure the limits matrix has the right size.
if size(limits) < [2, 3, 1]
    error('LIMITS must be a 2x3xN matrix.');
end

%% Plot the cubes.
for n = 1 : size(limits, 3)
    % Compute the vertices.
    vertices = [limits(1,1,n), limits(1,2,n), limits(1,3,n);
        limits(2,1,n), limits(1,2,n), limits(1,3,n);
        limits(2,1,n), limits(2,2,n), limits(1,3,n);
        limits(1,1,n), limits(2,2,n), limits(1,3,n);
        limits(1,1,n), limits(1,2,n), limits(2,3,n);
        limits(2,1,n), limits(1,2,n), limits(2,3,n);
        limits(2,1,n), limits(2,2,n), limits(2,3,n);
        limits(1,1,n), limits(2,2,n), limits(2,3,n)];

    % Define all combinations of the corners to form all 6 faces.
    faces = [1, 2, 3, 4; ...
        1, 2, 6, 5; ...
        2, 3, 7, 6; ...
        3, 4, 8, 7; ...
        4, 1, 5, 8; ...
        5, 6, 7, 8];

    % Plot all faces.
    patch('Faces', faces, 'Vertices', vertices, varargin{:});
end

end
