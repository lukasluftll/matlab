function cube(center, edgeLength, varargin)
% CUBE Plot a 3D cube.
%   CUBE(CENTER, EDGELENGTH) plots N axis-aligned cubes that are centered
%   at the 3D coordinates given by the rows of the Nx3 matrix CENTER.
%
%   CUBE(CENTER, EDGELENGTH, VARARGIN) plots N cubes with the properties 
%   indicated by the name-value pair arguments VARARGIN. 
%   For possible name-value pairs, see the documentation of PATCH.
%
%   Example:
%      cube([0, 5, 3; 4, 1, 1], 1, 'FaceColor', 'red', 'FaceAlpha', 0.5);
%
%   See also PATCH.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check if the user provided enough input arguments.
narginchk(2, inf);

% Make sure the center matrix has 3 columns.
if size(center, 2) ~= 3
    error('CENTER must be a Nx3 matrix.');
end

%% Plot the cubes.
for row = 1 : size(center, 1)
    % Compute the vertices.
    %#ok<*AGROW>
    vertices = center(row,:) - [edgeLength/2, edgeLength/2, edgeLength/2];
    vertices = [vertices; vertices + repmat([edgeLength, 0, 0], 1, 1)]; 
    vertices = [vertices; vertices + repmat([0, edgeLength, 0], 2, 1)];
    vertices = [vertices; vertices + repmat([0, 0, edgeLength], 4, 1)];

    % Define all combinations of the corners to form all 6 faces.
    faces = [1, 2, 4, 3; ...
        1, 2, 6, 5; ...
        2, 4, 8, 6; ...
        4, 3, 7, 8; ...
        7, 5, 1, 3; ...
        5, 6, 8, 7];

    % Plot all faces.
    patch('Faces', faces, 'Vertices', vertices, varargin{:});
end

end
