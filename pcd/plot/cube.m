function cube(center, edgeLength, varargin)
% CUBE Plot a 3D cube.
%   CUBE(CENTER, EDGELENGTH) plots an axis-aligned cube.
%
%   CUBE(CENTER, EDGELENGTH, VARARGIN) plots a cube with the properties 
%   indicated by a color indicated by the name-value pair arguments
%   varargin. For possible name-value pairs, see the documentation of 
%   function PATCH.
%
%   Example
%      cube([0, 5, 3], 1, 'FaceColor', 'yellow');
%
% Copyright 2016 Alexander Schaefer

% TODO
% + Specify face color.
% - Vectorize input.
% - Use face-vertex arguments of patch.

%% Validate input.
% Check if the user provided enough input arguments.
if nargin < 2
    error('Not enough input arguments.');
end

%% Plot the cube.
% Compute the corner coordinates.
corner = center - [edgeLength/2, edgeLength/2, edgeLength/2];
corner = [corner; corner + repmat([edgeLength, 0, 0], 1, 1)];
corner = [corner; corner + repmat([0, edgeLength, 0], 2, 1)];
corner = [corner; corner + repmat([0, 0, edgeLength], 4, 1)];

% Define all combinations of the corners to form all 6 faces.
combination = [1, 2, 4, 3; ...
    1, 2, 6, 5; ...
    2, 4, 8, 6; ...
    4, 3, 7, 8; ...
    7, 5, 1, 3; ...
    5, 6, 8, 7];

% Plot all faces.
for c = 1 : size(combination, 1)
    face = corner(combination(c,:), :);
    patch(face(:,1), face(:,2), face(:,3), varargin{:});
end

end
