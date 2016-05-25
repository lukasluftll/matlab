function cube(center, edgeLength, color)
% cube Plot a 3D cube.
%   cube(center, edgeLength) plots an axis-aligned cube.
%
%   cube(center, edgeLength, color) plots a cube with a color indicated by
%   the MATLAB ColorSpec specification. Its faces are aligned with the
%   x, y, and z axis.
%
%   Example
%      cube([0, 5, 3], 1, 'y');
%
% Copyright 2016 Alexander Schaefer

%% Validate input.
% If there is no center given, assume it is the origin.
if (nargin < 1)
    center = zeros(1, 3);
end

% If there is no edge length given, assume unity.
if (nargin < 2)
    edgeLength = 1;
end

% If there is no color specification given, plot green faces.
if (nargin < 3)
    color = 'green';
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
    patch(face(:,1), face(:,2), face(:,3), color);
end

end
