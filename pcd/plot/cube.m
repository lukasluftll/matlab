function cube(origin, edgeLength, color)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if (nargin < 3)
    color = 'green';
end

corners = origin - [edgeLength/2, edgeLength/2, edgeLength/2];
corners = [corners; corners + repmat([edgeLength, 0, 0], 1, 1)];
corners = [corners; corners + repmat([0, edgeLength, 0], 2, 1)];
corners = [corners; corners + repmat([0, 0, edgeLength], 4, 1)];

combinations = [1, 2, 4, 3; ...
    1, 2, 6, 5; ...
    2, 4, 8, 6; ...
    4, 3, 7, 8; ...
    7, 5, 1, 3; ...
    5, 6, 8, 7];

for (c = 1 : size(combinations, 1))
    face = corners(combinations(c,:), :);
    patch(face(:,1), face(:,2), face(:,3), color);
end

end
