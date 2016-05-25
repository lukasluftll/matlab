% Segments a field scan into ground and other objects.

% Define parameters.
markerSize = 50;

% Load the downsampled scan of a field.
ws = load('data/field.mat');

% Display the input point cloud.
inputfig = figure('Name', 'Input point cloud', 'NumberTitle', 'Off');
inputax = pcshow(ws.field, 'MarkerSize', markerSize);
xlabel(inputax, 'x [m]');
ylabel(inputax, 'y [m]');
zlabel(inputax, 'z [m]');
whitebg(inputfig, [0.2, 0.2, 0.2]);

% Segment the ground using the cone approach.
coneAngle = 60 * pi/180;
groundIndex = coneseg(ws.field, coneAngle);

% Colorize the points according to the classification result.
green = uint8([34, 139, 34]);
brown = uint8([165, 42, 42]);
ws.field.Color = repmat(green, size(ws.field.Color, 1), 1);
ws.field.Color(groundIndex, :) = repmat(brown, length(groundIndex), 1);

% Display the colorized point cloud.
segfig = figure('Name', 'Segmented field', 'NumberTitle', 'Off');
segax = pcshow(ws.field, 'MarkerSize', markerSize);
xlabel(segax, 'x [m]');
ylabel(segax, 'y [m]');
zlabel(segax, 'z [m]');
whitebg(segfig, [0.5, 0.5, 0.5]);

% Link the cameras of the two figures together.
cameralink = linkprop([inputax, segax], {'CameraUpVector', ...
    'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', cameralink);
