function rayplot(azimuth, elevation, radius, ret)
% RAYPLOT Visualize rays of Lidar scanner.
%   RAYPLOT(AZIMUTH, ELEVATION, RADIUS) creates a visualization of the rays
%   originating from a Lidar scanner. Returned rays are plotted in red,
%   no-return rays are plotted in light gray.
%
%   AZIMUTH and ELEVATION are HEIGHTxWIDTH matrices, where HEIGHT and WIDTH 
%   describe the size of the point cloud. The unit is rad.
%
%   RADIUS is a HEIGHTxWIDTH matrix that contains the length of the
%   respective ray. For no-return rays, RADIUS is NaN. 
%
%   RAYPLOT(AZIMUTH, ELEVATION, RADIUS, RET) additionally takes the logical
%   HEIGHTxWIDTH matrix RET as an input. RET specifies whether the a ray 
%   was reflected or not. If not, the corresponding RADIUS tells RAYPLOT 
%   the length of the gray line.
%
%   Example:
%      pc = pcdread('castle.pcd');
%      rayplot(pc.azimuth, pc.elevation, pc.radius)

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(3, 4);

% If RET is not given, create it.
if nargin < 4
    ret = isfinite(radius);
end

% Check whether the spherical coordinate matrices and the reflection matrix
% all have the same number of dimensions and the same size.
if ~(ismatrix(azimuth) && ismatrix(elevation) && ismatrix(radius) ...
        && ismatrix(ret))
    error('AZIMUTH, ELEVATION, RADIUS, and RET must be 2D matrices.')
end
if any(size(azimuth) ~= size(elevation) | size(azimuth) ~= size(radius) ...
        | size(azimuth) ~= size(ret))
    error('AZIMUTH, ELEVATION, RADIUS, and RET must have the same size.')
end

%% Preprocess input data.
% Remove all rows that contain NaN angles.
finiteAngles = isfinite(azimuth) & isfinite(elevation);
azimuth = azimuth(finiteAngles);
elevation = elevation(finiteAngles);
radius = radius(finiteAngles);
ret = ret(finiteAngles);

% Set all NaN radius values to maximum radius.
radius(~ret) = max(radius(:));

% Convert the spherical coordinates to Cartesian coordinates.
[x, y, z] = sph2cart(azimuth, elevation, radius);
p = [x, y, z];

%% Plot rays.
% Plot the returned rays.
pret = kron(p(ret(:),:), [0; 1]);
retray = plot3(pret(:,1), pret(:,2), pret(:,3), 'Color', 'red');
if ~isempty(retray)
    retray.Color(4) = 0.5;
end

% Plot the no-return rays.
pnan = kron(p(~ret(:),:), [0; 1]);
hold on
nanray = plot3(pnan(:,1), pnan(:,2), pnan(:,3), 'Color', 'k');
if ~isempty(nanray)
    nanray.Color(4) = 0.03;
end

%% Plot decoration.
% Plot the sensor origin.
plot3(0, 0, 0, 'Color', 'k', 'Marker', '.', 'MarkerSize', 50);

% Plot the Cartesian axes.
plot3([0, max(x(:))], [0, 0], [0, 0], 'Color', 'r', 'LineWidth', 3);
plot3([0, 0], [0, max(y(:))], [0, 0], 'Color', 'g', 'LineWidth', 3);
plot3([0, 0], [0, 0], [0, max(z(:))], 'Color', 'b', 'LineWidth', 3);

% Plot the point cloud.
pcshow(pointCloud(pret), 'MarkerSize', 30);
hold off

% Label the axes.
xlabel('x')
ylabel('y')
zlabel('z')

axis equal
grid on

end
