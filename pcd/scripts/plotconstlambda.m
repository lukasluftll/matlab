% Load the castle point cloud, divide it into voxels, and visualize the
% lambda value of each voxel.

%% Read data.
% Load the point cloud data.
data = pcdread('data/castle.pcd');

% Create pointCloud object.
cloud = pointCloud(cat(3, data.x, cat(3, data.y, data.z)));

% Define the voxel resolution used to compute the normal distributed ray 
% decay rates.
res = 1;

% Compute the axis-aligned volume spanned by the point cloud.
vol = [min(data.x(:)), min(data.y(:)), min(data.z(:)), ...
    max(data.x(:)), max(data.y(:)), max(data.z(:))];

% Make sure the points in the maximum plane of the volume are part of the
% volume.
vol(4:6) = vol(4:6) + eps(vol(4:6));

% Compute the axis-aligned volume spanned by the point cloud.
vol = [floor(vol(1:3) / res), ceil(vol(4:6) / res)] * res;

%% Visualize decay rate for each voxel.
% Plot the point cloud.
pcshow(cloud, 'MarkerSize', 75);
hold on

% Calculate the ray decay rate.
lambda = raydecay(data.azimuth, data.elevation, data.radius, res, vol);

% Fit lambda into [0; 1].
lambda = lambda / max(lambda(:));

% Visualize the decay rate.
alphaplot(lambda, vol);

% Set the visualization parameters.
axis equal; xlabel('x'); ylabel('y'); zlabel('z'); grid on
hold off
