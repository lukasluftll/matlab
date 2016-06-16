%% Read data.
% Load the point cloud data.
data = pcdread('data/castle.pcd');

% Create pointCloud object.
cloud = pointCloud(cat(3, data.x, cat(3, data.y, data.z)));

% Define the voxel resolution used to compute the normal distributed ray 
% decay rates.
res = 1;

% Compute the axis-aligned volume spanned by the point cloud.
vol = [floor([min(data.x(:)), min(data.y(:)), min(data.z(:))]/res), ...
        ceil([max(data.x(:)), max(data.y(:)), max(data.z(:))]/res)] * res;

% Make sure the points in the maximum plane of the volume are part of the
% volume.
vol(4:6) = vol(4:6) + eps(vol(4:6));

%% Visualize decay rate for each voxel.
% Calculate the ray decay rate.
lambda = raydecay(data.azimuth, data.elevation, data.radius, ...
    res, vol);

% Visualize the decay rate.
alphaplot(lambda, vol);
axis equal
xlabel('x'); ylabel('y'); zlabel('z')

% Overlay the point cloud.
hold on
pcshow(cloud);
hold off
