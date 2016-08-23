% Visualizes the result of performing an NDT on the point cloud depicting
% the Lidar scan of of castle.

%% Define parameters.
% Define the resolution of the grid used for NDT.
res = 1.5;

% Define the radius used for NDT.
radius = 5;

%% Perform NDT.
% Read the point cloud.
data = pcdread('data/castle.pcd');
cloud = pointCloud(cat(3, data.x, cat(3, data.y, data.z)));

% Perform NDT on the grid.
xgv = cloud.XLimits(1)-0.5 : res : cloud.XLimits(2)+0.5;
ygv = cloud.YLimits(1)-0.5 : res : cloud.YLimits(2)+0.5;
zgv = cloud.ZLimits(1)-0.5 : res : cloud.ZLimits(2)+0.5;
[x, y, z] = meshgrid(xgv, ygv, zgv);
center = [x(:), y(:), z(:)];
[mu, sigma] = ndt(cloud, center, radius);

% Sum up the probabilty densities of all NDT distributions for all points 
% of the grid.
density = reshape(ndpdf(center, mu, sigma), size(x));

%% Plot NDT.
% Plot a surface of constant density.
ndtfig = figure('Name', 'Isosurface', 'NumberTitle', 'Off');
isosurface(xgv, ygv, zgv, density, 0.1);
ndtax = gca;
axis equal

% Plot the point cloud.
cloudfig = figure('Name', 'Point cloud', 'NumberTitle', 'Off');
cloudax = pcshow(cloud);

% Link the cameras of the two figures together.
cameralink = linkprop([ndtax, cloudax], {'CameraUpVector', ...
    'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', cameralink);
