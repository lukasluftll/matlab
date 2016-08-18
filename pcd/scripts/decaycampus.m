% Build a lidar decay-rate map out of many lidar scans taken on the campus 
% of University of Freiburg.

% Resolution of the merged point cloud map.
pcres = 0.150;

% Resolution of the resulting decay map.
lambdares = 0.200;

% Get the PCD file names.
folder = '~/ros/datasets/lifenav/2013-10-31-campus/pcd_sph';
file = dir([folder, '/*.pcd']);

% Create the map point cloud.
map = pointCloud(zeros(0, 3));

% Create a progress bar.
waitbarHandle = waitbar(0, 'Building map ...');

% Loop over all PCD files.
for i = 1 : 100 : numel(file)
    % Read the PCD file.
    pcd = pcdread([folder, '/', file(i).name]);
    
    % Compute the Cartesian coordinates of the ray endpoints with reference
    % to the world frame.
    sp = trquat2tform(...
        [pcd.sensor_x(:), pcd.sensor_y(:), pcd.sensor_z(:)], ...
        [pcd.sensor_qw(:), pcd.sensor_qx(:), pcd.sensor_qy(:), ...
        pcd.sensor_qz(:)]);
    rlim = [2,120];
    ls = laserscan(sp,pcd.azimuth(:),pcd.elevation(:),pcd.radius(:),rlim);
    
    % Merge the point cloud with the map.
    map = pcmerge(map, ls2pc(ls), pcres);
    
    % Advance the progress bar.
    waitbar(i/numel(file), waitbarHandle);
end

waitbar(0, waitbarHandle, 'Building decay rate map ...');

% Compute the grid vectors of the decay map.
xgv = map.XLimits(1) : lambdares : map.XLimits(2);
ygv = map.YLimits(1) : lambdares : map.YLimits(2);
zgv = map.ZLimits(1) : lambdares : map.ZLimits(2);

% Create the lidar decay rate map.
r = voxelmap(zeros([numel(xgv), numel(ygv), numel(zgv)] - 1));
l = voxelmap(r.data);
for i = i : 100 : numel(file)
    pcd = pcdread([folder, '/', file(i).name]);
    
    % Compute the Cartesian coordinates of the ray endpoints with reference
    % to the world frame.
    sp = trquat2tform(...
        [pcd.sensor_x(:), pcd.sensor_y(:), pcd.sensor_z(:)], ...
        [pcd.sensor_qw(:), pcd.sensor_qx(:), pcd.sensor_qy(:), ...
        pcd.sensor_qz(:)]);
    rlim = [2,120];
    ls = laserscan(sp,pcd.azimuth(:),pcd.elevation(:),pcd.radius(:),rlim);
    
    [~,ri,li] = decaymap(ls, xgv, ygv, zgv);
    r.data = r.data + ri.data;
    l.data = l.data + li.data;
    
    % Advance the progress bar.
    waitbar(i/numel(file), waitbarHandle);
end

% Close the progress bar.
close(waitbarHandle);

lambda = voxelmap(r.data ./ l.data);
plot(lambda);
