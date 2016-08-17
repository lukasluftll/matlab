% Build a lidar decay-rate map out of many lidar scans taken on the campus 
% of University of Freiburg.

% Get the PCD file names.
file = dir(['data/pcd_sph', '/*.pcd']);

% Create the map point cloud.
map = pointCloud(zeros(0, 3));

% Create a progress bar and set up automatic destruction after use.
waitbarHandle = waitbar(0, 'Building map ...');
cleanupObj = onCleanup(@() close(waitbarHandle));

% Loop over all PCD files.
for i = 1 : numel(file)
    % Read the PCD file.
    pcd = pcdread([folder, '/', file(i).name]);
    pc = pointCloud([pcd.x(:), pcd.y(:), pcd.z(:)]);
    
    % If odometry information is given, transform the point cloud
    % accordingly.
    if isfield(pcd, 'odometry')
        pc = pctransform(pc, ht2affine3d(pcd.odometry));    
    end

    % Merge the point cloud with the map.
    map = pcmerge(map, pc, step);
    
    % Advance the progress bar.
    waitbar(i/numel(file), waitbarHandle);
end