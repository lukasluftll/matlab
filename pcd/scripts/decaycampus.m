% Build a lidar decay-rate map out of many lidar scans taken on the campus 
% of University of Freiburg.

% Set resolution of the merged point cloud map.
step = 0.1;

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
        [pcd.sensor_qw(:), pcd.sensor_qw(:), pcd.sensor_qy(:), ...
        pcd.sensor_qz(:)]);
    rlim = [2,120];
    ls = laserscan(sp,pcd.azimuth(:),pcd.elevation(:),pcd.radius(:),rlim);
    
    % Merge the point cloud with the map.
    map = pcmerge(map, ls2pc(ls), step);
    
    % Advance the progress bar.
    waitbar(i/numel(file), waitbarHandle);
end

% Close the progress bar.
close(waitbarHandle);
