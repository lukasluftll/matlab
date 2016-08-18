% Build a lidar decay-rate map out of many lidar scans taken on the campus 
% of University of Freiburg with the robot Viona.

%% Set paramters.
% Step that determines the fraction of PCD files to use.
step = 100;

% Resolution of the merged point cloud map.
pcres = 0.150;

% Resolution of the resulting decay map.
lambdares = 0.200;

% Sensor reading range.
rlim = [2, 120];

%% Compute extent decay-rate map.
% Get the PCD file names.
folder = 'pcd/data/campus/pcd_sph';
file = dir([folder, '/*.pcd']);

% Iterate over all PCD files.
waitbarHandle = waitbar(0, 'Computing map size ...');
lim = inf(3, 2) .* [+ones(3,1), -ones(3,1)];
for i = 1 : 100 : numel(file)
    % Read laser scan data from file.
    ls = lsread([folder, '/', file(i).name], rlim);
    
    % Convert spherical data to Cartesian point cloud and denoise cloud.
    pc = pcdenoise(removeInvalidPoints(ls2pc(ls)));
    
    % Store the point cloud's Cartesian limits.
    lim = [min([lim(1,1), pc.XLimits(1)]), max([lim(1,2), pc.XLimits(2)]);
        min([lim(2,1), pc.YLimits(1)]), max([lim(2,2), pc.YLimits(2)]);
        min([lim(3,1), pc.ZLimits(1)]), max([lim(3,2), pc.ZLimits(2)])];
    
    % Advance the progress bar.
    waitbar(i/numel(file), waitbarHandle);
end 

%% Create point cloud of all scans.
% Iterate over all PCD files.
waitbar(0, waitbarHandle, 'Building point cloud map ...');
map = pointCloud(zeros(0, 3));
for i = 1 : step : numel(file)
    % Read laser scan data from file.
    ls = lsread([folder, '/', file(i).name], rlim);
    
    % Merge this point cloud with the map.
    map = pcmerge(map, ls2pc(ls), pcres);
    
    % Advance the progress bar.
    waitbar(i/numel(file), waitbarHandle);
end

%% Create decay rate map.
% Compute the grid vectors of the decay map.
xgv = lim(1,1) : lambdares : lim(1,2)+lambdares;
ygv = lim(2,1) : lambdares : lim(2,2)+lambdares;
zgv = lim(3,1) : lambdares : lim(3,2)+lambdares;

% Create the lidar decay rate map.
waitbar(0, waitbarHandle, 'Building decay rate map ...');
r = voxelmap(zeros([numel(xgv), numel(ygv), numel(zgv)] - 1));
l = voxelmap(r.data);
for i = i : step : numel(file)
    % Read laser scan data from file.
    ls = lsread([folder, '/', file(i).name], rlim);
    
    % Build the local decay map.
    [~,ri,li] = decaymap(ls, xgv, ygv, zgv);
    
    % Merge the local decay map information into the global map.
    r.data = r.data + ri.data;
    l.data = l.data + li.data;
    
    % Advance the progress bar.
    waitbar(i/numel(file), waitbarHandle);
end

% Display the global decay rate map.
lambda = voxelmap(r.data ./ l.data);
plot(lambda);

% Close the progress bar.
close(waitbarHandle);
