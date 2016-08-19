% Build a lidar decay-rate map out of many lidar scans taken on the campus 
% of University of Freiburg with the robot Viona.

%% Set paramters.
% Step that determines the fraction of PCD files to use.
step = 1;

% Resolution of the merged point cloud map.
pcres = 0.100;

% Resolution of the resulting decay map.
lambdares = 1.000;

% Sensor reading range.
rlim = [2, 120];

%% Create folder for results.
if ~exist('pcd/results', 'dir')
    mkdir('pcd', 'results');
end

%% Compute extent of decay-rate map.
% Get the PCD file names.
folder = 'pcd/data/campus/pcd_sph';
file = dir([folder, '/*.pcd']);

% Create progress bar.
waitbarHandle = waitbar(0, 'Computing map size ...');

% Iterate over all PCD files and get their maximum extent.
lim = inf(3, 2) .* [+ones(3,1), -ones(3,1)];
for i = 1 : step : numel(file)
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

%% Create one point cloud of all scans.
% Update progress bar label.
waitbar(0, waitbarHandle, 'Building point cloud map ...');

% Iterate over all PCD files.
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

% Update progress bar label.
waitbar(0, waitbarHandle, 'Building decay rate map ...');

% Create the lidar decay rate map.
gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;
r = voxelmap(zeros(gridsize), xgv, ygv, zgv);
l = r.copy;
for i = 1 : step : numel(file)
    % Read laser scan data from file.
    ls = lsread([folder, '/', file(i).name], rlim);
    
    % Build the local decay map.
    [lll,ri,li] = decaymap(ls, xgv, ygv, zgv);
    
    % Integrate the local decay map information into the global map.
    r.add(ri);
    l.add(li);
    
    % Compute the global decay rate map.
    lambda = r ./ l;
    
    % Save the decay rate map to file.
    save('pcd/results/decaycampus.mat', 'lambda');
    
    % Advance the progress bar.
    waitbar(i/numel(file), waitbarHandle);
end

%% Plot decay rate map.
plot(log(lambda));

% Close the progress bar.
close(waitbarHandle);
