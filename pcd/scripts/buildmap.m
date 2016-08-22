% Build a map out of many lidar scans.

%% Parameters.
% Dataset name.
dataset = 'campus';

% Sensor model to use to build map: 'decay' | 'ref'.
model = 'decay';

% Dataset folder with PCD files.
folder = ['pcd/data/', dataset, '/pcd_sph'];

% Step that determines the fraction of PCD files to use.
step = 1000;

% Resolution of the merged point cloud map.
pcMapRes = 0.1;

% Resolution of the resulting lidar map.
lidarMapRes = 0.5;

% Sensor reading range.
rlim = [2, 120];

%% Select mapping algorithm.
switch model
    case 'decay'
        mapFun = @decaymap;
    case 'ref'
        mapFun = @refmap;
    otherwise
        error(['Map model ', model, ' not supported.'])
end

%% Prepare output file.
% Create folder for results.
if ~exist('pcd/results', 'dir')
    mkdir('pcd', 'results');
end

% Define the name of the output MAT file.
outfile = ['pcd/results/', model, 'map_', dataset, '.mat'];

% Save parameters to file.
save(outfile, 'dataset', 'model', 'folder', 'step', 'pcMapRes', ...
    'lidarMapRes', 'rlim', '-v7.3');

%% Merge point clouds.
% Get the PCD file names.
infile = dir([folder, '/*.pcd']);

% Create progress bar.
waitbarHandle = waitbar(0, 'Merging point cloud map ...');

% Iterate over all PCD files.
pcMap = pointCloud(zeros(0, 3));
for i = 1 : step : numel(infile)
    % Read laser scan data from file.
    ls = lsread([folder, '/', infile(i).name], rlim);
    
    % Merge this point cloud with the map.
    pcMap = pcmerge(pcMap, ls2pc(ls), pcMapRes);
    
    % Advance the progress bar.
    waitbar(i/numel(infile), waitbarHandle);
end

% Save the point cloud map to file.
save(outfile, 'pcMap', '-append');

% Denoise the map.
pcMap = pcdenoise(pcMap);

% Compute the extent of the denoised point cloud.
lim = [pcMap.XLimits(1), pcMap.XLimits(2);
    pcMap.YLimits(1), pcMap.YLimits(2);
    pcMap.ZLimits(1), pcMap.ZLimits(2)];

%% Create lidar map.
% Compute the grid vectors of the map.
xgv = lim(1,1) : lidarMapRes : lim(1,2)+lidarMapRes;
ygv = lim(2,1) : lidarMapRes : lim(2,2)+lidarMapRes;
zgv = lim(3,1) : lidarMapRes : lim(3,2)+lidarMapRes;

% Update the label of the progress bar.
waitbar(0, waitbarHandle, ['Computing ', model, ' map ...']);

% Use multiple workers to create the lidar map.
gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;
iInfile = 1 : step : numel(infile);
spmd
    % For each worker, preallocate the result matrices.
    num = voxelmap(zeros(gridsize), xgv, ygv, zgv);
    denom = num.copy;
    
    % Iterate over the worker's share of all laser scans.
    for i = iInfile(labindex : numlabs : end)
        % Read laser scan data from file.
        ls = lsread([folder, '/', infile(i).name], rlim);

        % Build the local lidar map.
        warning off
        [~,ai,bi] = mapFun(ls, xgv, ygv, zgv);
        warning on

        % Integrate the local map information into the global map.
        num.add(ai);
        denom.add(bi);

        % Advance the progress bar.
        if labindex == 1
            %waitbar(i/numel(infile), waitbarHandle);
        end
    end
end

%% Save the resulting lidar map.
lidarMap = sum([num{:}]) ./ sum([denom{:}]);
save(outfile, 'lidarMap', '-append');
        
% Close the progress bar.
close(waitbarHandle);
