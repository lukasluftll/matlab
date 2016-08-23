% Build a map out of many lidar scans.

%% Parameters.
% Define the folder where to keep the results.
resultFolder = 'pcd/result';

% Name of the file that contains the merged point cloud.
pcFile = [resultFolder, '/pcmap_campus.mat'];

% Sensor model to use to build map: 'decay' | 'ref'.
model = 'decay';

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
% Load the file that contains the merged point cloud.
load(pcFile);

% Define the name of the output MAT file.
lidarMapFile = [resultFolder, '/', model, 'map_', dataset, '.mat'];

% Save parameters to file.
save(lidarMapFile, 'dataset', 'folder', 'pcMapRes', 'model', 'lidarMapRes', ...
    'rlim', '-v7.3');

%% Compute map extent.
% Denoise the map.
pcMap = pcdenoise(pcMap);

% Compute the extent of the denoised point cloud.
lim = [pcMap.XLimits(1), pcMap.XLimits(2);
    pcMap.YLimits(1), pcMap.YLimits(2);
    pcMap.ZLimits(1), pcMap.ZLimits(2)];

%% Create lidar map.
% Get the PCD file names.
pcdFile = dir([folder, '/*.pcd']);

% Compute the grid vectors of the map.
xgv = lim(1,1) : lidarMapRes : lim(1,2)+lidarMapRes;
ygv = lim(2,1) : lidarMapRes : lim(2,2)+lidarMapRes;
zgv = lim(3,1) : lidarMapRes : lim(3,2)+lidarMapRes;

% Create the progress bar.
progressbar(['Computing ', model, ' map ...']);
  
% Iterate over all laser scans, compute local maps and merge them into a 
% global map.
gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;
num = zeros(gridsize);
denom = zeros(gridsize);
nPcdFile = numel(pcdFile);
parfor i = 1 : nPcdFile
    % Read laser scan data from file.
    ls = lsread([folder, '/', pcdFile(i).name], rlim);

    % Build the local lidar map.
    [~,numi,denomi] = mapFun(ls, xgv, ygv, zgv);

    % Integrate the local map information into the global map.
    num = num + numi;
    denom = denom + denomi;
    
    % Display the progress of the first worker.
    task = getCurrentTask;
    if task.ID == 1
        progressbar(i/nPcdFile);
    end
end

%% Save the resulting lidar map.
lidarMap = voxelmap(num./denom, xgv, ygv, zgv);
save(lidarMapFile, 'lidarMap', '-append');
