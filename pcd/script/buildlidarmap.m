% Build a map out of many lidar scans.

%% Parameters.
% Define the folder from where to read and where to keep the results.
resultFolder = 'pcd/result';

% Name of the file that contains the merged point cloud.
pcFile = [resultFolder, '/pcmap_demo.mat'];

% Sensor model to use to build map: 'decay' | 'ref' | 'lf'.
model = 'ref';

% Resolution of the resulting lidar map.
lidarMapRes = 0.5;

% Sensor reading range.
rlim = [2, 120];

% Variance for likelihood field computation.
sigma = 1;

%% Prepare output file.
% Load the file that contains the merged point cloud.
load(pcFile);

% Print caption.
hline(75)
disp(['Computing ', model, 'map for ', dataset, ' ...'])

% Define the name of the output MAT file.
lidarMapFile = [resultFolder, '/', model, 'map_', dataset, '.mat'];

%% Compute map extent.
% Get the PCD file names.
pcdFile = dir([folder, '/*.pcd']);
nPcdFile = numel(pcdFile);

% Iterate over all scans and compute an axis-aligned bounding box of all
% sensor poses.
maplim = repmat([Inf, -Inf], 3, 1);
for i = 1 : nPcdFile
    ls = lsread([folder, '/', pcdFile(i).name], rlim);
    lsmin = min(tform2trvec(ls.sp)).';
    lsmax = max(tform2trvec(ls.sp)).';
    maplim = [min([maplim(:,1), lsmin], [], 2), ...
        max([maplim(:,2), lsmax], [], 2)];
end

% Compute the bounding box of a scan that consists of maximum range
% readings only.
nScan4Rev = numel(pcdFile);
lsmr = laserscan.empty(0, nScan4Rev);
for i = 1 : nScan4Rev
    lsmr(i) = lsread([folder, '/', pcdFile(i).name], rlim);
end
%%% TODO: Check point clouds of laser scans.
lsmr = lsconcat(lsmr); 
lsmr.radius = repmat(lsmr.rlim(2), size(lsmr.radius));
pcmr = ls2pc(lsmr);

% Add the extent of the maximum range point cloud to the bounding box of
% all sensor poses to ensure the map covers all rays.
maplim = maplim + [pcmr.XLimits; pcmr.YLimits; pcmr.ZLimits];

%% Create lidar map.
% Define the mapping function to use.
switch lower(model)
    case 'decay'
        mapFun = @decaymap;
    case 'ref'
        mapFun = @refmap;
    case 'lf'
        mapFun = @lfmap;
    otherwise
        error('Mapping function not supported.')
end
    
% Compute the grid vectors of the map.
xgv = maplim(1,1) : lidarMapRes : maplim(1,2);
ygv = maplim(2,1) : lidarMapRes : maplim(2,2);
zgv = maplim(3,1) : lidarMapRes : maplim(3,2);

% For the decay and the reflectivity model, loop over all scans and compute
% the numerator and denominator of each map cell's value.
if strcmpi(model, 'decay') || strcmpi(model, 'ref')    
    % Iterate over all laser scans, compute local maps and merge them into 
    % a global map.
    gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;
    num = zeros(gridsize);
    denom = zeros(gridsize);
    parprogress(nPcdFile);
    parfor i = 1 : nPcdFile
        % Read laser scan data from file.
        ls = lsread([folder, '/', pcdFile(i).name], rlim);

        % Build the local lidar map.
        warning('off') %#ok<WNOFF>
        [~,numi,denomi] = mapFun(ls, xgv, ygv, zgv);
        warning('on') %#ok<WNON>

        % Integrate the local map information into the global map.
        num = num + numi;
        denom = denom + denomi;

        % Update progress bar.
        parprogress;
    end
    parprogress(0);
    
    % Compute the lidar map.
    lidarMap = voxelmap(num./denom, xgv, ygv, zgv);
else
    % Compute the likelihood field.
    lidarMap = mapFun(pcMap, sigma, xgv, ygv, zgv);
end

%% Save map.
save(lidarMapFile, 'dataset', 'folder', 'pcMapRes', 'model', ...
    'lidarMapRes', 'rlim', 'sigma', 'lidarMap', '-v7.3');
display(['Result written to ', lidarMapFile, '.'])
