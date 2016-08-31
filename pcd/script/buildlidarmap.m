% Build a map out of many lidar scans.

%% Parameters.
% Define the folder from where to read and where to keep the results.
resultFolder = 'pcd/result';

% Dataset name.
dataset = 'demo';

% Name of the file that contains the merged point cloud.
pcFile = [resultFolder, '/pcmap_', dataset, '.mat'];

% Sensor model to use to build the map: 'decay' | 'ref' | 'lf'.
model = 'decay';

% Resolution of the resulting lidar map.
mapResFine = 1;
mapResCoarse = 5;

% Sensor reading range.
rlim = [2, 120];

% Variance for likelihood field computation.
sigma = 1;

%% Prepare output file.
% Load the file that contains the merged point cloud.
load(pcFile);

% Print caption.
hline(75)
disp(['Computing ', model, 'map for ', dataset, ' dataset ...'])

% Define the name of the output MAT file.
lidarMapFile = [resultFolder, '/', model, 'map_', dataset, '.mat'];

%% Compute extent of lidar map.
% Get the volume where the map resolution is high.
mapLimFine = [pcMap.XLimits; pcMap.YLimits; pcMap.ZLimits];

% Get the PCD file names.
pcdFile = dir([folder, '/*.pcd']);
nPcdFile = numel(pcdFile);

% Iterate over all scans and compute an axis-aligned bounding box of all
% sensor poses.
mapLimCoarse = repmat([Inf, -Inf], 3, 1);
for i = 1 : nPcdFile
    ls = lsread([folder, '/', pcdFile(i).name], rlim);
    spmin = min(tform2trvec(ls.sp)).';
    spmax = max(tform2trvec(ls.sp)).';
    mapLimCoarse = [min([mapLimCoarse(:,1), spmin], [], 2), ...
        max([mapLimCoarse(:,2), spmax], [], 2)];
end

% Compute the bounding box of a scan that consists of maximum range
% readings.
elevationLim = [min(ls.elevation), max(ls.elevation)];
zlim = rlim(2) * sin(elevationLim);
xylim = [-1,+1] * rlim(2);

% Add the extent of the maximum range scan to the bounding box of the
% sensor poses to get the extent of the coarse-grained map.
mapLimCoarse = mapLimCoarse + [xylim; xylim; zlim];

%% Create lidar map.
% Compute the grid vectors of the map.
xgv = [mapLimCoarse(1,1) : mapResCoarse : mapLimFine(1,1), ...
    mapLimFine(1,1) : mapResFine : mapLimFine(1,2), ...
    mapLimFine(1,2) : mapResCoarse : mapLimCoarse(1,2)+mapResCoarse];
ygv = [mapLimCoarse(2,1) : mapResCoarse : mapLimFine(2,1), ...
    mapLimFine(2,1) : mapResFine : mapLimFine(2,2), ...
    mapLimFine(2,2) : mapResCoarse : mapLimCoarse(2,2)+mapResCoarse];
zgv = [mapLimCoarse(3,1) : mapResCoarse : mapLimFine(3,1), ...
    mapLimFine(3,1) : mapResFine : mapLimFine(3,2), ...
    mapLimFine(3,2) : mapResCoarse : mapLimCoarse(3,2)+mapResCoarse];

% Remove too small spacings.
xgv([false, diff(xgv)<mapResFine]) = [];
ygv([false, diff(ygv)<mapResFine]) = [];
zgv([false, diff(zgv)<mapResFine]) = [];

% For the decay model and the reflectivity model, loop over all scans, 
% compute local maps and merge them to form a global map.
% For the endpoint model, simply compute the Gaussian of the distance to
% the nearest obstacle for each voxel.
gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;
parprogress(nPcdFile);
pNan = [];
switch lower(model)
    case 'decay'
        r = zeros(gridsize);
        l = zeros(gridsize);
        rtot = 0;
        ltot = 0;
        parfor i = 1 : nPcdFile
            % Read laser scan data from file.
            ls = lsread([folder, '/', pcdFile(i).name], rlim);

            % Build the local decay rate map.
            warning('off', 'pcd:mapping:rlim')
            [~,ri,li] = decaymap(ls, xgv, ygv, zgv);
            warning('on', 'pcd:mapping:rlim')
                        
            % Integrate the local map information into the global map.
            r = r + ri.data;
            l = l + li.data;
            rtot = rtot + ri.prior;
            ltot = ltot + li.prior;

            % Update progress bar.
            parprogress;
        end
        
        % Compute the decay rate map.
        lidarMap = voxelmap(single(r./l), xgv, ygv, zgv, rtot/ltot);
    case 'ref'
        % Iterate over all laser scans, compute local reflectivity maps and
        % merge them into a global map.
        h = zeros(gridsize);
        m = zeros(gridsize);
        htot = 0;
        mtot = 0;
        parfor i = 1 : nPcdFile
            % Read laser scan data from file.
            ls = lsread([folder, '/', pcdFile(i).name], rlim);
            
            % Build the local reflectivity map.
            warning('off', 'pcd:mapping:rlim')
            [~,hi,mi] = refmap(ls, xgv, ygv, zgv);
            warning('on', 'pcd:mapping:rlim')
            
            % Integrate the local map information into the global map.
            h = h + hi.data;
            m = m + mi.data;
            htot = htot + hi.prior;
            mtot = mtot + mi.prior;
            
            % Update the progress bar.
            parprogress;
        end
        
        % Compute the reflectivity map.
        lidarMap = voxelmap(single(h./(h+m)), xgv, ygv, zgv, ...
            htot/(htot+mtot));
    case 'lf'
        % Compute the likelihood field.
        lidarMap = lfmap(pcMap, sigma, xgv, ygv, zgv);
        lidarMap.data = single(lidarMap.data);
        
        % Count the number of no-returns and the number of returned rays.
        nRet = 0;
        nNret = 0;
        parfor i = 1 : nPcdFile
            % Read laser scan data from file.
            ls = lsread([folder, '/', pcdFile(i).name], rlim);
            
            % Sum up the number of returns and no-returns.
            nRet = nRet + sum(ls.ret);
            nNret = nNret + sum(~ls.ret);
            
            % Update the progress bar.
            parprogress;
        end
        
        % Compute the unconditioned probability of no-returns.
        pNan = nNret / (nRet+nNret);
    otherwise
        error(['Sensor model ', model, ' not supported.'])
end
parprogress(0);

%% Save map.
save(lidarMapFile, 'dataset', 'folder', 'pcRes', 'model', 'mapResFine', ...
    'mapResCoarse', 'rlim', 'sigma', 'lidarMap', 'pNan', '-v7.3');
display(['Result written to ', lidarMapFile, '.'])
