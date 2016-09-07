% Build fine-grained map out of many lidar scans for visualization.

%% Fetch parameters.
lidarparams
mapRes = 0.2;
xgv = -40 : mapRes : 0;
ygv =  10 : mapRes : 40;
zgv =  -5 : mapRes : 5;

%% Create lidar map.
% Print caption.
hline(75, '#')
disp(['Computing ', model, 'map for ', dataset, ' dataset ...'])

% Compute lidar map.
parprogress(numel(mappingFile));
gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;
pnr = [];
h = [];
m = [];
switch lower(model)
    case 'decay'
        % Loop over all scans, compute local maps, and merge them to form 
        % a global map.
        r = zeros(gridsize);
        l = zeros(gridsize);
        rtot = 0;
        ltot = 0;
        parfor i = 1 : numel(mappingFile)
            % Read laser scan data from file.
            ls = lsread([dataFolder, '/', mappingFile(i).name], rlim);

            % Build the local decay rate map.
            warning('off', 'pcd:mapping:rlim')
            [~,ri,li] = decaymap(ls, xgv, ygv, zgv);
            warning('on', 'pcd:mapping:rlim')
                        
            % Integrate the local map information into the global map.
            r = r + ri.data;
            l = l + li.data;
            rtot = rtot + ri.prior;
            ltot = ltot + li.prior;

            parprogress;
        end
        
        % Compute the global decay rate map.
        lidarMap = voxelmap(single(r./l), xgv, ygv, zgv, rtot/ltot);
    case 'ref'
        % Loop over all scans, compute local maps, and merge them to form 
        % a global map.
        h = zeros(gridsize);
        m = zeros(gridsize);
        htot = 0;
        mtot = 0;
        parfor i = 1 : numel(mappingFile)
            % Read laser scan data from file.
            ls = lsread([dataFolder, '/', mappingFile(i).name], rlim);
            
            % Build the local reflectivity map.
            warning('off', 'pcd:mapping:rlim')
            [~,hi,mi] = refmap(ls, xgv, ygv, zgv);
            warning('on', 'pcd:mapping:rlim')
            
            % Integrate the local map information into the global map.
            h = h + hi.data;
            m = m + mi.data;
            htot = htot + hi.prior;
            mtot = mtot + mi.prior;

            parprogress;
        end
        
        % Compute the global reflectivity map.
        lidarMap = voxelmap(single(h./(h+m)), xgv, ygv, zgv, ...
            htot/(htot+mtot));
    case 'lf'
        % Compute the likelihood field.
        lidarMap = lfmap(pcMap, sigma, xgv, ygv, zgv);
        lidarMap.data = single(lidarMap.data);
        
        % Count the number of no-returns and the number of returned rays.
        nRet = 0;
        nNret = 0;
        parfor i = 1 : numel(mappingFile)
            % Read laser scan data from file.
            ls = lsread([dataFolder, '/', mappingFile(i).name], rlim);
            
            % Sum up the number of returns and no-returns.
            nRet = nRet + sum(ls.ret);
            nNret = nNret + sum(~ls.ret);
            
            parprogress;
        end
        
        % Compute the unconditioned probability of no-returns.
        pnr = nNret / (nRet+nNret);
    otherwise
        error(['Sensor model ', model, ' not supported.'])
end
parprogress(0);

%% Save map.
file = [resultFolder, '/', model, 'map_', dataset, '_vis.mat'];
save(file, 'dataset', 'model', 'dataFolder', 'mappingFile', ...
    'pcMapFile', 'pcRes', 'mapRes', 'rlim', 'sigma', 'lidarMap', 'pnr', ...
    'h', 'm', '-v7.3');
display(['Result written to ', file, '.'])
