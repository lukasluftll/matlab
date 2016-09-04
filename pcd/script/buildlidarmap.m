% Build a map out of many lidar scans.

%% Fetch parameters.
lidarparams

%% Prepare output file.
% Load the file that contains the merged point cloud.
load(pcMapFile, 'pcMap');

% Print caption.
hline(75, '#')
disp(['Computing ', model, 'map for ', dataset, ' dataset ...'])

%% Compute extent of lidar map.
% Iterate over all scans and compute an axis-aligned bounding box of all
% ray endpoints.
maplim = repmat([Inf, -Inf], 3, 1, numel(mappingFile));
disp('Determining map size ...')
parprogress(numel(mappingFile));
parfor i = 1 : numel(mappingFile)
    ls = lsread([dataFolder,'/',mappingFile(i).name], rlim); %#ok<*PFGV>
    
    % Set the length of all rays to maximum sensor range.
    ls.radius = repmat(ls.rlim(2), size(ls.radius));
    
    % Compute the minimum and maximum coordinates of the endpoints of this
    % scan.
    pc = ls2pc(ls);
    maplim(:,:,i)= [pc.XLimits; pc.YLimits; pc.ZLimits];
    
    parprogress;
end
parprogress(0);

% Compute the global minimum and maximum coordinates of the scan endpoints.
maplim = [min(maplim(:,1,:), [], 3), max(maplim(:,2,:), [], 3)];

% Compute the grid vectors of the map.
xgv = maplim(1,1)-mapRes-rkli : mapRes : maplim(1,2)+mapRes+rkli;
ygv = maplim(2,1)-mapRes-rkli : mapRes : maplim(2,2)+mapRes+rkli;
zgv = maplim(3,1)-mapRes : mapRes : maplim(3,2)+mapRes;

%% Create lidar map.
% For the decay model and the reflectivity model, loop over all scans, 
% compute local maps and merge them to form a global map.
% For the endpoint model, simply compute the Gaussian of the distance to
% the nearest obstacle for each voxel.
gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;
disp('Computing map ...')
parprogress(numel(mappingFile));
pnr = [];
switch lower(model)
    case 'decay'
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
        parfor i = 1 : numel(mappingFile)
            % Read laser scan data from file.
            ls = lsread([dataFolder, '/', mappingFile(i).name], rlim);
            
            % Sum up the number of returns and no-returns.
            nRet = nRet + sum(ls.ret);
            nNret = nNret + sum(~ls.ret);
            
            % Update the progress bar.
            parprogress;
        end
        
        % Compute the unconditioned probability of no-returns.
        pnr = nNret / (nRet+nNret);
    otherwise
        error(['Sensor model ', model, ' not supported.'])
end
parprogress(0);

%% Save map.
save(lidarMapFile, 'dataset', 'model', 'pcMapFile', 'dataFolder', ...
    'pcRes', 'mapRes', 'rlim', 'sigma', 'lidarMap', 'pnr', '-v7.3');
display(['Result written to ', lidarMapFile, '.'])
