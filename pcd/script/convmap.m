% Check convergence of a selected map value.

%% Fetch parameters.
lidarparams

%% Compute map convergence.
% Print caption.
hline(75, '#')
disp(['Computing ',model,'map convergence for ',dataset,' dataset ...'])

% Load the reflectivity map.
load([resultFolder,'/refmap_',dataset,'.mat'], 'lidarMap', 'h', 'm');

% Find the voxel traversed by the greatest number of rays.
[~,v] = max(h(:) + m(:));

% Store the final map value of the selected voxel.
pfin = lidarMap.data(v);

%% Create lidar map.
parprogress(numel(mappingFile));
gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;
pv = [];
switch lower(model)
    case 'decay'
        % Loop over all scans, compute local maps, and merge them to form 
        % a global map.
        r = zeros(gridsize);
        l = zeros(gridsize);
        for i = 1 : numel(mappingFile)
            % Read laser scan data from file.
            ls = lsread([dataFolder, '/', mappingFile(i).name], rlim);

            % Build the local decay rate map.
            warning('off', 'pcd:mapping:rlim')
            [~,ri,li] = decaymap(ls, xgv, ygv, zgv);
            warning('on', 'pcd:mapping:rlim')
                        
            % Integrate the local map information into the global map.
            r = r + ri.data;
            l = l + li.data;            
            
            % If the scan adds new information to the voxel, store the new
            % voxel value.
            pv(end+1) = r(v) / l(v); %#ok<*SAGROW>

            parprogress;
        end
    case 'ref'
        % Loop over all scans, compute local maps, and merge them to form 
        % a global map.
        h = zeros(gridsize);
        m = zeros(gridsize);
        for i = 1 : numel(mappingFile)
            % Read laser scan data from file.
            ls = lsread([dataFolder, '/', mappingFile(i).name], rlim);
            
            % Build the local reflectivity map.
            warning('off', 'pcd:mapping:rlim')
            [~,hi,mi] = refmap(ls, xgv, ygv, zgv);
            warning('on', 'pcd:mapping:rlim')
            
            % Integrate the local map information into the global map.
            h = h + hi.data;
            m = m + mi.data;            
            pv(end+1) = h(v) / (h(v)+m(v));

            parprogress;
        end
        
        % Compute the global reflectivity map.
        lidarMap = voxelmap(single(h./(h+m)), xgv, ygv, zgv, ...
            htot/(htot+mtot));
    otherwise
        error(['Sensor model ', model, ' not supported.'])
end
parprogress(0);

%% Save map.
save(convFile, 'dataset', 'model', 'dataFolder', 'mappingFile', ...
    'pv', 'pfin', '-v7.3');
display(['Result written to ', lidarMapFile, '.'])
