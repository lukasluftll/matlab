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
pv = [];
switch lower(model)
    case 'decay'
        % Loop over all scans, compute local maps, and merge them to form 
        % a global map.
        r = zeros(size(lidarMap.data));
        l = zeros(size(lidarMap.data));
        for i = 1 : numel(mappingFile)
            % Read laser scan data from file.
            ls = lsread([dataFolder, '/', mappingFile(i).name], rlim);
            
            % Iterate over each ray of the laser scan.
            for j = 1 : ls.count
                ray = laserscan(ls.sp(:,:,j), ls.azimuth(j), ...
                    ls.elevation(j), ls.radius(j), ls.rlim);

                % Build the local decay rate map.
                warning('off', 'pcd:mapping:rlim')
                [~,ri,li] = decaymap(ray, ...
                    lidarMap.xgv, lidarMap.ygv, lidarMap.zgv);
                warning('on', 'pcd:mapping:rlim')

                % Integrate the local map information into the global map.
                r = r + ri.data;
                l = l + li.data;
                
                % If the ray adds new information to the voxel, store the
                % new voxel value.
                if li.data(v) > 0
                    pv(end+1) = r(v) / l(v); %#ok<*SAGROW>
                end
            end
            
            parprogress;
        end
    case 'ref'
        % Loop over all scans, compute local maps, and merge them to form 
        % a global map.
        h = zeros(size(lidarMap.data));
        m = zeros(size(lidarMap.data));
        for i = 1 : numel(mappingFile)
            % Read laser scan data from file.
            ls = lsread([dataFolder, '/', mappingFile(i).name], rlim);
            
            % Iterate over each ray of the laser scan.
            for j = 1 : ls.count
                ray = laserscan(ls.sp(:,:,j), ls.azimuth(j), ...
                    ls.elevation(j), ls.radius(j), ls.rlim);
                
                % Build the local reflection map.
                warning('off', 'pcd:mapping:rlim')
                [~,hi,mi] = refmap(ls, ...
                    lidarMap.xgv, lidarMap.ygv, lidarMap.zgv);
                warning('on', 'pcd:mapping:rlim')

                % Integrate the local map information into the global map.
                h = h + hi.data;
                m = m + mi.data;
                
                % If the ray adds new information to the voxel, store the
                % new voxel value.
                if h(v)+m(v) > 0
                    pv(end+1) = h(v) / (h(v)+m(v));
                end
            end

            parprogress;
        end
    otherwise
        error(['Sensor model ', model, ' not supported.'])
end
parprogress(0);

%% Save map.
save(convFile, 'dataset', 'model', 'dataFolder', 'mappingFile', ...
    'pv', 'pfin', '-v7.3');
display(['Result written to ', lidarMapFile, '.'])
