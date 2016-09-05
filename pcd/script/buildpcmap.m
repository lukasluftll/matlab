% Build a map out of many lidar scans.

%% Fetch parameters.
lidarparams

%% Merge point clouds.
% Print caption.
hline(75, '#')
display(['Merging ', dataset, ' point cloud ...'])

% Assemble the map and denoise it.
pcMap = pcdenoise(pcmap(dataFolder, mappingFile, pcRes));

% Save the point cloud map to file.
save(pcMapFile, 'dataset', 'dataFolder', 'mappingFile', 'pcRes', ...
    'pcMap', '-v7.3')
display(['Result written to ', pcMapFile, '.'])
