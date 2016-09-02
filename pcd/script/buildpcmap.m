% Build a map out of many lidar scans.

%% Fetch parameters.
lidarparams

%% Prepare output file.
% Print caption.
hline(75, '#')
display(['Merging ', dataset, ' point cloud ...'])

%% Merge point clouds.
% Build the map and denoise it.
pcMap = pcdenoise(pcmap(mappingFolder, pcRes));

% Save the point cloud map to file.
save(pcMapFile, 'dataset', 'mappingFolder', 'pcRes', 'pcMap', '-v7.3')
display(['Result written to ', pcMapFile, '.'])
