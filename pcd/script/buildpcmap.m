% Build a map out of many lidar scans.

%% Prepare output file.
% Print caption.
hline(75)
display(['Merging ', dataset, ' point cloud ...'])

% Create folder for results.
[resultFolderPath, resultFolderName] = fileparts(resultFolder);
if ~exist(resultFolder, 'dir')
    mkdir(resultFolderPath, resultFolderName);
end

% Define the name of the output MAT file.
pcMapFile = [resultFolder, '/pcmap_', dataset, '.mat'];

%% Merge point clouds.
% Get the PCD file names.
pcdFile = dir([folder, '/*.pcd']);
    
% Build the map and denoise it.
pcMap = pcdenoise(pcmap(folder, pcRes));

% Save the point cloud map to file.
save(pcMapFile, 'dataset', 'folder', 'pcRes', 'pcMap', '-v7.3')
display(['Result written to ', pcMapFile, '.'])
