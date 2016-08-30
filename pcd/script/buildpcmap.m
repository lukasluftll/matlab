% Build a map out of many lidar scans.

%% Parameters.
% Dataset name.
dataset = 'demo';

% Dataset folder containing the PCD files.
folder = ['pcd/data/', dataset, '/pcd_sph'];

% Resolution of the merged point cloud map.
pcRes = 0.1;

%% Prepare output file.
% Print caption.
hline(75)
display(['Merging ', dataset, ' point cloud ...'])

% Create folder for results.
resultFolder = 'pcd/result';
[resultFolderPath, resultFolderName] = fileparts(resultFolder);
if ~exist(resultFolder, 'dir')
    mkdir(resultFolderPath, resultFolderName);
end

% Define the name of the output MAT file.
pcMapFile = [resultFolder, '/pcmap_', dataset, '.mat'];

%% Merge point clouds.
% Get the PCD file names.
pcdFile = dir([folder, '/*.pcd']);
    
% Build the map.
pcMap = pcmap(folder, pcRes);

% Save the point cloud map to file.
save(pcMapFile, 'dataset', 'folder', 'pcRes', 'pcMap', '-v7.3')
display(['Result written to ', pcMapFile, '.'])
