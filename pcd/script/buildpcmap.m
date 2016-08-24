% Build a map out of many lidar scans.

%% Parameters.
% Dataset name.
dataset = 'campus';

% Dataset folder with PCD files.
folder = ['pcd/data/', dataset, '/pcd_sph'];

% Resolution of the merged point cloud map.
pcMapRes = 0.1;

%% Prepare output file.
% Create folder for results.
resultFolder = 'pcd/result';
[resultFolderPath, resultFolderName] = fileparts(resultFolder);
if ~exist(resultFolder, 'dir')
    mkdir(resultFolderPath, resultFolderName);
end

% Define the name of the output MAT file.
pcMapFile = [resultFolder, '/pcmap_', dataset, '.mat'];

% Save parameters to file.
save(pcMapFile, 'dataset', 'folder', 'pcMapRes', '-v7.3');

%% Merge point clouds.
% Get the PCD file names.
pcdFile = dir([folder, '/*.pcd']);

% Build the map.
pcMap = pcmap(folder, pcMapRes);
    
% Save the point cloud map to file.
save(pcMapFile, 'pcMap', '-append')
