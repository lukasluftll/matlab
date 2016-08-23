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

% Preallocate resulting map.
pcMap = pointCloud(zeros(0, 3));

% Iterate over all PCD files.
progressbar('Merging point cloud map ...')
for i = 1 : numel(pcdFile)
    % Read laser scan data from file.
    ls = lsread([folder, '/', pcdFile(i).name]);
    
    % Merge the point cloud made from the laser scan with the map.
    pcMap = pcmerge(pcMap, ls2pc(ls), pcMapRes);
    
    % Advance the progress bar.
    progressbar(i/numel(pcdFile));
end

% Save the point cloud map to file.
save(pcMapFile, 'pcMap', '-append');
