function map = pcmap(folder, step)
% PCMAP Build map from multiple PCD files.
%   MAP = PCMAP(FOLDER, STEP) reads all PCD files in directory FOLDER and 
%   concatenates them to a single pointCloud object MAP. 
%   If the PCD files come with DAT files that provide odometry information, 
%   the MAP point cloud is built w.r.t. the odometry frame.
%
%   STEP specifies the grid step parameter used for subsampling. During 
%   subsampling, the volume occupied by both point clouds is voxelized. 
%   STEP is the edge length of the voxels. Whenever two points reside in 
%   one voxel, they are merged to a single point.

% Copyright 2016 Alexander Schaefer

%% Valiate input.
% Check the number of input arguments.
narginchk(2, 2)

% Check the FOLDER input argument.
if ~ischar(folder)
    error('FOLDER must be a string.')
end
if ~isdir(folder)
    error('FOLDER does not refer to a folder.')
end

%% Merge point clouds.
% Get the PCD file names.
file = dir([folder, '/*.pcd']);

% Create the map point cloud.
map = pointCloud(zeros(0, 3));

% Create a progress bar and set up automatic destruction after use.
waitbarHandle = waitbar(0, 'Building map ...');
cleanupObj = onCleanup(@() close(waitbarHandle));

% Loop over all PCD files.
for i = 1 : 10 : numel(file)
    % Read the PCD file.
    pcd = pcdread([folder, '/', file(i).name]);
    pc = pointCloud([pcd.x(:), pcd.y(:), pcd.z(:)]);
    
    % If odometry information is given, transform the point cloud
    % accordingly.
    if isfield(pcd, 'odometry')
        pc = pctransform(pc, ht2affine3d(pcd.odometry));
    end
    hold on; plotht(pcd.odometry); axis equal; pcd.odometry
    % Merge the point cloud with the map.
    map = pcmerge(map, pc, step);
    
    % Advance the progress bar.
    waitbar(i/numel(file), waitbarHandle);
end

end
