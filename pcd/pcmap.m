function map = pcmap(folder, file, res, mode)
% PCMAP Build map from multiple PCD files.
%   MAP = PCMAP(FOLDER, RES) reads all PCD files in the folder FOLDER and 
%   merges them into a single pointCloud object MAP.
%
%   RES specifies the grid resolution parameter used for subsampling. 
%   During subsampling, the volume occupied by both point clouds is 
%   voxelized. RES is the edge length of the voxels. Whenever two points 
%   reside in one voxel, they are merged to a single point.
%
%   MAP = PCMAP(FOLDER, RES, MODE) defines the mode used when merging the
%   point cloud. The MODE string can take different values:
%
%   'direct' (default) - The point coordinates in the PCD files are 
%                        specified in the world coordinate system, so they 
%                        are merged without applying any coordinate 
%                        transformation.
%
%   'odometry'         - The point coordinates in the PCD files are 
%                        specified in the robot coordinate system, so they 
%                        are transformed into the world coordinate system 
%                        by applying the world-to-robot transformation 
%                        provided by the odometry property of the 
%                        accompanying DAT files.
%
%   PCMAP shows its progress in the command window.
%
%   Example:
%      pc = pcmap('pcd/data', 0.1, 'direct')
%
%   See also POINTCLOUD, PCREAD, PCDREAD.

% Copyright 2016 Alexander Schaefer

%% Valiate input.
% Check the number of input arguments.
narginchk(3, 4)

% Check the FILE input argument.
validateattributes(folder, {'char'}, {'nonempty'}, '', 'FOLDER')
if exist(folder, 'dir') ~= 7
    error('Folder %s does not exist.', folder)
end

if isempty(file)
    file = dir([folder, '/*.pcd']);
    if isempty(file)
        error('Folder %s does not contain any PCD files.', file)
    end
else
    for i = 1 : numel(file)
        if exist([folder, '/', file(i).name], 'file') ~= 2
            error('File %s does not exist.', file(i).name);
        end
    end
end

% Check the RES input argument.
validateattributes(res, {'numeric'}, {'numel',1, 'nonnegative'}, '', 'RES')

% Define the MODE input argument, if not specified.
if nargin < 4
    mode = 'direct';
end

validateattributes(mode, {'char'}, {'nonempty'}, '', 'MODE')

% Check the MODE input argument.
switch lower(mode)
    case 'direct'
    case 'odometry'
    otherwise
        error('MODE must be ''direct'' or ''odometry''.')
end    

%% Merge point clouds.
% Use multiple workers to create local maps.
disp('Merging local maps ...')
parprogress(numel(file));
spmd
    % Create the map point cloud for this worker.
    mapw = pointCloud(zeros(0, 3));
    
    % Compute the number of files per worker.
    fpw = max([1, round(numel(file)/numlabs)]);
    
    % Merge all of this worker's point clouds.
    for i = 1 + (labindex-1)*fpw : min([numel(file), labindex*fpw])
        % Read the PCD file.
        pcd = pcdread([folder, '/', file(i).name]);

        % Build a point cloud from the data in the PCD file depending on
        % whether the points are specified in Cartesian or spherical
        % coordinates.
        if all(isfield(pcd, {'x', 'y', 'z'})) % Cartesian PCD.
            pc = pointCloud([pcd.x(:), pcd.y(:), pcd.z(:)]);
        elseif all(isfield(pcd, {'sensor_x', 'sensor_y', 'sensor_z', ...
                'sensor_qw', 'sensor_qx', 'sensor_qy', 'sensor_qz', ...
                'azimuth', 'elevation', 'radius'})) % Spherical PCD.
            sp = trquat2tform([pcd.sensor_x(:), pcd.sensor_y(:), ...
                pcd.sensor_z(:)], ...
                [pcd.sensor_qw(:), pcd.sensor_qx(:), pcd.sensor_qy(:), ...
                pcd.sensor_qz(:)]);
            pc = ls2pc(laserscan(sp, pcd.azimuth(:), pcd.elevation(:), ...
                pcd.radius(:)));
        else
            error(['File ', file(i).name, ' has invalid format.'])
        end

        % Transform the point cloud, if necessary.
        switch lower(mode)
            case 'direct'
            case 'odometry'
                if ~isfield(pcd, 'odometry')
                    error(['No odometry available for file ', ...
                        file(i).name, '.'])
                end
                pc = pctransform(pc, ht2affine3d(pcd.odometry));    
        end

        % Merge the point cloud with the local map of the worker.
        mapw = pcmerge(mapw, pc, res);
        parprogress;
    end
end
parprogress(0);

% Merge the local maps of the workers to form a global map.
disp('Merging global map ...')
map = mapw{1};
parprogress(numel(mapw)-1);
for i = 2 : numel(mapw)
    map = pcmerge(map, mapw{i}, res);
    parprogress;
end
parprogress(0);

end
