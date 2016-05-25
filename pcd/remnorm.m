% Calculates surface normals using the intensity values contained in a
% point cloud. The normals are vectors that point from the surface to the
% sensor position from which the brightest remission was received.

%% Read data.
% Find a series of PCD files in the data directory.
datadir = '~/ros/datasets/bonirob/hangar_pcd/';
files = dir([datadir, 'velodyne-02*.pcd']);

waitbarHandle = waitbar(...
    0, ['Reading ', int2str(numel(files)), ' files ...']);

% Read the files.
for (i = 1 : numel(files))
    cloud = pcdread([datadir, files(i).name]);
    
    % Unorganize the point cloud for easier handling.
    cloud.pointCloud = pointCloud(reshape(...
        cloud.pointCloud.Location, [cloud.pointCloud.Count, 3]));
    cloud.intensity = reshape(...
        cloud.intensity, [cloud.pointCloud.Count, 1]);
    clouds{i} = cloud;
    
    % Store the sensor viewpoint of the cloud.
    viewpoint(i, :) = cloud.viewpoint.T(4, 1:3); %#ok<*SAGROW>
    
    % Advance the progress bar.
    waitbar(i / numel(files), waitbarHandle);
end

close(waitbarHandle);

% Use the cloud in the middle of the dataset as reference cloud whose 
% normals are computed.
refcloud = clouds{round(numel(clouds) / 2)};

%% Compute intensity normals.
% Initialize the normals matrix.
normals = zeros(size(refcloud.pointCloud.Location));

% Initialize the vector which contains the number of comparison 
% intensities for each point.
compPoints = zeros(size(refcloud.intensity));

% Determine the indices of the points whose normals are computed.
selection = 1 : refcloud.pointCloud.Count;

% Iterate over all points of the reference cloud in parallel.
parfor (p = selection)
    refpoint = refcloud.pointCloud.Location(p, :); %#ok<*PFBNS>
    
    % If the point is not a number, skip it.
    if (isnan(refpoint))
        continue;
    end
    
    % For each other point cloud in the dataset, determine the intensity of
    % the point that is the nearest neighbor of the considered point in the 
    % reference cloud. If this neighbor point's intensity is greater than 
    % the maximum intensity found so far, store the beam from the point to 
    % the sensor viewpoint as normal vector.
    maxIntensity = 0;
    for (c = 1 : numel(clouds))
        compCloud = clouds{c};
        
        [nearestIndex, dist] = findNearestNeighbors(...
            compCloud.pointCloud, refpoint, 1);
        
        % Create shortcuts.
        compPoint = compCloud.pointCloud.Location(nearestIndex, :);
        compIntensity = compCloud.intensity(nearestIndex);
        compViewpoint = viewpoint(c, :);
        
        % If the compared point's intensity is greater than the intensity
        % stored last, calculate and save the normal.
        if (compIntensity > maxIntensity && dist <= 0.5)
            maxIntensity = compIntensity;
            
            % Compute the normal.
            normals(p, :) = compViewpoint - compPoint;
            
            % Normalize the normal.
            normals(p, :) = normals(p, :) / norm(normals(p, :));
            
            compPoints(p) = compPoints(p) + 1;
        end
    end
end

% Zero out the normals that were computed using less than a minimum number
% of intensity values.
mincomp = 6;
normals(compPoints < mincomp, :) = zeros(...
    size(normals(compPoints < mincomp, :)));

%% Create map.
map = pccolor(clouds{1});
for (i = 2 : numel(clouds))
    map = pcmerge(map, pccolor(clouds{i}), 0.5);
end

%% Plot intensity normals.
intensityFig = figure('Name', 'Normals calculated using intensity', ...
    'NumberTitle', 'Off');

% Plot the map.
intensityax = pcshow(map);

% Plot the normals.
hold on;
nx = refcloud.pointCloud.Location(selection, 1);
ny = refcloud.pointCloud.Location(selection, 2);
nz = refcloud.pointCloud.Location(selection, 3);
nu = normals(selection, 1);
nv = normals(selection, 2);
nw = normals(selection, 3);
quiver3(nx, ny, nz, nu, nv, nw, '.');

% Plot the sensor viewpoints.
sx = viewpoint(:, 1);
sy = viewpoint(:, 2);
sz = viewpoint(:, 3);
plot3(sx, sy, sz, 'Marker', 'x', 'MarkerSize', 10);
hold off;

%% Compute local neighborhood normals.
neighborFig = figure(...
    'Name', 'Normals calculated using local neighborhood', ...
    'NumberTitle', 'Off');

% Plot the reference point cloud.
neighborax = pcshow(map);

% Find the nearest neighbors of the reference cloud points in the map.
nearestIndices = zeros(size(refcloud.intensity));
for (p = selection)
    point = refcloud.pointCloud.Location(p, :);
    
    % If the point is not a number, skip it.
    if (isnan(point))
        continue;
    end
    
    [nearestIndex, dist] = findNearestNeighbors(map, point, 1);
    nearestIndices(p) = nearestIndex;
end

% Compute the normals of the map.
normals = pcnormals(map);

% Delete normals corresponding to NaN points in the reference cloud.
selection = selection(nearestIndices > 0);
nearestIndices = nearestIndices(nearestIndices > 0);

%% Plot local neighborhood normals.
% Plot the normals.
hold on;
nx = refcloud.pointCloud.Location(selection, 1);
ny = refcloud.pointCloud.Location(selection, 2);
nz = refcloud.pointCloud.Location(selection, 3);
nu = normals(nearestIndices, 1);
nv = normals(nearestIndices, 2);
nw = normals(nearestIndices, 3);
quiver3(nx, ny, nz, nu, nv, nw, '.');

% Plot the sensor viewpoints.
plot3(sx, sy, sz, 'Marker', 'x', 'MarkerSize', 10);
hold off;

% Link the cameras of the two figures together.
cameralink = linkprop([intensityax, neighborax], {'CameraUpVector', ...
    'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', cameralink);
