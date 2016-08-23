% Creates a map from the laser scans collected while driving 
% a robot around on campus.

% Define parameters.
res = 0.3;
markerSize = 70;

% Create the progress bar.
waitbarHandle = waitbar(0, 'Registering point clouds ...');

% Find all PCD files in the data directory.
datadir = '~/ros/datasets/bonirob/campus_pcd/';
files = dir([datadir, 'velodyne_freiburg*.pcd']);

% Read the first point cloud.
campus = pcdread([datadir, files(1).name]);
campus = pccolor(campus);
campus = pcdownsample(campus, 'gridAverage', res);

% Save the minimum and maximum coordinate values.
limits = [campus.XLimits, campus.YLimits, campus.ZLimits];

% Display the first point cloud and configure the figure.
fig = figure;
ax = pcshow(campus, 'MarkerSize', markerSize);
whitebg(fig, [0.2, 0.2, 0.2]);

% Incrementally register every of the following point clouds 
% to the point clouds that are already registered.
for (i = 2 : numel(files))
    % Read the point cloud data.
    cloud = pcdread([datadir, files(i).name]);
    
    % Colorize the point cloud.
    cloud = pccolor(cloud);
    
    % Downsample the moving point cloud.
    cloud = pcdownsample(cloud, 'gridAverage', res);
    
    % Align the moving point cloud to the fixed point cloud.
    [tform, cloud] = pcregrigid(cloud, campus, ...
        'Metric', 'pointToPlane', ...
        'MaxIterations', 100, ...
        'Tolerance', [0.001, 0.001]);
    
    % Merge the fixed and the moving point cloud.
    campus = pcmerge(campus, cloud, res);
    
    % Update the limits.
    limits = max(limits, [campus.XLimits, campus.YLimits, campus.ZLimits]);
    
    % Display the result of the merge.
    pcshow(campus, 'MarkerSize', markerSize, 'Parent', ax);
    axis(limits);
    drawnow;
    
    % Advance the progress bar.
    waitbar(i / numel(files), waitbarHandle);
end

% Close the progress bar.
close(waitbarHandle);
