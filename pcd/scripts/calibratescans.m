% Estimates the pose of the second laser scanner w.r.t. the first one.
%
% The calibration works as follows:
% * Local maps are built for the first laser scanner using odometry only.
% * These local maps or clusters are then registered using 
%   the ICP algorithm, building a global map.
% * The same steps are repeated for the second laser.
% * The two resulting global maps are registered via ICP. In this way,
%   ICP provides the pose of the second laser scanner w.r.t. the first.

%% Preparation.
% Define parameters.
res = 0.100;
markerSize = 100;
clusterSize = 5;
palette = 'parula';
metric = 'pointToPoint';
datadir = '~/ros/datasets/bonirob/laser_calibration/';

% Create the progress bar.
waitbarHandle = waitbar(0);
set(findall(waitbarHandle, 'type', 'text'), 'Interpreter', 'none');

%% Build global maps from both scanners.
scanName = {'velodyne_freiburg', 'velodyne_bosch'};
% igm ... iterator over global maps
for (igm = 1 : numel(scanName))     
    % Update the progress bar text.
    waitbar(0, waitbarHandle, ['Building map for ', scanName{igm}, ' ...']);
    
    % Find all PCD files in the data directory.
    files = dir([datadir, scanName{igm}, '*.pcd']);
    
    % Display the first point cloud and configure the figure.
    fig = figure('Name', ['Map for ', scanName{igm}], 'NumberTitle', 'off');
    whitebg(fig, [0.2, 0.2, 0.2]);

    % Create a local map and register it to the global map.
    rmse = [];
    % ilm ... iterator over local maps
    for (ilm = 1 : clusterSize : numel(files)) 
        % Build the local map.
        for (ic = 1 : clusterSize)
            if (ic == 1)
                % Read the scan from file and colorize the point cloud.
                cluster = pcdread([datadir, files(ilm).name]);
                cluster = pccolor(cluster, palette);
            else
                % Read and colorize the scan and merge it with the 
                % other point clouds of the local map.
                cloud = pcdread([datadir, files(ilm+ic-1).name]);
                cloud = pccolor(cloud, palette);
                cluster = pcmerge(cluster, cloud, res);
            end
        end

        % Build the global map.
        if (ilm == 1)
            % Display the first local map.
            map{igm} = cluster;  %#ok<*SAGROW>
            ax = pcshow(map{igm}, 'MarkerSize', markerSize);
            axis(ax, [map{igm}.XLimits, map{igm}.YLimits, map{igm}.ZLimits]);
        else
            % Register the local map to the global map.
            [tform, cluster, e] = pcregrigid(cluster, map{igm}, ...
                'Metric', metric, ...
                'MaxIterations', 100, ...
                'Tolerance', [0.001, 0.001], ...
                'InlierRatio', 0.90, ...
                'Extrapolate', true, ...
                'Verbose', true); %#ok<ASGLU>
            rmse = [rmse; e]; %#ok<AGROW>

            % Merge the local and the global map.
            map{igm} = pcmerge(map{igm}, cluster, res);
        
            % Display the result of the merge.
            pcshow(map{igm}, 'MarkerSize', markerSize, 'Parent', ax);
            axis(ax, [map{igm}.XLimits, map{igm}.YLimits, map{igm}.ZLimits]);
        end
        
        % Update the figure and the progress bar.
        drawnow;
        waitbar(ilm / numel(files), waitbarHandle);
    end
    
    % Display the root mean squared error.
    figure('Name', ['RMSE for ', scanName{igm}], 'NumberTitle', 'off');
    histogram(rmse);
end

%% Register global maps.
% Update the progress bar.
waitbar(0, waitbarHandle, 'Registering maps ...');

% Register the maps and output the transform from the first scanner to
% the second.
[tform, map{end}, rmse] = pcregrigid(map{end}, map{1}, ...
    'Metric', metric, ...
    'MaxIterations', 1000, ...
    'Tolerance', [0.001, 0.001], ...
    'InlierRatio', 0.90, ...
    'Verbose', true);
transform = tform.T' %#ok<NOPTS>
RMSE = rmse %#ok<NOPTS>

% Save the transform to file.
source = scanName{1}; target = scanName{end};
save('calibration_result.mat', 'tform', 'source', 'target')

% Display the registered maps.
fig = figure('Name', 'Registered maps', 'NumberTitle', 'off');
ax = pcshowpair(pointCloud(map{1}.Location), ...
    pointCloud(map{end}.Location), ...
    'MarkerSize', markerSize);
whitebg(fig, [0.2, 0.2, 0.2]);

% Close the progress bar.
waitbar(1, waitbarHandle);
close(waitbarHandle);
