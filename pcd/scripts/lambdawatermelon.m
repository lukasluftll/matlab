% Computes the mean Lidar ray decay rate per voxel over all scans in the 
% watermelon dataset and visualizes them.

% Read all files in the folder that contains the watermelon dataset.
datasetPath = 'pcd/data/watermelon_pcd';
file = dir(datasetPath);

% Create the progress bar. 
waitbarHandle = waitbar(0, 'Parsing dataset ...');

% Remove all files that are no PCD files.
remove = [];
for i = 1 : numel(file)
    [~, ~, extension] = fileparts([datasetPath, '/', file(i).name]);
    if ~strcmpi(extension, '.pcd')
        remove(end+1) = i; %#ok<SAGROW>
    end
    
    % Advance the progress bar.
    waitbar(i / numel(file), waitbarHandle);
end
file(remove) = [];

% Create new progress bar.
close(waitbarHandle);
waitbarHandle = waitbar(0, 'Computing decay rates ...');

% Define the grid vectors.
res = 1;
pcd = pcdread([datasetPath, '/', file(1).name]);
xgv = min(pcd.x(:)) : res : max(pcd.x(:));
ygv = min(pcd.y(:)) : res : max(pcd.y(:));
zgv = min(pcd.z(:)) : res : max(pcd.z(:));

% Create the matrices that contain the total number of remissions per voxel
% and the total length of all rays that traversed each voxel.
rTot = zeros(numel(xgv)-1, numel(ygv)-1, numel(zgv)-1);
lTot = zeros(size(rTot));

% For each voxel, determine the number of remissions and the total ray 
% length.
tstart = tic;
for i = 1 : numel(file)
    pcd = pcdread([datasetPath, '/', file(i).name]);
    [~, r, l] = raydecay(...
        pcd.azimuth, pcd.elevation, pcd.radius, xgv, ygv, zgv);
    rTot = rTot + r;
    lTot = lTot + l;

    % Advance the progress bar.
    waitbar(i / numel(file), waitbarHandle);
    
    % Update the execution time estimate.
    set(get(findobj(waitbarHandle, 'type', 'axes'), 'title'), 'string', ...
        ['Computing decay rates; ', ...
        int2str((toc(tstart)/i) * (numel(file)-i)/60), ' min remaining']);
end
close(waitbarHandle);

% Compute the mean decay rate per voxel.
lambda = rTot ./ lTot;
lambda(lTot == 0) = NaN;

% Fit the decay rate into [0; 1].
lambda = lambda / max(lambda(:));

% Create a voxel plot that visualizes the mean decay rate per voxel.
alphaplot(lambda, xgv, ygv, zgv);
