% Evaluate endpoint sensor model for leek field dataset.

%% Define parameters.
% Step size when shifting the scan.
shiftres = 0.1;

% Minimum and maximum shifting offset in x and y direction.
shiftlim = [-10, 100; -10, 20];

% PCD file containing the map.
mapfile = 'pcd/data/leek.pcd';

% PCD file containing the scan.
sensorfile = 'pcd/data/sensmiddle.pcd';

% Orientation of the scan.
rpysens = deg2rad([0,0,0]);

%% Create elevation map of field.
% Read PCD file to point cloud.
pcfield = removeInvalidPoints(pcd2pc(pcdread('pcd/data/leek.pcd')));

% Apply rotation to map.
pcfield = pctransform(pcfield, ht2affine3d(eul2tform([pi,0,0])));

%% Read sensor measurements.
% Read PCD file to point cloud.
pcsens = removeInvalidPoints(pcd2pc(pcdread(sensorfile)));

% Rotate the scan.
pcsens = pctransform(pcsens, ht2affine3d(eul2tform(rpysens)));

%% Shift scan.
% Define the offset vectors.
x = shiftlim(1,1) : shiftres : shiftlim(1,2);
y = shiftlim(2,1) : shiftres : shiftlim(2,2);
nx = numel(x);
ny = numel(y);

% Initialize the matrix that stores the mean z difference.
d = NaN(nx, ny);

% Initialize the matrix that stores the fraction of NaN differences.
nanfrac = NaN(nx, ny);

% Vary the position of the scan and compute the z difference for
% each offset.
progressbar(nx)
pfield = pcfield.Location;
psens = pcsens.Location;
parfor ix = 1 : nx
    for iy = 1 : ny
        [~,dist] = knnsearch(pfield, psens);
        d(ix,iy) = mean(dist);
    end
    
    % Advance the progress bar.
    progressbar
end

%% Plot results.
% Plot the mean z distance for all offsets.
fig = figure('Name', 'Mean distance to nearest neighbor');
surf(x, y, d', 'EdgeColor', 'none')
labelaxes
