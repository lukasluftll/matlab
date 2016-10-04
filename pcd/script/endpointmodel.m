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
pcfield = pcd2pc(pcdread('pcd/data/leek.pcd'));

% Apply rotation to map.
pcfield = pctransform(pcfield, ht2affine3d(eul2tform([pi,0,0])));

% Create the elevation map and fill the gaps in the map.
em = fillnan(elevationmap(pcfield, 0.05), [5,5]);

%% Read sensor measurements.
% Read PCD file to point cloud.
pcsens = pcd2pc(pcdread(sensorfile));

% Rotate the scan.
pcsens = pctransform(pcsens, ht2affine3d(eul2tform(rpysens)));

%% Shift scan.
% Define the offset vectors.
x = shiftlim(1,1) : shiftres : shiftlim(1,2);
y = shiftlim(2,1) : shiftres : shiftlim(2,2);

% Initialize the matrix that stores the mean z difference.
d = NaN(numel(x), numel(y));

% Initialize the matrix that stores the fraction of NaN differences.
nanfrac = NaN(numel(x), numel(y));

% Vary the position of the scan and compute the z difference for
% each offset.
progressbar(numel(x))
for ix = 1 : numel(x)
    for iy = 1 : numel(y)
        d(ix,iy) = mean(pccompare(pcfield, pcsens));
    end
    
    % Advance the progress bar.
    progressbar
end

%% Plot results.
% Plot the mean z distance for all offsets.
fig = figure('Name', 'Mean distance to nearest neighbor');
surf(x, y, d', 'EdgeColor', 'none')
labelaxes
