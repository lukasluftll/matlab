% Evaluate eleavation map sensor model for leek field dataset.

% Define parameters.
shiftres = 0.1;
shiftlim = [-10, 100; -10, 20];
sensorfile = 'pcd/data/sensmiddle.pcd';
roifield = [-Inf,+Inf,-Inf,+Inf,-Inf,+Inf]; %[1.6 87.7 0 13.27 -Inf +Inf];
roisens = [-Inf,+Inf,-Inf,+Inf,-Inf,+Inf]; %[-29 30 -6.11 8 -Inf +Inf]
rpysens = [deg2rad(10),0,0];

% Create elevation map of field.
pcfield = pcd2pc(pcdread('pcd/data/leek.pcd'));
pcfield = pctransform(pcfield, ht2affine3d(eul2tform([pi,0,0])));
pcfield = select(pcfield, findPointsInROI(pcfield, roifield)); 
emnan = elevationmap(pcfield, 0.05);

% Compute mean elevation.
emean = mean(emnan.elevation(:), 'omitnan');

% Fill gaps in elevation map.
em = emnan.fillnan([5,5]);

% Read sensor measurements.
pcsens = pcd2pc(pcdread(sensorfile));
pcsens = select(pcsens, findPointsInROI(pcsens, roisens)); 
pcsens = pctransform(pcsens, ht2affine3d(eul2tform(rpy)));
psens = pcsens.Location;

% Shift the measurements and compute the height difference at each point.
x = shiftlim(1,1) : shiftres : shiftlim(1,2);
y = shiftlim(2,1) : shiftres : shiftlim(2,2);
nx = numel(x);
ny = numel(y);
d = NaN(nx, ny);
nanfrac = NaN(nx, ny);
n = size(psens, 1);
progressbar(nx)
parfor ix = 1 : nx
    for iy = 1 : ny
        % Compute the disparity between elevation map and scan.
        zoffset = -mean(em.diff(psens+repmat([x(ix),y(iy),0],n,1)), ...
            'omitnan');
        p = psens + repmat([x(ix),y(iy),zoffset], n, 1);
        dz = constrain(em.diff(p), [0,+Inf]);
        if sum(isnan(dz)) < 0.1*n
            d(ix,iy) = mean(dz, 'omitnan');
        end
        
        % Compute the fraction of unmatched points.
        nanfrac(ix,iy) = sum(isnan(dz)) / numel(dz);
    end
    progressbar
end

% Plot result.
figure('Name', 'Sensor model in middle of field')
surf(x, y, d', 'EdgeColor', 'none')
daspect([1 1 0.01])
labelaxes

% Colorize field map according to visualize the sensor model output.
% c = floor(normm(d) * 63.9999) + 1;
% palette = uint8(round(colormap * 255));
% roi = [lim(1,1), lim(1,2)-eps, lim(2,1), lim(2,2)-eps, -Inf, +Inf];
% pcfieldsel = select(pcfield, findPointsInROI(pcfield, roi));
% minxy = repmat(lim(1:2,1).', pcfieldsel.Count, 1);
% i = floor((pcfieldsel.Location(:,1:2) - minxy) / res) + 1;
% i = sub2ind(size(c), i(:,1), i(:,2));
% pcfieldsel.Color = palette(c(i),:);
% figure('Name', 'Field map showing robot location probability')
% pcshow(pcfieldsel, 'MarkerSize', 40)

% Show the fraction of missing correspondences.
% figure('Name', 'NaN count')
% surf(x, y, nanfrac', 'EdgeColor', 'none')
% labelaxes

% Plot point cloud of field and the scan point cloud at the most probable
% location.
[~,imin] = min(d(:));
[xmin,ymin] = ind2sub(size(d), imin);
offset = repmat([x(xmin),y(ymin),0], n, 1);
offset(:,3) = -mean(em.diff(psens+offset), 'omitnan');
figure('Name', 'Field')
pcshowpair(pcfield, pointCloud(psens+offset), 'MarkerSize', 50)
axis equal
labelaxes

% Vary yaw angle.
yaw = -pi : 0.01 : pi;
dyaw = NaN(numel(yaw), 1);
for iyaw = 1 : numel(yaw)
    pcsensrot = pctransform(pcsens, ...
        ht2affine3d(eul2tform([yaw(iyaw),0,0])));

    % Compute the disparity between elevation map and scan.
    p = pcsensrot.Location + offset;
    dz = constrain(em.diff(p), [0,+Inf]);
    dyaw(iyaw) = mean(dz, 'omitnan');
end
figure('Name', 'Matching result of best position over yaw angle')
plot(yaw, dyaw)

% Vary pitch angle.
pitch = -pi : 0.01 : pi;
dpitch = NaN(numel(pitch), 1);
for ipitch = 1 : numel(pitch)
    pcsensrot = pctransform(pcsens, ...
        ht2affine3d(eul2tform([0,pitch(ipitch),0])));
    
    % Compute the disparity between elevation map and scan.
    p = pcsensrot.Location + offset;
    dz = constrain(em.diff(p), [0;+Inf]);
    dpitch(ipitch) = mean(dz, 'omitnan');
end
figure('Name', 'Matching result of best position over pitch angle')
plot(pitch, dpitch)

% Vary roll angle.
roll = -pi : 0.01 : pi;
droll = NaN(numel(roll), 1);
for iroll = 1 : numel(roll)
    pcsensrot = pctransform(pcsens, ...
        ht2affine3d(eul2tform([0,0,roll(iroll)])));
    
    % Compute the disparity between elevation map and scan.
    p = pcsensrot.Location + offset;
    dz = constrain(em.diff(p), [0;+Inf]);
    droll(iroll) = mean(dz, 'omitnan');
end
figure('Name', 'Matching result of best position over roll angle')
plot(roll, droll)
