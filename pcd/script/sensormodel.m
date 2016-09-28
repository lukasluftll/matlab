% Evaluate sensor model for BoniRob localization on leek field

% Create elevation map of field.
pcfield = pcd2pc(pcdread('pcd/data/leek.pcd'));
pcfield = pctransform(pcfield, ht2affine3d(eul2tform([pi,0,0])));
em = elevationmap(pcfield, 0.05);

% Read sensor measurements.
pcsens = pcd2pc(pcdread('pcd/data/sensmiddle.pcd'));
psens = permute(pcsens.Location, [2,3,1]);

% Shift the measurements and evaluate the height difference at each point.
res = 0.1;
x = -10 : res : 100;
y = -10 : res : 20;
nx = numel(x);
ny = numel(y);
dsens = NaN(numel(x), numel(y));
nsens = pcsens.Count;
progressbar(nx)
parfor ix = 1 : nx
    for iy = 1 : ny       
        offset = repmat([x(ix),y(iy),0], nsens, 1) %#ok<PFBNS>
        offset(:,3) = em - (psens+offset);
        dsens(ix,iy) = -mean(abs(em-(psens+offset)), 'omitnan');
    end
    progressbar
end

figure('Name', 'Sensor model in middle of field')
surf(y, x, dsens, 'EdgeColor', 'none')
labelaxes
