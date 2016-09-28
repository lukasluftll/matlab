% Evaluate sensor model for BoniRob localization on leek field

% Create elevation map of field.
pcfield = pcd2pc(pcdread('pcd/data/leek.pcd'));
pcfield = pctransform(pcfield, ht2affine3d(eul2tform([pi,0,0])));
em = elevationmap(pcfield, 0.05);

% Read sensor measurements.
pcstart = pcd2pc(pcdread('pcd/data/sensstart.pcd'));
pcmiddle = pcd2pc(pcdread('pcd/data/sensmiddle.pcd'));
pstart = permute(pcstart.Location, [2,3,1]);
pmiddle = permute(pcmiddle.Location, [2,3,1]);

% Shift the measurements and evaluate the height difference at each point.
res = 0.1;
x = -10 : res : 100;
y = -10 : res : 20;
nx = numel(x);
ny = numel(y);
dstart = NaN(numel(x), numel(y));
dmiddle = dstart;
nmiddle = pcmiddle.Count;
progressbar(nx)
parfor ix = 1 : nx
    for iy = 1 : ny
        %offset = repmat([x(ix),y(iy),0], pcstart.Count, 1); %#ok<*PFBNS>
        %dstart(ix,iy) = -mean(abs(em-(pstart+offset)), 'omitnan');
       
        offset = repmat([x(ix),y(iy),0], nmiddle, 1) %#ok<PFBNS>
        offset(:,3) = em - (pmiddle+offset);
        dmiddle(ix,iy) = -mean(abs(em-(pmiddle+offset)), 'omitnan');
    end
    progressbar
end

% Plot the results.
%figure('Name', 'Sensor model at border of field')
%surf(y, x, dstart, 'EdgeColor', 'none')
%labelaxes

figure('Name', 'Sensor model in middle of field')
surf(y, x, dmiddle, 'EdgeColor', 'none')
labelaxes
