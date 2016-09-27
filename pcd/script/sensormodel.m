pcfield = pcd2pc(pcdread('pcd/data/leek.pcd'));
em = elevationmap(pcfield, 0.05);

pcsens = pcd2pc(pcdread('pcd/data/sens.pcd'));
pcsens = pctransform(pcsens, ht2affine3d(eul2tform([pi,0,0])));

res = 0.1;
x = -100 : res : 20;
y = -20 : res : 10;
nx = numel(x);
ny = numel(y);
d = NaN(numel(x), numel(y));
z = mean(pcsens.Location(:,:,3), 'omitnan');

progressbar(nx)
parfor ix = 1 : nx
    for iy = 1 : ny
        zOpt = fminbnd(@(z) em.diff(pctransform(pcsens,ht2affine3d(...
            trvec2tform([x(ix),y(iy),z])))), -3, +3);
        pc = pctransform(pcsens,ht2affine3d(trvec2tform([x(ix),y(iy),zOpt])));
        d(ix,iy) = em.diff(pc); %#ok<*PFBNS>
    end
    progressbar
end

surf(d, 'EdgeColor', 'none')
labelaxes
