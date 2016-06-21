
res = 0.2;
radius = 1;
cloud = pcread('teapot.ply');
xgv = cloud.XLimits(1)-0.5 : res : cloud.XLimits(2)+0.5;
ygv = cloud.YLimits(1)-0.5 : res : cloud.YLimits(2)+0.5;
zgv = cloud.ZLimits(1)-0.5 : res : cloud.ZLimits(2)+0.5;
[x, y, z] = meshgrid(xgv, ygv, zgv);
center = [x(:), y(:), z(:)];
[mu, sigma] = ndt(cloud, center, radius);
density = reshape(ndpdf(mu, sigma, center), size(x));
isosurface(xgv, ygv, zgv, density, 170);
axis equal
