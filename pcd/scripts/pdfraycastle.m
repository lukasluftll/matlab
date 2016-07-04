pcd = pcdread('data/castle.pcd');

xgv = min(pcd.x(:)) : res : max(pcd.x(:));
ygv = min(pcd.y(:)) : res : max(pcd.y(:));
zgv = min(pcd.z(:)) : res : max(pcd.z(:));

lambda = raydecay(pcd.azimuth, pcd.elevation, pcd.radius, xgv, ygv, zgv);
p = prod(pdfray(origin, ray, lambda, xgv, ygv, zgv));
