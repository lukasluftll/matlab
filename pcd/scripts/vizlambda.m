cloud = pcdread('data/castle.pcd');
cloud = cloud.pointCloud;

res = 10;

ray = pc2sph(cloud.Location, 'vlp16');
lambda = raydecay(cloud, res);
[mu, sigma] = pc2gd(cloud, res);
