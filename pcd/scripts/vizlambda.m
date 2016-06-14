% Load the point cloud data.
cloud = pcdread('data/castle.pcd');

% Define the voxel resolution used to compute the normal distributed ray 
% decay rates.
calcres = 5;

% Calculate the ray decay rate.
lambda = raydecay(cloud.azimuth, cloud.elevation, cloud.radius, calcres);

% Perform normal distribution transformation.
[mu, sigma] = ndt(cloud, calcres);

% Define the voxel resolution used to visualize the decay rate
% distribution.
vizres = 1;

% Visualize the sum of the normal distributions of the decay rates of 
% all voxels.
mvnpdf(mu, sigma)
calcres = 5;

ray = pc2sph(cloud.Location, 'vlp16');
lambda = raydecay(cloud, calcres);
[mu, sigma] = pc2gd(cloud, calcres);
