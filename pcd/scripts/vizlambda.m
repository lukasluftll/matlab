% Load the point cloud data.
data = pcdread('data/castle.pcd');

% Create pointCloud object.
cloud = pointCloud(cat(3, data.x, cat(3, data.y, data.z)));

% Define the voxel resolution used to compute the normal distributed ray 
% decay rates.
calcres = 1;

% Compute the axis-aligned volume spanned by the point cloud.
vol = [min(data.x(:)), min(data.y(:)), min(data.z(:)), ...
        max(data.x(:)), max(data.y(:)), max(data.z(:))];

% Make sure the points in the maximum plane of the volume are part of the
% volume.
vol(4:6) = vol(4:6) + eps(vol(4:6));
vol = [0, 0, 0, 2, 1, 1];

% Calculate the ray decay rate.
lambda = raydecay(data.azimuth, data.elevation, data.radius, ...
    calcres, vol);

% Perform normal distribution transformation.
[mu, sigma] = ndt(cloud, calcres, vol);
muLin = reshape(mu(:), [], 3);
sigmaLin = reshape(sigma(:), [], 3, 3);
sigmaLin = permute(sigmaLin, [2, 3, 1]);
            
% Define the voxel resolution used to visualize the decay rate
% distribution.
vizres = 1;

% Iterate over all voxels and compute sum of all normal distributions for 
% each voxel.
for x = 1 : size(lambda, 1)
    for y = 1 : size(lambda, 2)
        for z = 1 : size(lambda, 3)
            % Compute the center of the voxel.
            c = vol(1:3) + ([x, y, z] - 0.5*ones(1, 3)) * calcres;
            p(x,y,z) = mvnpdf(c, muLin, sigmaLin) .* lambda(:);
        end
    end
end
