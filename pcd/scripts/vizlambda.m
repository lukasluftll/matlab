%% Read data.
% Load the point cloud data.
data = pcdread('data/castle.pcd');

% Create pointCloud object.
cloud = pointCloud(cat(3, data.x, cat(3, data.y, data.z)));

% Define the voxel resolution used to compute the normal distributed ray 
% decay rates.
res = 5;

% Compute the axis-aligned volume spanned by the point cloud.
vol = [floor([min(data.x(:)), min(data.y(:)), min(data.z(:))]/res), ...
        ceil([max(data.x(:)), max(data.y(:)), max(data.z(:))]/res)] * res;

% Make sure the points in the maximum plane of the volume are part of the
% volume.
vol(4:6) = vol(4:6) + eps(vol(4:6));

%% Compute decay rates and normal distributions.
% Calculate the ray decay rate.
lambda = raydecay(data.azimuth, data.elevation, data.radius, ...
    res, vol);

% Perform normal distribution transformation.
[mu, sigma] = ndt(cloud, res, vol);

% Reduce dimensionality of mean and covariance matrices.
mu = reshape(mu, 3, [])';
sigma = reshape(sigma, 3, 3, []);

% Use only covariance matrices that are positive definite.
i = 1;
while i <= size(sigma, 3)
    if ~any(any(isnan(sigma(:,:,i))))
        if min(eig(sigma(:,:,i))) > 0
            i = i + 1;
            continue
        end
    end
    
    sigma(:,:,i) = [];
    mu(i,:) = [];
end

%% Visualize mixture of normal distributions.
% Iterate over all voxels and compute the sum of all normal distributions 
% for each voxel.
cumlambda = zeros(size(lambda));
for x = 1 : size(lambda, 1)
    for y = 1 : size(lambda, 2)
        for z = 1 : size(lambda, 3)
            % Compute the center of the voxel.
            c = vol(1:3) + ([x, y, z] - 0.5*ones(1, 3)) * res;
            cumlambda(x,y,z) = sum(mvnpdf(c, mu, sigma) .* lambda(x,y,z));
        end
    end
end

cumlambda(isnan(cumlambda)) = 0;
pdfmax = max(cumlambda(:));
for x = 1 : size(cumlambda, 1)
    for y = 1 : size(cumlambda, 2)
        for z = 1 : size(cumlambda, 3)
            limits = vol(1:3) + ([x, y, z] - ones(1, 3))*res;
            limits = [limits, limits + res];
            cuboid(limits, 'FaceAlpha', cumlambda(x,y,z)/pdfmax);
        end
    end
end

axis equal
