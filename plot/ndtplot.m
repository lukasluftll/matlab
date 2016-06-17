function ndtplot(cloud, ndtres, plotres)
% NDTPLOT Plot normal distributions transform of point cloud.
%   ... evaluated at the center of the voxel ...

%% Validate input.
narginchk(1,3)

%%
% Compute the axis-aligned volume spanned by the point cloud.
vol = reshape([cloud.XLimits; cloud.YLimits; cloud.ZLimits], 1, 6);
    
% Make sure the points in the maximum plane of the volume are part of the
% volume.
vol(4:6) = vol(4:6) + eps(vol(4:6));

% Perform the normal distributions transform.
[mu, sigma] = ndt(cloud, ndtres);

% Discard all NaN values.
mu = reshape(mu(isfinite(mu)), 3, [])';
sigma = reshape(sigma(isfinite(sigma)), 3, 3, []);

% Discard all voxels whose corresponding covariance matrices are not
% positive definite.
i = 1;
while i <= size(mu, 1)
    if any(eig(sigma(:,:,i)) < 1e-12)
        sigma(:,:,i) = [];
        mu(i,:) = [];
    else
        i = i + 1;
    end
end

% Compute the volume that spans from the beginning of the first voxel to 
% the end of the last.
plotvol = [floor(vol(1:3)/plotres), ceil(vol(4:6)/plotres)] * plotres;

% Compute the centers of the voxels.
cx = plotvol(1) : plotres : plotvol(4) + plotres/2;
cy = plotvol(2) : plotres : plotvol(5) + plotres/2;
cz = plotvol(3) : plotres : plotvol(6) + plotres/2;

% Compute a Nx3 matrix that contains all center points.
[y, x, z] = meshgrid(cy, cx, cz);
c = [x(:), y(:), z(:)];

% Evaluate the mixture of all Gaussians at the center points.
density = zeros(size(c, 1), 1);
for i = 1 : size(mu, 1)
    density = density + mvnpdf(c, mu(i,:), sigma(:,:,i));
end

density = reshape(density, size(x));

% Fit lambda into [0; 1].
density = density / max(density(:));

%% Visualize .
% Visualize the decay rate.
alphaplot(density, plotvol);
 
% Set the visualization parameters.
axis equal; xlabel('x'); ylabel('y'); zlabel('z'); grid on

end
