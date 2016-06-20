function ndtplot(mu, sigma, vol, res)
% NDTPLOT Plot normal distributions transform of point cloud.

%% Validate input.
narginchk(2, 4)

if size(mu, 2) ~= 3
    error('')
end

if size(sigma, 1) ~= 3 || size(sigma, 2) ~= 3
    error('')
end

if size(mu, 1) ~= size(sigma, 3)
    error('')
end

if nargin < 3
    vol = [min(mu) max(mu)];
end

if any(vol(1:3) <= vol(4:6))
    error('')
end

if nargin < 4
    res = max(diff(reshape(vol, 3, 2), 1, 2)) / 100;
end

%%
% Discard all NaN values.
mu = mu(all(isfinite(mu), 2),:);
sigma = sigma(:,:,all(all(isfinite(sigma), 2)));

% Discard all normal distributions whose covariance matrices are not
% positive definite.
i = 1;
while i <= size(sigma, 3)
    if any(eig(sigma(:,:,i)) < 1e-12)
        sigma(:,:,i) = [];
        mu(i,:) = [];
    else
        i = i + 1;
    end
end

xgv = vol(1) : res : vol(4);
ygv = vol(2) : res : vol(5);
zgv = vol(3) : res : vol(6);

% Compute a Nx3 matrix that contains ...
[x, y, z] = ndgrid(xgv, ygv, zgv);
c = [x(:), y(:), z(:)];

% Evaluate the sum of all Gaussians at .
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
