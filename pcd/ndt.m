function [mu, sigma] = ndt(cloud, center, radius)
% NDT Normal distributions transform of point cloud.
%   [MU, SIGMA] = NDT(CLOUD, CENTER, RADIUS) performs normal distributions 
%   transforms on all points of the point cloud CLOUD that are contained in 
%   the spheres defined by CENTER and RADIUS.
%
%   CLOUD is a pointCloud object.
%
%   CENTER is an Nx3 matrix. CENTER(n,:) contains the x, y, and z 
%   coordinates of the center of the n-th sphere.
%
%   RADIUS is a scalar that defines the sphere radius. A point lies
%   inside a sphere if its distance to the sphere center is at most RADIUS.
%
%   MU is an Nx3 matrix. MU(n,:) is the mean position of all cloud points
%   inside the n-th sphere. The mean of a sphere that contains no points is 
%   NaN.
%
%   SIGMA is a 3x3xN matrix. SIGMA(:,:,n) is the 3x3 covariance matrix of 
%   the cloud points inside the n-th sphere. The covariance of a sphere 
%   that contains no points is NaN.
%
%   Example:
%      pc = pcread('teapot.ply');
%      [mu, sigma] = ndt(pc, [1 1 1; 2 1 1], 3)
%
%   See also NDTPLOT, POINTCLOUD, NAN.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(3, 3);

% Check the dimensions of the sphere center matrix.
if size(center, 2) ~= 3
    error('CENTER must be a Nx3 matrix.')
end

% Check if the radius is positive.
if radius <= 0
    error('RADIUS must be positive.')
end

%% Compute mean and covariance values.
% Construct the return matrices.
mu = [];
sigma = [];

% Make sure the point coordinates are contained in a 2D matrix.
location = reshape(cloud.Location, [], 3);

% Define the axis-aligned region-of-interest boxes that contain the 
% individual spheres.
roi = kron(center, [1, 1]) + repmat([-radius, +radius], size(center));

% Perform NDT for every sphere.
for i = 1 : size(center, 1)
    % Compute the indices of the points inside the ROI boxes.
    index = findPointsInROI(cloud, roi(i,:));

    % Get the Cartesian coordinates of the points inside the ROI box.
    roicloud = location(index,:);
    
    % Compute which points inside the ROI box are also located inside the 
    % sphere.
    dsq = sum((roicloud - repmat(center(i,:), size(roicloud,1), 1)).^2, 2);
    spherecloud = roicloud(dsq<=radius^2,:);
    
    % If there are less than 2 points inside the sphere, computing the 
    % covariance is impossible.
    if size(spherecloud, 1) < 2
        continue
    end
    
    % Compute the covariance of all points inside the sphere.
    sigmaTmp = cov(spherecloud);
    
    % Add the mean and covariance to the respective matrices only if the
    % covariance is finite and positive definite.
    if all(isfinite(sigmaTmp(:)))
        if all(eig(sigmaTmp()) >= 1e-12)
            mu = [mu; mean(spherecloud)]; %#ok<*AGROW>
            sigma = cat(3, sigma, sigmaTmp);
        end
    end
end

end
