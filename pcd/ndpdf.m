function y = ndpdf(mu, sigma, x)
% NDPDF Probability density function of multivariate normal distributions.
%   Y = NDPDF(MU, SIGMA, X) takes one or multiple multivariate normal 
%   distributions characterized by their mean MU and their covariance SIGMA
%   and for each point in X computes the sum of their probability
%   densities.
%
%   MU is a MxD matrix. The rows contain the means of the D-dimensional
%   normal distributions.
%
%   SIGMA is a DxDxM matrix. Its pages contain the DxD covariance matrices
%   of the normal distributions.
%
%   X is a NxD matrix. Its rows contain the N points where the sum of all
%   normal distributions is evaluated.
%
%   Y is a Nx1 matrix. Its elements contain the sum of the probability
%   density of all normal distributions evaluated at the point defined by
%   the corresponding row of X.
%
%   Example:
%      mu = [1, 1, 1; 5, 0, 0];
%      sigma = cat(3, eye(3), [2, -1, 0; -1, 2, -1; 0, -1, 2]);
%      x = [0, 0, 0; 4, 0, 1];
%      y = ndpdf(mu, sigma, x)
%
%   See also NDT, PDF.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(3, 3);

% Check if input argument dimensionality.
dim = size(mu, 2);
if size(sigma, 1) ~= dim || size(sigma, 2) ~= dim || size(x, 2) ~= dim
    error('MU, SIGMA, and X must all be D-variate.')
end

% Check if mu and sigma contain the same number of distributions.
if size(mu, 1) ~= size(sigma, 3)
    error('MU and SIGMA must specify the same number of distributions.')
end

% Remove all distributions whose covariance is not finite or not positive
% definite.
i = 1;
while i <= size(sigma, 3)
    if all(all(isfinite(sigma(:,:,i))))
        if all(eig(sigma(:,:,i))) >= 1e-12
            i = i + 1;
            continue
        end
    end
    mu(i,:) = [];
    sigma(:,:,i) = [];
end

%% Evaluate sum of normal distributions at given points.
y = zeros(size(x, 1), 1);
for i = 1 : size(mu, 1)
    y = y + mvnpdf(x, mu(i,:), sigma(:,:,i));
end

end
