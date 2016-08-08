function y = ndpdf(x, mu, Sigma, w)
% NDPDF Total probability density of multiple normal distributions.
%   Y = NDPDF(X, MU, SIGMA) takes one or multiple multivariate normal 
%   distributions characterized by their mean MU and covariance SIGMA and 
%   for each point in X computes the sum of their probability densities.
%
%   Y = NDPDF(X, MU, SIGMA, W) weights the distributions according to the 
%   factors specified by W.
%
%   MU is a MxD matrix. The rows contain the means of the M D-dimensional
%   normal distributions.
%
%   SIGMA is a DxDxM matrix. Its pages contain the DxD covariance matrices
%   of the M normal distributions.
%
%   W is a Mx1 vector. Its elements contain the weighting factors of the
%   corresponding normal distributions.
%
%   X is a NxD matrix. Its rows contain the N points where the sum of all
%   normal distributions is evaluated.
%
%   Y is a Nx1 vector. Its elements contain the sum of the probability
%   densities of all normal distributions evaluated at the corresponding 
%   rows of X.
%
%   Example:
%      mu = [1, 1, 1; 5, 0, 0];
%      sigma = cat(3, eye(3), [2, -1, 0; -1, 2, -1; 0, -1, 2]);
%      x = [0, 0, 0; 4, 0, 1];
%      y = ndpdf(x, mu, sigma)
%
%   See also NDT, PDF.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(3, 4);

% If no weights are given, set the weighting vector to 1.
if nargin < 4
    w = ones(size(mu, 1), 1);
end

% Check if input argument dimensionality.
dim = size(mu, 2);
if size(Sigma, 1) ~= dim || size(Sigma, 2) ~= dim || size(x, 2) ~= dim
    error('MU, SIGMA, and X must all be D-variate.')
end
if size(w, 2) ~= 1
    error('W must be a Mx1 vector.')
end

% Check if mu, sigma, and w contain the same number of distributions.
if size(mu, 1) ~= size(Sigma, 3) || size(mu, 1) ~= size(w, 1)
    error('MU, SIGMA, W must specify the same number of distributions.')
end

% Remove all distributions whose covariances are not finite or not positive
% definite.
remove = false(size(Sigma, 3), 1);
for i = 1 : size(Sigma, 3)
    % Check if the covariance matrix is symmetric positive definite.
    SigmaStable = isspd(Sigma(:,:,i));
    
    % If it is not positive definite, remove the corresponding normal 
    % distribution. If it is close to singular, change it slightly to 
    % ensure numerical stability. 
    if isempty(SigmaStable)
        remove(i) = true;
    else
        Sigma(:,:,i) = SigmaStable;
    end
end
Sigma(:,:,remove) = [];
mu(remove,:) = [];
w(remove) = [];

%% Evaluate sum of normal distributions at given points.
y = zeros(size(x, 1), 1);
for i = 1 : size(mu, 1)
    y = y + w(i) * mvnpdf(x, mu(i,:), Sigma(:,:,i));
end

end
