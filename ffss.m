function pred = ffss(pf, prev, u, sigma)

%% Validate input.
% Check number of input arguments.
narginchk(3, 3)

% Check size of previous particles matrix.
sp = size(pf.Particles);
if ~all(size(prev) == sp)
    error(['PREV must be a matrix of size ', num2str(sp(1)), ...
        'x', num2str(sp(2)), '.'])
end

% Check size of covariance matrix.
if ~all(size(sigma) == [sp(2), sp(2)])
    error(['SIGMA must be a matrix of size ', num2str(sp(2)), ...
        'x', num2str(sp(2)), '.'])
end

% Check if covariance matrix is symmetric positive-definite.
if ~isspd(sigma)
    error('SIGMA must be a symmetric positive-definite matrix.')
end

%% Sample from multivariate Gaussian.
% Sample new particles from multivariate Gaussian.
pred = mvnrnd(prev, sigma);

% Make sure RPY angles stay in interval [-pi; pi].
pred(:,4:6) = wrapToPi(pred(:,4:6));

end
