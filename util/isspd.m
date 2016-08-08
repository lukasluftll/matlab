function y = isspd(x)
% ISSPD Checks if matrix is symmetric positive definite.
%   Y = ISSPD(X) checks whether the given matrix X is symmetric positive
%   definite. If so, it is returned as Y. If not, Y is empty.
%
%   If X is close to singular, ISSPD ensures the smaller eigenvalues of X 
%   are at least 0.001 of the largest eigenvalue. In this way, ISSPD 
%   guarantees the numerical stability of the matrix.
%
%   Example:
%      isspd(zeros(3))
%
%   See also EIG.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(1, 1)

%% Create return value.
y = [];

%% Check if matrix is symmetric.
if size(x, 1) ~= size(x, 2)
    return
end

%% Check if matrix values are finite.
if any(~isfinite(x(:)))
    return
end

%% Check for positive definiteness.
% The matrix is positive definite if all eigenvectors are positive.
[eigenvector, eigenvalue] = eig(x);
if any(diag(eigenvalue) <= 0)
    return
end

%% Ensure numerical stability.
% Compute the minimum stable eigenvalue.
mineig = max(diag(eigenvalue)) * 1e-3;

% Stabilize unstable eigenvalues and reconstruct the matrix.
if any(diag(eigenvalue) < mineig)
    eigenvalue([1,5,9]) = max([eigenvalue([1,5,9]); repmat(mineig, 1, 3)]);
    y = eigenvector * eigenvalue / eigenvector;
    return
end

% If all checks passed, return the original matrix.
y = x;

end
