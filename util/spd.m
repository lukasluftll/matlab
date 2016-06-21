function y = spd(x)
% SPD Checks if matrix is symmetric positive definite.
%   Y = SPD(X) checks whether the given matrix X is symmetric positive
%   definite. If so, it is returned as Y. If not, Y is empty.
%
%   If X is close to singular, SPD changes the smallest eigenvalue of X so
%   it is 0.001 of the largest eigenvalue. In this way, SPD ensures the
%   numerical stability of the matrix.
%
%   Example:
%      spd(zeros(3))
%
%   See also EIG.

% Copyright 2016 Alexander Schaefer

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
% Check if the matrix is close to singular.
[mineig, minidx] = min(diag(eigenvalue));
[maxeig, maxidx] = max(diag(eigenvalue));
if mineig < maxeig * 1e-3
    % Change the minimum eigenvector and reconstruct the matrix.
    eigenvalue(minidx,minidx) = eigenvalue(maxidx,maxidx) * 1e-3;
    y = eigenvector * eigenvalue / eigenvector;
    return
end

% If all checks passed, return the original matrix.
y = x;

end
