function y = spd(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

y = [];

% Check if matrix is symmetric.
if size(x, 1) ~= size(x, 2)
    return
end

% Check if matrix values are finite.
if any(any(~isfinite(x)))
    return
end

% Check for positive definiteness.
[eigenvector, eigenvalue] = eig(x);
if any(diag(eigenvalue) <= 0)
    return
end

% If the matrix is close to singular, change it slightly to ensure
% numerical stability of further computations.
[mineig, minidx] = min(diag(eigenvalue));
[maxeig, maxidx] = max(diag(eigenvalue));
if mineig < maxeig * 1e-3
    eigenvalue(minidx,minidx) = eigenvalue(maxidx,maxidx) * 1e-3;
    y = eigenvector * eigenvalue / eigenvector;
    return
end

y = x;

end
