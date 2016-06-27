function p = pdfray(origin, normray, lambda, xgv, ygv, zgv)
%

%% Validate input.

%%
% Compute the indices of the grid cells that the ray traverses.
[idx, t] = trav(origin, ray, xgv, ygv, zgv);

idx(t > 1,:) = [];
t(t > 1) = [];
t = t / norm(ray);

% Compute the lengths of the rays apportioned to each voxel.
l = diff(t);

N = ones(size(t));
i = 2;
while t(i) < norm(ray)
    N(i) = N(i-1) * exp(-lambda(i-1)*l(i-1));
end

p = lambda(i) * N(i) * exp(-lambda(i) * t(i));

end
