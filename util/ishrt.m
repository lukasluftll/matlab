function t = ishrt(tform)
% ISHRT Verify homogeneous rotation-translation matrices.
%   H = ISHRT(TFORM) checks whether each transformation matrix in TFORM is
%   a 4x4 affine homogeneous transformation matrix that specifies rotation 
%   and translation only.
%
%   TFORM is a 4x4xN matrix, whose pages contain N transformation matrices.
%
%   T is an N-element logical column vector. If TFORM(:,:,n) is a
%   homogeneous rotation-translation matrix, T(n) is true; otherwise, it is
%   false.
%
%   Example:
%      ishrt(eye(4))

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check matrix size.
if ndims(tform) > 3
    error('TFORM must be 2D or 3D.')
end
if size(tform, 1) ~= 4 || size(tform, 2) ~= 4
    error('TFORM must be of size 4x4xN.')
end

%% Perform check.
% Define tolerance factor for double-to-double equality comparisons.
kEps = 100;

% Check last row of all transformation matrices.
t = all(tform(4,:,:) == repmat([0,0,0,1], 1, 1, size(tform, 3)), 2);

% Check whether all values are finite.
t(t) = t(t) & all(all(isfinite(tform(:,:,t))));
    
% Check whether rotation matrices are valid.
parfor i = 1 : numel(t)
    rot = tform(1:3,1:3,i);
    t(i) = t(i) && abs(det(rot)-1) < eps(kEps) ...
        && max(max(abs(rot.' - inv(rot)))) < eps(kEps*max(rot(:)));
end

% Convert result to column vector.
t = t(:);

end
