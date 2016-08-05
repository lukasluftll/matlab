function h = ishtform(tform)
% ISHTFORM Check whether matrix is homogeneous transformation matrix.
%   H = ISHTFORM(TFORM) returns true if matrix TFORM is a homogeneous
%   transformation matrix.
%
%   A homogeneous transformation matrix has the following properties:
%      - Its size is 4x4.
%      - It contains finite values only.
%      - The upper left 3x3 submatrix is a rotation matrix.
%      - The last row is [0, 0, 0, 1].
%
%   Example:
%      h = ishtform(eye(4))

% Copyright 2016 Alexander Schaefer

% Check finite values.
finiteOk = all(isfinite(tform(:)));

% Check size.
sizeOk = size(tform) == [4, 4];

% Check correctness of rotation matrix.
rotOk = tform(1:3,1:3)' == inv(tform(1:3,1:3));

% Check last row.
rowOk = tform(4,:) == [0, 0, 0, 1];

% Determine whether or not matrix is homogeneous transformation matrix.
h = finiteOk && sizeOk && rotOk && rowOk;

end

