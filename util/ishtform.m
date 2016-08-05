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

%% Check matrix properties.
h = false;

% Check finite values.
if ~all(isfinite(tform(:)))
    return
end

% Check size.
if ~ismatrix(tform)
    return
end
if ~all(size(tform) == [4, 4])
    return
end

% Check correctness of rotation matrix.
if ~all(tform(1:3,1:3)' == inv(tform(1:3,1:3)))
    return
end

% Check last row.
if ~all(tform(4,:) == [0, 0, 0, 1])
    return
end

h = true;

end

