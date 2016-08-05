function t = ishrt(tform)
% ISHRT Check whether matrix is homogeneous rotation-translation matrix.
%   H = ISHRT(TFORM) returns true if matrix TFORM is an affine homogeneous 
%   transformation matrix that specifies rotation and translation only,
%   no scaling, sheering, etc.
%
%   Example:
%      ishrt(eye(4))

% Copyright 2016 Alexander Schaefer

t = false;

% Check whether all values are finite.
if ~all(isfinite(tform(:)))
    return
end

% Check correct matrix size.
if ~ismatrix(tform)
    return
end
if ~all(size(tform) == [4, 4])
    return
end

% Check correctness of rotation matrix.
rot = tform(1:3,1:3);
if det(rot) ~= 1
    return
end
if ~all(rot' == inv(rot))
    return
end

% Check last row.
if ~all(tform(4,:) == [0, 0, 0, 1])
    return
end

t = true;

end
