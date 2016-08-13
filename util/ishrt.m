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
if ~all(size(tform) == 4)
    return
end

% Check correctness of rotation matrix.
epsFactor = 100;
rot = tform2rotm(tform);
if abs(det(rot)-1) > eps(epsFactor)
    return
end
if max(abs(rot'-inv(rot))) > eps(epsFactor * max(rot(:)))
    return
end

% Check last row.
if ~all(tform(4,:) == [0, 0, 0, 1])
    return
end

t = true;

end
