function tv = ht2tv(ht)
% HT2TV Homogeneous transformation matrix to translation vector.
%   TV = HT2TV(HT) transforms the homogeneous transformation matrices HT to
%   translation vectors TV.
%
%   HT is a 4x4xN matrix whose pages contain N transformation matrices.
%
%   TV is a Nx3 matrix, whose n-th row contains the translation vector of
%   the transformation matrix on the n-th page of HT.
%
%   In contrast to function TFORM2TRVEC, HT2TV can deal with an empty input
%   matrix.
%
%   Example:
%      ht2tv(eye(4))
%
%   See also TFORM2TRVEC, TRVEC2TFORM.

% Copyright 2016 Alexander Schaefer

%% Validate input.
narginchk(1, 1)

% If HT is empty, return.
tv = zeros(0, 3);
if isempty(ht)
    return
end

% Check contents of HT.
if ~isnumeric(ht)
    error('HT must contain numeric values only.')
end

% Check size of HT.
if ndims(ht) > 3 || size(ht, 1) ~= 4 || size(ht, 2) ~= 4
    error('HT must be a 4x4xN matrix.')
end

%% Get transformation vector.
tv = tform2trvec(ht);

end
