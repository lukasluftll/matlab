function hrtchk(tform)
% HRTCK Validate 4x4 homogenous transformation matrix.
%   HRTCHK(TFORM) throws an error if TFORM is not a 4x4 homogeneous 
%   transformation matrix.
%
%   Example:
%      hrtchk(ones(4))
%
%   See also ISHRT.

% Copyright 2016 Alexander Schaefer

% Check the validity of the matrix size.
if ~ismatrix(tform) || any(size(tform) ~= 4)
    error('TFORM must be of size 4x4.')
end

% Check the validity of the matrix values.
if ~ishrt(tform)
    error('TFORM must be a 4x4 homogeneous coordinate transform.')
end

end
