function ht = trquat2tform(tr, quat)
% TRQUAT2TFORM Translation and quaternion to homogeneous transformation.
%   HT = TRQUAT2TFORM(TR, QUAT) returns a 4x4xN matrix whose pages contain
%   N homogeneous transformation matrices. 
%
%   TR is an Nx3 matrix whose n-th row defines the translational part of 
%   the n-th homogeneous transformation.
%
%   QUAT is Nx4 matrix whose n-th row contains a unit quaternion that 
%   defines the rotational part of the n-th homogeneous transformation.
%
%   Example:
%      trquat2tform([1 2 3], [1 0 0 0])
%
%   See also TRVEC2TFORM, QUAT2TFORM.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(2, 2)

% Check size of input arguments.
if ~ismatrix(tr) || ~ismatrix(quat)
    error('TR and QUAT must be 2D matrices.')
end
if size(tr, 2) ~= 3
    error('TR must be of size Nx3.')
end
if size(quat, 2) ~= 4
    error('QUAT must be of size Nx4.')
end

% Check whether TR and QUAT contain the same number of translations and
% rotations.
if size(tr, 1) ~= size(quat, 1)
    error('TR and QUAT must contain the same number of rows.')
end

%% Build transformation matrix.
ht = quat2tform(quat);
ht(1:3,4,:) = reshape(tr.', 3, 1, []);

end
