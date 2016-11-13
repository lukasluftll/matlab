function m = quatmean(q, w)
% QUATMEAN Weighted average of set of quaternions.
%   M = QUATMEAN(Q) computes the average orientation resulting from the
%   given quaternions Q.
%   Q is a 4xN matrix, whose columns are quaternions.
%   M is a 4-element column vector containing the average quaternion.
%
%   M = QUATMEAN(Q,W) computes the weighted average. 
%   W is a N-element row vector that contains the weights for the
%   individual quaternions in Q. The n-th element represents the weight of
%   quaternion Q(:,n).
%   
%   Example:
%      q = eul2quat([0,0,0; 1,2,3; 3,2,0])'
%      quat2eul(quatmean(q, [0.3,0.3,0.4])')
%
%   See also ANGLEMEAN, MEAN.

% Copyright 2016 Alexander Schaefer
%
% QUATMEAN implements the quaternion averaging algorithm presented by 
% Markley et al.:
% F. Landis Markley, Yang Cheng, John Lucas Crassidis, and Yaakov Oshman. 
% Averaging Quaternions.
% Journal of Guidance, Control, and Dynamics, 30(4):1193-1197, 2007.  

%% Validate input and output.
% Validate numbers of output and input arguments.
nargoutchk(0, 1)
narginchk(1, 2)

% Validate the quaternion matrix.
validateattributes(q, {'numeric'}, {'real', 'rows', 4}, '', 'Q')
nq = size(q,2);

% Validate the weight vector.
if nargin < 2
    w(1:nq) = 1/nq;
end
validateattributes(w, {'numeric'}, ...
    {'real', 'positive', 'rows', 1, 'cols', nq}, '', 'W')

%% Average quaternions.
M = zeros(4);
for i = 1 : size(q,2)
    M = M + w(i) * q(:,i) * q(:,i)';
end

[m,~] = eigs(M,1);

end
