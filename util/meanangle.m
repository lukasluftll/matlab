function m = meanangle(theta)
% MEANANGLE Mean of set of angles in radians.
%   M = MEANANGLE(THETA) computes the mean of THETA. M minimizes the
%   squared arc length error.
%
%   THETA is a n-D matrix that contains a set of angles in radians.
%
%   M is the mean angle wrapped to [0; 2*pi].
%
%   Example:
%      meanangle([-0.5 : 0.1 : 0.2])
%
%   See also MEAN.

% Copyright 2016 Alexander Schaefer
%
% MEANANGLE implements the orientation averaging algorithm proposed by 
% Olson:
% Edwin Olson. On computing the average orientation of vectors and lines.
% 2011 IEEE International Conference on Robotics and Automation,
% Shanghai, China.

%% Validate input and output.
nargoutchk(0, 1)
narginchk(1, 1)
validateattributes(theta, {'numeric'}, {'real'}, '', 'THETA')

%% Compute mean.
% Notation is inspired by Olson's paper.

% Map all angles to [0; 2*pi] and sort them.
theta = sort(wrapTo2Pi(theta(:)))';

% Compute the moments for all angle arrangements.
N = numel(theta);
M1 = sum(theta) + (0:N-1) * 2*pi;
M2 = sum(theta.^2) + [0, cumsum(4*pi*(theta(1:end-1) + pi))];

% Compute the mean and the variance for all arrangements.
m = M1 / N;
sigma = M2 - 2*M1.*m + N*m.^2;

% Return the mean that minimizes the squared error.
[~,imin] = min(sigma);
m = wrapTo2Pi(m(imin));

end
