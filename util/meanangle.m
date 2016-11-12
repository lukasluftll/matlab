function m = meanangle(theta, method)
% MEANANGLE Mean of set of angles in radians.
%   M = MEANANGLE(THETA) computes the mean of THETA. 
%
%   THETA is a n-D matrix that contains a set of angles in radians.
%
%   M is the mean angle wrapped to [-pi; +pi].
%
%   M = MEANANGLE(THETA, METHOD) specifies the algorithm used to compute
%   the mean:
%   'arcmin'    - Minimizes the squared arc lengths between the mean angle 
%                 and the given angles when drawn as points on the unit
%                 circle.
%   'vectorsum' - Compute the unit vector corresponding to each angle, sum
%                 up the vectors, and compute the arc tangent. This method
%                 minimizes the squared chord lengths between the mean 
%                 angle and the given angles when drawn as points on the 
%                 unit circle.
%   Default is 'arcmin'.
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
% Check the numbers of output and input arguments.
nargoutchk(0, 1)
narginchk(1, 2)

% If given no data, return no data.
if isempty(theta)
    m = [];
    return
end

% Check the given angles and ensure they are organized in a row vector.
validateattributes(theta, {'numeric'}, {'real'}, '', 'THETA')
theta = reshape(theta, 1, []);

% Check the given method.
if nargin < 2
    method = 'arcmin';
end
validatestring(method, {'arcmin', 'vectorsum'});

%% Compute mean.
switch lower(method)
    case 'arcmin'
        % Use Olson's method.
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
        m = wrapToPi(m(imin));
    case 'vectorsum'
        % Use the vector sum method.
        m = atan2(sum(sin(theta)), sum(cos(theta)));
end

end
