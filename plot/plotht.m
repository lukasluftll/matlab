function plotht(ht, l, varargin)
% PLOTHT Plot homogeneous transformations.
%   PLOTHT(HT) plots N Cartesian coordinate systems that represent the
%   rotation and translation of N 4x4 homogeneous transformation matrices.
%   These transformation matrices are the pages of the 4x4xN matrix HT. 
%
%   PLOTHT(HT, L) plots the axis vectors with length L.
%
%   PLOTHT(HT, L, VARARGIN) additionally specifies the plot linestyle for 
%   the axis vectors. Any marker in LINESPEC is drawn at the base instead 
%   of an arrow on the tip. Use a marker of '.' to specify no marker at 
%   all. See PLOT for available linestyles.
%
%   Example:
%      plotht(eye(4))
%
%   See also PLOT, QUIVER3.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(1, +Inf)

% Check it HT is a concatenation of homogeneous transformation matrices.
if ~all(ishrt(ht))
    error(['HT(:,:,' num2str(find(~ishrt(ht), 1)), ...
        ') must be a 4x4 homogeneous translation-rotation matrix.'])
end

% If the axis vector length is not given, define it.
l(nargin<2) = 1;

% Check the axis vector length argument.
validateattributes(l, {'numeric'}, {'scalar', 'positive'}, '', 'L')

%% Plot unit vectors.
% Compute the origins of the coordinate systems.
o = tform2trvec(ht);

% Compute the unit vectors.
e = pagetimes(tform2rotm(ht), repmat(eye(3), 1, 1, size(ht,3))) * l;
e = permute(e, [3,1,2]);

% Plot the unit vectors.
quiver3(o(:,1),o(:,2),o(:,3),e(:,1,1),e(:,2,1),e(:,3,1),0,'r',varargin{:})
hold on
quiver3(o(:,1),o(:,2),o(:,3),e(:,1,2),e(:,2,2),e(:,3,2),0,'g',varargin{:})
quiver3(o(:,1),o(:,2),o(:,3),e(:,1,3),e(:,2,3),e(:,3,3),0,'b',varargin{:})

end
