function cuboid(vol, varargin)
% CUBOID Plot a 3D axis-aligned cuboid.
%   CUBOID(VOL) plots axis-aligned cuboids with the limits specified 
%   by VOL.
%
%   VOL is a Nx6 matrix. Its n-th row indicates the limits of the n-th
%   cuboid: 
%      VOL(n,:) = [xmin, ymin, zmin, xmax, ymax, zmax].
%
%   CUBOID(VOL, VARARGIN) plots a cuboid with the properties indicated by 
%   the name-value pair arguments VARARGIN. 
%   For possible name-value pairs, see the documentation of PATCH.
%
%   Example:
%      vol = [0, 2, 3, 4, 15, 11; -5, -5, -5, -3, -3, -3];
%      cuboid(vol, 'FaceColor', 'red', 'FaceAlpha', 0.5);
%
%   See also PATCH, CAT.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check if the user provided the required input argument.
narginchk(1, inf)

% Make sure the limits matrix has the right size.
if size(vol, 2) ~= 6
    error('VOL must have exactly 6 columns.')
end

%% Calculate the cube vertices and faces.
% Define the indices of VOL that form the vertices of the first cuboid.
vi = [1, 2, 3;
    4, 2, 3;
    4, 5, 3;
    1, 5, 3;
    1, 2, 6;
    4, 2, 6;
    4, 5, 6;
    1, 5, 6];

% Compute the indices of VOL that form the vertices of all cuboids.
vi = repmat(vi, size(vol, 1), 1) + kron((0:size(vol, 1)-1)', 6*ones(8, 3));

% Compute the vertices.
volT = vol'; 
vertices = volT(vi);

% Define the combinations of the rows of the vertices that form the faces 
% of the first cuboid.
faces = [1, 2, 3, 4; ...
    1, 2, 6, 5; ...
    2, 3, 7, 6; ...
    3, 4, 8, 7; ...
    4, 1, 5, 8; ...
    5, 6, 7, 8];

% Compute the combinations of the rows of the vertices that form the faces
% of all cuboids.
faces = repmat(faces, size(vol, 1), 1) ...
    + kron((0:size(vol, 1)-1)', 8*ones(6, 4));

%% Plot all faces.
patch('Faces', faces, 'Vertices', vertices, varargin{:});

end
