function d = pccompare(pc, pcref)
% PCCOMPARE Compute distance between point clouds.
%   D = PCCOMPARE(PC, PCREF) computes for each point of point cloud PC the
%   Euclidean distance to the nearest point of PCREF.
%
%   PC and PCREF are pointCloud objects.
%
%   D is an N-element row vector. Its n-th element contains the distance
%   between the n-th point of PC and the nearest point of PCREF.
%
%   Example:
%      pc = pointCloud(rand(3));
%      pcref = pointCloud(rand(3));
%      d = pccompare(pc, pcref)
%
%   See also POINTCLOUD, KNNSEARCH.

% Copyright 2016 Alexander Schaefer

%% Validate input.
nargoutchk(0, 1)
narginchk(2, 2)
validateattributes(pcref, {'pointCloud'}, {'scalar'}, '', 'PCREF')
validateattriubtes(pc, {'pointCloud'}, {'scalar'}, '', 'PC')

%% Compute point-to-point distances.
[~,d] = knnsearch(pcref.Location, pc.Location);

end

