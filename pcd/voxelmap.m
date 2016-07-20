classdef voxelmap
    % VOXELMAP Object for storing a 3D voxel map.
    %   M = VOXELMAP(DATA) creates a voxel map object.
    %   
    %   DATA is a IxJxK matrix that contains the map data.
    %
    %   M = VOXELMAP(DATA, XGV, YGV, ZGV) creates a voxel map object with a
    %   defined voxel grid.
    %   
    %   XGV, YGV, ZGV are vectors that define the rasterization of the 
    %   voxel grid. A voxel with index [i,j,k] contains all points [x,y,z]
    %   that satisfy the inequality:
    %
    %      (XGV(i) <= x < XGV(i+1))
    %      && (YGV(j) <= y < YGV(j+1)) 
    %      && (ZGV(k) <= z < ZGV(k+1))
    %
    %   DATA is a IxJxK matrix that contains the map data, where 
    %   I = numel(XGV)-1, J = numel(YGV)-1, and K = numel(ZGV)-1.
    %
    %   Example:
    %      data = [0.1, 0.2, 0.3; 0, 0, 0; 0, 0, 0];
    %      data = cat(3, data, zeros(3), 0.5 * ones(3));
    %      m = voxelmap(data)
    %
    %   VOXELMAP methods:
    %      plot - Visualize the map using transparency values.
    %
    %   See also LFMAP, REFMAP, DECAYMAP.
    
    % Copyright 2016 Alexander Schaefer
    
    properties
    end
    
    methods
    end
    
end

