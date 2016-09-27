classdef elevationmap
    % ELEVATIONMAP Object for storing an elevation data grid.
    %   EM = ELEVATIONMAP(PC, RES) rasterizes the x-y plane with grid 
    %   resolution RES and sets the elevation of each grid tile to the 
    %   maximum z coordinate of all points of point cloud PC that belong
    %   to this tile.
    %
    %   PC is a pointCloud object. RES is a positive scalar.
    %
    %   Example:
    %      em = elevationmap(pcread('teapot.ply'), 0.1)
    %
    %   ELEVATIONMAP methods:
    %   EVAL      - Elevation at x-y coordinate
    %   LIMITS    - Extension of map
    %   DIFF      - Disparity between elevation map and point cloud
    %   PLOT      - Visualize elevation map
    %
    %   See also POINTCLOUD.
    
    % Copyright 2016 Alexander Schaefer
    
    properties ( SetAccess = private )
        % Minimum x and y coordinates of the map; 1x2 vector.
        support
        
        % Extension of the map in numbers of map tiles in x and 
        % y direction; 1x2 vector.
        extension
        
        % Edge length of each tile; scalar.
        resolution
        
        % Elevation data of the map; IxJ matrix with I and J being the
        % number of tiles in x and y direction.
        elevation
    end
    
    methods ( Access = private )
        function i = idx(obj, point)
            nargoutchk(0, 1)
            narginchk(2, 2)

            d = point(:,1:2) - repmat(obj.support, size(point,1), 1);
            i = floor(d / obj.resolution) + 1;
            i = sub2ind(obj.extension, i(:,1), i(:,2));
        end                
    end
    
    methods ( Access = public )
        function obj = elevationmap(pc, res)
            nargoutchk(0, 1)
            narginchk(2, 2)
            
            validateattributes(pc, {'pointCloud'}, {'numel',1}, '', 'PC')            
            validateattributes(pc.XLimits, {'numeric'}, {'finite'}, '', 'PC.XLIMITS')
            validateattributes(pc.YLimits, {'numeric'}, {'finite'}, '', 'PC.YLIMITS')
            validateattributes(res, {'numeric'}, {'scalar', 'nonnegative'}, '', 'RES')

            obj.resolution = res;
            
            obj.support = floor([pc.XLimits(1),pc.YLimits(1)]/obj.resolution)*obj.resolution;
            
            obj.extension = floor(([pc.XLimits(2),pc.YLimits(2)]-obj.support)/obj.resolution) + 1;
            
            obj.elevation = ones(obj.extension) * min([0, pc.ZLimits(1)]);
            
            pc = removeInvalidPoints(pc);
            point = reshape(pc.Location(:), pc.Count, 3, 1);
            ie = obj.idx(point);
            for i = 1 : size(point,1)
                obj.elevation(ie(i)) = max(obj.elevation(ie(i)), point(i,3));
            end
        end
        
        function d = diff(obj, pc)
            nargoutchk(0, 1)
            narginchk(2, 2)
            
            validateattributes(pc, {'pointCloud'}, {'numel',1}, '', 'PC')
            
            pc = removeInvalidPoints(pc);
            point = reshape(pc.Location(:), pc.Count, 3, 1);
            ie = obj.idx(point);
            d = 0;
            for i = 1 : size(point,1)
                d = d + abs(point(i,3)-obj.elevation(ie(i)));
            end
            d = d / size(point,1);
        end 
        
        function plot(obj)
            xgv = obj.support(1) + kron(0:obj.extension(1), [1,1]) * obj.resolution;
            ygv = obj.support(2) + kron(0:obj.extension(2), [1,1]) * obj.resolution;
            
            [x, y] = ndgrid(xgv(2:end-1), ygv(2:end-1));
            
            z = kron(obj.elevation, ones(2));
            
            surf(x, y, z, 'EdgeColor', 'none')
            shading interp
            labelaxes
            axis equal
            grid
        end
    end 
end
