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
        function i = idx(obj, p)
            % IDX Compute tile index given x-y coordinates.
            %   I = IDX(OBJ, P) computes the linear index I of the tile to
            %   which the given x-y coordinates P belong.
            %
            %   P is an Mx2 vector.
            %
            %   I is an M-element column vector. If P(m,:) lies outside the
            %   map, I(m) is NaN.
            
            %% Validate input and output.
            nargoutchk(0, 1)
            narginchk(2, 2)

            %% Compute tile index.
            % Determine which coordinates lie inside the map.
            lim = obj.limits;
            rin = lim(1,1) <= p(:,1) & p(:,1) < lim(1,2) ...
                & lim(2,1) <= p(:,2) & p(:,2) < lim(2,2);
            
            % Set the tile indices corresponding to points outside the map
            % to NaN.
            i(~rin) = NaN;
            
            % Compute the indices corresponding to the valid coordinates.
            d = p(rin,:) - repmat(obj.support, sum(rin), 1);
            isub = floor(d / obj.resolution) + 1;
            i(rin) = sub2ind(obj.extension, isub(:,1), isub(:,2));
            
        end                
    end
    
    methods ( Access = public )
        function obj = elevationmap(pc, res)
            % ELEVATIONMAP Constructor.
            
            %% Validate input and output.
            % Check number of input and output arguments.
            nargoutchk(0, 1)
            narginchk(2, 2)
            
            % Check the validity of the point cloud.
            validateattributes(pc, {'pointCloud'}, {'numel',1}, '', 'PC')            
            validateattributes(pc.XLimits, {'numeric'}, {'finite'}, ...
                '', 'PC.XLIMITS')
            validateattributes(pc.YLimits, {'numeric'}, {'finite'}, ...
                '', 'PC.YLIMITS')
            
            % Check the validity of the resolution.
            validateattributes(res, {'numeric'}, ...
                {'scalar', 'nonnegative'}, '', 'RES')

            %% Compute property values.
            % Store the resolution.
            obj.resolution = res;
            
            % Compute the minimum x and y coordinates of the map.
            obj.support = floor([pc.XLimits(1),pc.YLimits(1)] / res) * res;
            
            % Compute the number of map tiles in x and y direction.
            obj.extension = floor(...
                ([pc.XLimits(2),pc.YLimits(2)]-obj.support) / res) + 1;
            
            % Set the elevation data to the minimum elevation.
            obj.elevation = ones(obj.extension) * min([0, pc.ZLimits(1)]);
            
            % Remove all NaN and infinite points from the point cloud.
            pc = removeInvalidPoints(pc);
            
            % Make sure the point cloud is represented by an Mx3 vector.
            point = reshape(pc.Location(:), pc.Count, 3, 1);
            
            % Compute the indices of the tiles to which each point belongs.
            ie = obj.idx(point);
            
            % Set the elevation of each tile to the maximum z coordinate of
            % all points belonging to this tile.
            for i = 1 : size(point,1)
                obj.elevation(ie(i))=max(obj.elevation(ie(i)),point(i,3));
            end
        end
        
        function l = limits(obj)
            nargoutchk(0, 1)
            s = repmat(obj.support, 2, 1);
            l = (s + [0,0; obj.extension*obj.resolution]).';
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
