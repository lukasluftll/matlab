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
    
    properties ( SetAccess = private )
        data;
        xgv;
        ygv;
        zgv;
    end
    
    methods ( Access = public )
        % Construct a voxelmap object.
        function obj = voxelmap(data, xgv, ygv, zgv)
            % Check number of input arguments.
            if ~(nargin == 1 || nargin == 4)
                error('VOXELMAP takes 1 or 4 input arguments.')
            end
            
            % If the grid vectors are not given, construct them.
            if nargin == 1
                xgv = 0 : size(data, 1);
                ygv = 0 : size(data, 2);
                zgv = 0 : size(data, 3);
            end
            
            % Check the grid vectors.
            gvchk(xgv, ygv, zgv);
            
            % Store the input.
            obj.data = data;
            obj.xgv = xgv;
            obj.ygv = ygv;
            obj.zgv = zgv;
        end

        % Plot the map.
        function plot(obj, mode)
            % PLOT Plot voxel map.
            %   PLOT(VM, MODE) visualizes the voxelmap object VM using
            %   semi-transparent voxels. 
            %
            %   MODE can assume the following values:
            %
            %      'scale' (default) - The highest value of VM is drawn as
            %                          completely intransparent voxel, the 
            %                          lowest value of VM is drawn as
            %                          completely transparent voxel.
            %      'direct'          - Data values are directly interpreted 
            %                          as transparencies, with 0
            %                          representing complete intransparency
            %                          and 1 representing complete
            %                          transparency.
            
            %% Validate input.
            % Check the number of input arguments.
            narginchk(1, 2);
            
            % If there is mode given, set it.
            if nargin < 2
                mode = 'scale';
            end
            
            %% Plot map.
            switch lower(mode)
                case 'scale'
                    % Compute minimum and maximum data elements.
                    lim = [min(obj.data(:)), max(obj.data(:))];
                    
                    % Compute the transparency values to plot.
                    t = (obj.data + lim(1)) * diff(lim);
                case 'direct'
                    % Crop data to interval [0; 1].
                    t = obj.data;
                    t(t < 0) = 0;
                    t(t > 1) = 1;
                otherwise
                    error(['MODE ', mode, ' not supported.'])
            end
            
            % Plot.
            alphaplot(t, obj.xgv, obj.ygv, obj.zgv);
            axis equal
            grid on
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
    end
end

