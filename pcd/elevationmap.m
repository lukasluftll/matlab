classdef elevationmap
    
    properties ( SetAccess = private )
        support
        extension
        resolution
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
