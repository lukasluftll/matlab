classdef elevationmap
    
    properties ( SetAccess = private )
        x
        y
        cmin
        res
        elevation
    end
    
    methods ( Access = private )
        function varargout = idx(obj, point)
            narginchk(2, 2)
            nargoutchk(1, 2)
            
            point = point(all(isfinite(point), 2),:);
            
            switch nargout
                case 1
                    [ix,iy] = obj.idx(point);
                    varargout = {sub2ind(size(obj.x), ix, iy)};
                case 2
                    i = floor(point(:,1:2)-repmat(obj.cmin,size(point,1),1) / obj.res) + 1;
                    varargout = {i(:,1), i(:,2)};         
            end
        end                
    end
    
    methods ( Access = public )
        function obj = elevationmap(pc, res)
            validateattributes(pc, {'pointCloud'}, {'numel',1}, '', 'PC')            
            validateattributes(res, {'numeric'}, {'numel',1, 'nonnegative'}, '', 'RES')
            
            obj.res = res;
            
            obj.cmin = floor([pc.XLimits(1), pc.YLimits(1)] / obj.res) * obj.res;
            
            lim = floor([pc.XLimits; pc.YLimits] / res) + 0.5;
            [obj.x, obj.y] = ndgrid((lim(1,1):lim(1,2)) * obj.res, ...
                (lim(2,1):lim(2,2)) * obj.res);
            
            obj.elevation = zeros(size(obj.x));
            i = obj.idx(pc.Location);
            obj.elevation(i) = max(obj.elevation(i), pc.Location(i,3), 'omitnan');
        end
        
        function plot(obj)
            surf(obj.x, obj.y, obj.elevation);
            demcmap(obj.elevation);
        end
    end
    
end
