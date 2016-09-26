classdef elevationmap
    
    properties ( SetAccess = private )
        sup
        ext
        res
        ele
    end
    
    methods ( Access = private )
        function i = idx(obj, point)
            nargoutchk(0, 1)
            narginchk(2, 2)

            d = point(:,1:2) - repmat(obj.sup, size(point,1), 1);
            i = floor(d / obj.res) + 1;
            i = sub2ind(size(obj.ele), i(:,1), i(:,2));
        end                
    end
    
    methods ( Access = public )
        function obj = elevationmap(pc, res)
            nargoutchk(0, 1)
            narginchk(2, 2)
            
            validateattributes(pc, {'pointCloud'}, {'numel',1}, '', 'PC')            
            validateattributes(res, {'numeric'}, {'scalar', 'nonnegative'}, '', 'RES')
            
            obj.res = res;
            
            obj.sup = floor([pc.XLimits(1),pc.YLimits(1)]/obj.res)*obj.res;
            
            obj.ext = floor(([pc.XLimits(2),pc.YLimits(2)]-obj.sup)/obj.res) + 1;
            
            obj.ele = ones(size(obj.ext)) * min([0, pc.ZLimits(1)]);
            
            point = reshape(pc.Location(:), pc.Count, 3, 1);
            for i = 1 : pc.Count
                ip = obj.idx(point(i,:));
                obj.ele(ip) = max(obj.ele(ip), point(i,3));
            end
        end
        
        function plot(obj)
            xgv = obj.sup(1) + (1:obj.ext(1)) * obj.res;
            ygv = obj.sup(2) + (1:obj.ext(2)) * obj.res;
            [x, y] = ndgrid(xgv, ygv);
            surf(x, y, obj.ele);
        end
    end
    
end
