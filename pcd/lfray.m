function [p, L] = lfray(ls, lf, por)

narginchk(2, 4)

if ~isa(ls, 'laserscan')
    error('LS must be a laserscan object.')
end
if ~isa(lf, 'voxelmap')
    error('LF must be a voxelmap object.')
end

if nargin < 3
    por = 0;
end

if por < 0 || por > 1
    error('PNIR must stay in [0;1].')
end

iret = ls.ret;

ray = dir2cart(ls, iret);

p = ones(ls.count, 1);
L = zeros(ls.count, 1);
for i = 1 : ls.count
    if iret(i)
        start = tform2trvec(ls.sp(:,:,i)) + ray .* repmat(ls.rlim(1), 1, 3);
        [iv,t] = trav(start, ray .* repmat(diff(ls.rlim),1,3), ...
            ls.xgv, ls.ygv, ls.zgv);
        
        iv = sub2ind(size(ls.data), iv(:,1), iv(:,2), iv(:,3));
        
        pir = sum(ls.data(iv).*diff(t)*ls.radius(i)) / (1-por);

        L(i) = log(lf.data(iv(end)) / pir);
    else
        p(i) = por;
    end
end

end
