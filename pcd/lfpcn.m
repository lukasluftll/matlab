function p = lfpcn(ls, lf)

ray = dir2cart(ls) .* repmat(ls.radius, 1, 3);
ray = ray(ls.ret,:);

for i = 1 : sum(ls.ret)
    [vi,t] = trav(tform2trvec(ls.sp(:,:,i)), ray(i,:), ...
        ls.xgv, ls.ygv, ls.zgv);
    vi = sub2ind(size(ls.data), vi(:,1), vi(:,2), vi(:,3));
    sum(ls.data(vi) .* diff(t));
end

end
