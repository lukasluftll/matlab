function sizechk(a, b)

if ndims(a)~=ndims(b) || any(size(a)~=size(b))
    error(['Sizes of ', inputname(a), ' and ', inputname(b), ...
        ' do not match.'])
end

end
