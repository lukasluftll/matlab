function t = hrtchk(tform)

if ~ismatrix(tform) || any(size(tform) ~= 4)
    error('TFORM must be of size 4x4.')
end

if ~ishrt(tform)
    error('TFORM must be a 4x4 homogeneous coordinate transform.')
end

end
