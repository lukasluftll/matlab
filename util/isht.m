function h = isht(tform)

sizeOk = size(tform) == [4, 4];
rotOk = tform(1:3,1:3)' == inv(tform(1:3,1:3));
rowOk = tform(4,:) == [0, 0, 0, 1];

h = all(isfinite(tform(:))) && sizeOk && rotOk && rowOk;

end

