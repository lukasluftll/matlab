function map = pcmap(path, step)

narginchk(1, 2)

if ~isdir(path)
    error('PATH does not refer to a folder.')
end

if nargin < 2
    step = 0;
end

if step < 0
    error('STEP must be positive.')
end

file = dir([path, '/*.pcd']);

pcd = pcdread([path, '/', file(1).name]);
map = pointCloud([pcd.x(:), pcd.y(:), pcd.z(:)]);

for i = 2 : numel(file)
    pcd = pcdread([path, '/', file(i).name]);
    pc = pointCloud([pcd.x(:), pcd.y(:), pcd.z(:)]);
    
    if step == 0
        pcmerge(map, pc);
    else
        pcmerge(map, pc, step)
    end
end

end
