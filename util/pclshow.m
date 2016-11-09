function pclshow(varargin)

nargoutchk(0, 0)
narginchk(1, +Inf)
valfun = @(x) validateattributes(x, {'pointCloud'}, {}, '', 'VARARGIN');
cellfun(valfun, varargin)

dirname = [tempdir, 'pclshow/'];
if exist(dirname, 'dir') == 7
    delete([dirname, '*.pcd']);
else
    [success,msg,msgid] = mkdir(dirname);
    if ~success
        error(msgid, msg)
    end
end

pcwritefun = @(x) pcwrite(x, [tempname(dirname), '.pcd']);
cellfun(pcwritefun, varargin);

script = [dirname, 'pclshow'];
fid = fopen(script, 'w');
if fid < 0
    error('Failed to create script for starting pcl_viewer.')
end
fprintf(fid, '%s', ['pcl_viewer ', dirname, '*.pcd']);
fclose(fid);

status = system(['chmod ugo+x ' script]);
if status ~= 0
    error('Failed to make script for starting pcl_viewer executable.')
end
status = system(script);
if status ~= 0
    error('Failed to launch pcl_viewer via command line.')
end

end
