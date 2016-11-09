function pclshow(varargin)

nargoutchk(0, 0)
narginchk(1, +Inf)
valfun = @(x) validateattributes(x, {'pointCloud'}, {}, '', 'VARARGIN');
cellfun(valfun, varargin)

dirname = [tempdir, 'pclshow/'];
if exist(dirname, 'dir') == 7
    delete([dirname, '*']);
else
    [success,msg,msgid] = mkdir(dirname);
    if ~success
        error(msgid, msg)
    end
end

pcwritefun = @(x) pcwrite(x, [tempname(dirname), '.pcd']);
cellfun(pcwritefun, varargin);

if isunix
    status = system(['unset LD_LIBRARY_PATH; pcl_viewer ', dirname, '*.pcd &']);
else
    status = system(['pcl_viewer ', dirname, '*.pcd &']);
end
if status ~= 0
    error('Failed to launch pcl_viewer.')
end

end
