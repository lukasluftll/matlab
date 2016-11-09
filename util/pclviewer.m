function pclviewer(varargin)
% PCLVIEWER Visualize point clouds with pcl_viewer.
%   PCLVIEWER(PC1, PC2, ...) starts the external application pcl_viewer to
%   display the pointCloud objects PC1, PC2, etc.
%
%   For this function to work, pcl_viewer must be installed on the system
%   and can be started by typing "pcl_viewer" into the console window.
%
%   Example:
%      pclviewer(pcread('teapot.ply'));
%
%   See also PCSHOW.

% Copyright 2016 Alexander Schaefer

%% Validate input and output.
% Validate the number of input and output arguments.
nargoutchk(0, 0)
narginchk(1, +Inf)

% Validate the input arguments.
valfun = @(x) validateattributes(x, {'pointCloud'}, {}, '', 'VARARGIN');
cellfun(valfun, varargin)

%% Write point clouds to file.
% Create an empty folder in a temporary directory.
dirname = [tempdir, 'pclshow/'];
if exist(dirname, 'dir') == 7
    delete([dirname, '*']);
else
    [success,msg,msgid] = mkdir(dirname);
    if ~success
        error(msgid, msg)
    end
end

% Write all point clouds to file.
pcwritefun = @(x) pcwrite(x, [tempname(dirname), '.pcd']);
cellfun(pcwritefun, varargin);

%% Start pcl_viewer.
cmd = ['pcl_viewer ', dirname, '*.pcd &'];
if isunix
    % If the operating system is Unix, MATLAB prepends its custom shared 
    % library paths to the default library path and in this way shadows 
    % some libraries on which pcl_viewer depends. However, the library 
    % versions MATLAB provides are not compatible with pcl_viewer. 
    % Therefore, unset the MATLAB library path to make the original 
    % libraries accessible again.
    status = system(['unset LD_LIBRARY_PATH; ', cmd]);
else
    status = system(cmd);
end
if status ~= 0
    error('Failed to launch pcl_viewer.')
end

end
