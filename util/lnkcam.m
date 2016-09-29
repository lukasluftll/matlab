function lnkcam
% LNKCAM Link plot cameras.
%   LNKCAM synchronizes the cameras of all open figure windows. This way,
%   when you move the camera of one figure window, all other cameras move
%   as well.
%
%   Example:
%      surf(rand(5))
%      figure
%      surf(rand(5))
%      lnkcam
%
%   See also LINKPROP.

% Copyright 2016 Alexander Schaefer

%% Validate input and output.
nargoutchk(0,0)
narginchk(0,0)

%% Link all open axes.
% Find all open figures.
fig = findobj('type', 'figure');
ax = [];
for i = 1 : numel(fig)
    ax = [ax, findobj(fig(i), 'type', 'axes')]; %#ok<AGROW>
end

% Link the cameras of all axes.
for i = 2 : numel(ax)
    % Create the link.
    cameralink = linkprop([ax(1), ax(i)], {'CameraUpVector', ...
        'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
    
    % Store the link.
    setappdata(fig(1), 'camlnk', cameralink);
end

end
