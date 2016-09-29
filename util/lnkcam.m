function lnkcam
    fig = findobj('type','figure');
    ax = [];
    for i = 1 : numel(fig)
        ax = [ax, findobj(fig(i),'type','axes')]; %#ok<AGROW>
    end

    % Link the cameras of the figures together.
    for i = 2 : numel(ax)
        cameralink = linkprop([ax(1), ax(i)], {'CameraUpVector', ...
        'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
        setappdata(fig(1), 'StoreTheLink', cameralink);
    end
end