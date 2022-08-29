
function surf2subnuclear(obj,hit,volview,transform)

%
% surf2subnuclear(obj,hit,volview,transform)
%
% Intended to be used as a callback for surface plots, which executes when
% a point is selected, updating the current_point of a SubNuclear volumeview 
% object. Right clicking will, in addition, add the point to the volume.
%
% To update points for a volumeview object, volview, upon clicking a graphical object having
% handle h in a figure fig, use:
% 
%   set(fig,'WindowButtonDownFcn',@(obj,hit)surf2subnuclear(obj,hit,volview,transform)); 
%
% and/or
%
%   set(h,'ButtonDownFcn',@(obj,hit)surf2subnuclear(obj,hit,volview,transform));
%
% where 'transform' is a function handle to the transform from the space of the plot 
% to the image voxel space. If no transform is given, then 
%
%   transform = @(x)volview.transforms(1).tr(x)
% 
% is assumed. 
%
% See also VOLUMEVIEW

% C. Kovach 2019

persistent plh
persistent numsurfpoints

if nargin <2
    if isvalid(plh)
        delete(plh)
    end
    return
    
end
if nargin < 4
    transform = @(x)snobj.transforms(1).itr(x);
end


if isfield(hit,'Button') && hit.Button > 1 ||  ismember('SelectionType',fieldnames(hit.Source))&& strcmp(hit.Source.SelectionType,'alt')
    if ~isempty(plh) && isvalid(plh)
        numsurfpoints=numsurfpoints+1;
        x = [plh.XData plh.YData plh.ZData];
        volview.addpoint(sprintf('surf. point %i',numsurfpoints),transform(x));
        pt = volview.points(end);
        plot3(x(1),x(2),x(3),'+','color',pt.plotcolor,'markersize',15,'linewidth',1);
        fprintf('\nAdding surface point %i',numsurfpoints);
    else
        fprintf('\nFirst left click to select a point')
    end
else
    if isempty(numsurfpoints)
        numsurfpoints = 0;
    end

    volview.current_point = transform(hit.IntersectionPoint);

    if isempty(plh) || ~isvalid(plh)

        hold on
        plh = plot3(hit.IntersectionPoint(1),hit.IntersectionPoint(2),hit.IntersectionPoint(3),'r+','markersize',20,'linewidth',2,'ButtonDownFcn',@(obj,hit)surf2subnuclear(obj,hit,volview,transform));

    else

        plh.XData = hit.IntersectionPoint(1);
        plh.YData = hit.IntersectionPoint(2);
        plh.ZData = hit.IntersectionPoint(3);

    end
end
