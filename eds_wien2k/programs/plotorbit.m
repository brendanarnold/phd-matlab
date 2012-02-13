function [newhaxes,handle]=plotorbit(haxes,vertices,orbitnum,params);

if ~any(ishandle(haxes))
    figure;
    haxes=axes;
else
    axes(haxes);
end

if ~any(ismember(fieldnames(params), 'LineWidth'))
    params.LineWidth = 1;
end
if ~any(ismember(fieldnames(params), 'Color'))
    params.Color = 'k';
end
if ~any(ismember(fieldnames(params), 'LineStyle'))
    params.LineStyle = '-';
end

hold on;
if isnumeric(vertices{orbitnum}) & ~isempty(vertices{orbitnum})
    handle=plot3(vertices{orbitnum}(:,1),vertices{orbitnum}(:,2),vertices{orbitnum}(:,3));
    if ~isempty(params)
        set(handle,'LineWidth',params.LineWidth,'Color',params.Color,'LineStyle',params.LineStyle);
    end
else
    handle=[];
end
newhaxes=haxes;
