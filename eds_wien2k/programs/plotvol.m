function hpatch=plotvol(haxes,bkps,bnormals,inclpoint,varargin);
%function hpatch=plotvol(haxes,bkps,bnormals,inclpoint,varargin);
%
%plot volume bounded by bplanes as semi-transparent wiremesh (e.g. BZ)


if ishandle(haxes)
    axes(haxes);
else
    figure;
end

if isempty(bkps)
    [bkps bnormals]=normdir(bnormals);
end
[BZfaces BZvertices]=vol_vertices(bkps,bnormals,inclpoint);      
hpatch=patch('Faces',BZfaces,'Vertices',BZvertices);
%set(hpatch,'FaceColor',[0 1 0],'FaceAlpha',0.07,'EdgeColor',[0 0 0]);
set(hpatch,'FaceColor', 'none','EdgeColor',[0 0 0], 'LineWidth', 2);
if ~isempty(varargin)
    set(hpatch, varargin{:});
end
