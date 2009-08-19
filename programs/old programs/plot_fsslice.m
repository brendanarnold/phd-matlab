function plot_FSslice(normal,kpara,FS);
%makes contour plot of FS in plane perp. to normal at a point kpara along
%the direction of normal

%make a new grid tilted wrt X,Y,Z
[tgridX tgridY tgridZ]=ndgrid(linspace(-FS.maxX,FS.maxX,FS.ndims(1)),linspace(-FS.maxY,FS.maxY,FS.ndims(2)),...
    linspace(-FS.maxZ,FS.maxZ,FS.ndims(3)));
[tgridX tgridY tgridZ]=rotate3grid(tgridX,tgridY,tgridZ,normal);

tE=interp3(permute(FS.cartX,[2 1 3]),permute(FS.cartY,[2 1 3]),permute(FS.cartZ,[2 1 3]),...
    FS.cartE,reshape(tgridX,[],1),reshape(tgridY,[],1),reshape(tgridZ,[],1));
tE=reshape(tE,FS.ndivs,FS.ndivs,FS.ndivs);

plothandle=findobj('Tag','TiltedFermiSurfacePlot');
if ishandle(plothandle)
    %do nothing
else
    plothandle=figure('Tag','TiltedFermiSurfacePlot');    
end
plothandle;


contour(tE(:,:,1));
%isosurface(FS.cartY,FS.cartX,FS.cartZ,tE,FS.FermiLevel);
