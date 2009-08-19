function E=plot_FS(filename);

if ~exist([filename '.mat'])
    ['Need to make ' filename '.mat using conv_bxsf']
else

FS=load([filename '.mat']);
%plot surface
plothandle=findobj('Tag','FermiSurfacePlot');
if ishandle(plothandle)
    a=1;
    %do nothing
else
    plothandle=figure('Tag','FermiSurfacePlot');    
end
plothandle;

p=patch(isosurface(permute(FS.cartX,[2 1 3]),permute(FS.cartY,[2 1 3]),permute(FS.cartZ,[2 1 3]),...
    permute(FS.cartE,[2 1 3]),FS.FermiLevel));
isonormals(permute(FS.cartX,[2 1 3]),permute(FS.cartY,[2 1 3]),permute(FS.cartZ,[2 1 3]),permute(FS.cartE,[2 1 3]),p);
set(p,'FaceColor','red','EdgeColor','none');
camlight;

p2=patch(isosurface(permute(FS.cartX-FS.rec_latt_vecs(1,1),[2 1 3]),permute(FS.cartY-FS.rec_latt_vecs(1,2),[2 1 3]),...
    permute(FS.cartZ-FS.rec_latt_vecs(1,3),[2 1 3]),permute(FS.cartE,[2 1 3]),FS.FermiLevel));
isonormals(permute(FS.cartX-FS.rec_latt_vecs(1,1),[2 1 3]),permute(FS.cartY-FS.rec_latt_vecs(1,2),[2 1 3]),...
    permute(FS.cartZ-FS.rec_latt_vecs(1,3),[2 1 3]),permute(FS.cartE,[2 1 3]),p2);
set(p2,'FaceColor','green','EdgeColor','none');
lighting gouraud; 
camlight;

xlabel('X'); ylabel('Y'); zlabel('Z'); 

%make plane perp to Z
%hold on;
%hcontsurf=surf(linspace(-maxX,maxX,ndivs),linspace(-maxY,maxY,ndivs),...
%    zeros(ndivs));
%rotate(hcontsurf,[1 1 0],48.24);
%contsurfX=get(hcontsurf,'XData'); contsurfY=get(hcontsurf,'YData'); contsurfZ=get(hcontsurf,'ZData');
%delete(hcontsurf);
%hold on;
%contourslice(cartY,cartX,cartZ,cartE,contsurfX,contsurfY,contsurfZ);

%slicedata=interp3(cartY,cartX,cartZ,cartE,contsurfX,contsurfY,contsurfZ);



end
