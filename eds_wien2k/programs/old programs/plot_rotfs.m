function plot_rotFS(normal,kpara,FS);
%makes contour plot of FS in plane perp. to normal at a point kpara along
%the direction of normal

%make a new grid of original data points but whose Z' axis is parallel to
%'normal' i.e. tgrid is a non-meshgrid object because coordinate system is
%not parallel to a matrix dimension.
[tgridX tgridY tgridZ]=ndgrid(linspace(-FS.maxX/2,FS.maxX/2,FS.cartdims(1)),linspace(-FS.maxY/2,FS.maxY/2,FS.cartdims(2)),...
    linspace(-FS.maxZ/2,FS.maxZ/2,FS.cartdims(3)));
[tgridX tgridY tgridZ]=rotate3grid(tgridX,tgridY,tgridZ,normal);

FS.cartX=FS.cartX-FS.maxX/2; FS.cartY=FS.cartY-FS.maxY/2; FS.cartZ=FS.cartZ-FS.maxZ/2;
tE=interp3(permute(FS.cartX,[2 1 3]),permute(FS.cartY,[2 1 3]),permute(FS.cartZ,[2 1 3]),permute(FS.cartE,[2 1 3]),...
    reshape(tgridX,[],1),reshape(tgridY,[],1),reshape(tgridZ,[],1));
tE=reshape(tE,size(FS.cartE));
isosurface(permute(FS.cartX,[2 1 3]),permute(FS.cartY,[2 1 3]),permute(FS.cartZ,[2 1 3]),permute(tE,[2 1 3]),FS.FermiLevel);
