function plot_FSslice2(normal,kpara,FS,bandnums,justFermicontour);
%makes contour plot of FS in plane perp. to normal at a point kpara along
%the direction of normal

if isempty(bandnums)
    bandnums=1:length(FS.Wien2k_bandnums);
end

%generate coordinates in slice plane by rotating the z=kpara plane such
%that a vector in the z direction becomes parallel to normal
[sliceplaneX sliceplaneY sliceplaneZ]=ndgrid(linspace(2*FS.cuboidextents(1,1)-FS.cuboidextents(2,1),...
    2*FS.cuboidextents(2,1)-FS.cuboidextents(1,1),200),linspace(2*FS.cuboidextents(1,2)-FS.cuboidextents(2,2),...
    2*FS.cuboidextents(2,2)-FS.cuboidextents(1,2),200),kpara);
[sliceplaneX sliceplaneY sliceplaneZ]=rotate3grid(sliceplaneX,sliceplaneY,sliceplaneZ,normal,-1);
wrappedcoords=wrapcoords(sliceplaneX,sliceplaneY,sliceplaneZ,FS);

figure;
for bandnum=1:length(bandnums)
    Eslice=interp3(permute(FS.cartX,[2 1 3]),permute(FS.cartY,[2 1 3]),permute(FS.cartZ,[2 1 3]),permute(FS.cartE{bandnums(bandnum)},[2 1 3]),...
        reshape(wrappedcoords.X,[],1),reshape(wrappedcoords.Y,[],1),reshape(wrappedcoords.Z,[],1));
    Eslice=reshape(Eslice,size(sliceplaneX));

    if ~justFermicontour
        contour(Eslice);
    end
    hold on;
    colchars='rgbkymcw';
    colchar=colchars(mod(FS.Wien2k_bandnums(bandnums(bandnum)),7)+1);
    switch FS.spindirns(bandnums(bandnum))
        case 1      
            stylechar='-';
        case 2
            stylechar='-.';
    end
    contour(Eslice,FS.FermiLevel(FS.spindirns(bandnums(bandnum))),[stylechar colchar],'LineWidth',1);
end


%contourslice(permute(FS.cartX,[2 1 3]),permute(FS.cartY,[2 1 3]),permute(FS.cartZ,[2 1 3]),permute(FS.cartE,[2 1 3]),...
%    permute(wrappedcoords.X,[2 1 3]),permute(wrappedcoords.Y,[2 1 3]),permute(wrappedcoords.Z,[2 1 3]),[FS.FermiLevel FS.FermiLevel]);

