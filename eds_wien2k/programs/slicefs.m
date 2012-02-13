function slice=sliceFS(normal,inplanedirns,kpara,dk,bandsdata,bandnums);
%function slice=sliceFS(normal,inplanedirns,kpara,dk,bandsdata,bandnums);
%
%returns 2d k-space slice of E(k) lin interpolated on to rect array
%plane perp. to normal at a point kpara along the direction of normal
%bandsdata.searchvol contains normals of volume (Cartesian) in which all of orbit must lie

normal=normal/norm(normal);
if isempty(bandnums)
    bandnums=1:length(bandsdata.Wien2k_bandnums);
end

kp=normal*kpara;
bplanes=bandsdata.searchvol(:,:);
[bkps bnrms]=normdir(bplanes);
plist=boundingpoly(kpara,normal,bkps,bnrms,bandsdata.ptinvol);
if isempty(plist)
    disp(['Slice plane does not intersect search volume']);
    hold on; plot3(kp(1),kp(2),kp(3),'^k');
end;

%express vertices of bounding polygon in in-plane coords rel. to point kp
for pnum=1:size(plist,1)
    planeplist(pnum,1)=dot(plist(pnum,:)-kp,inplanedirns{1});
    planeplist(pnum,2)=dot(plist(pnum,:)-kp,inplanedirns{2});
end
[inplanecoords.X inplanecoords.Y]=ndgrid(min(planeplist(:,1)):dk:max(planeplist(:,1)), ...
    min(planeplist(:,2)):dk:max(planeplist(:,2)));
if size(inplanecoords.X,1)<2
    debug=true;
end
slice.inplanedx=inplanecoords.X(2,1)-inplanecoords.X(1,1);
slice.inplanedy=inplanecoords.Y(1,2)-inplanecoords.Y(1,1);
[insidebounds onbounds]=inpolygon(inplanecoords.X,inplanecoords.Y,planeplist(:,1),planeplist(:,2));
slice.outsidebounds=~(insidebounds | onbounds);
slice.X=kp(1)+inplanecoords.X*inplanedirns{1}(1)+inplanecoords.Y*inplanedirns{2}(1);
slice.Y=kp(2)+inplanecoords.X*inplanedirns{1}(2)+inplanecoords.Y*inplanedirns{2}(2);
slice.Z=kp(3)+inplanecoords.X*inplanedirns{1}(3)+inplanecoords.Y*inplanedirns{2}(3);
slice.inplaneX=(slice.X-kp(1))*inplanedirns{1}(1)+(slice.Y-kp(2))*inplanedirns{1}(2)+(slice.Z-kp(3))*inplanedirns{1}(3);
slice.inplaneY=(slice.X-kp(1))*inplanedirns{2}(1)+(slice.Y-kp(2))*inplanedirns{2}(2)+(slice.Z-kp(3))*inplanedirns{2}(3);
wrappedcoords=wrapcoords(slice.X,slice.Y,slice.Z,bandsdata);

for bandnum=1:length(bandnums)
    slice.E{bandnum}=interp3(permute(bandsdata.cartX,[2 1 3]),permute(bandsdata.cartY,[2 1 3]),permute(bandsdata.cartZ,[2 1 3]),permute(bandsdata.cartE{bandnums(bandnum)},[2 1 3]),...
        reshape(wrappedcoords.X,[],1),reshape(wrappedcoords.Y,[],1),reshape(wrappedcoords.Z,[],1),'linear');
    slice.E{bandnum}=reshape(slice.E{bandnum},size(slice.X));
end
slice.kpara=kpara;

