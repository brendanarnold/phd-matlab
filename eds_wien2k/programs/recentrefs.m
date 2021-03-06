function newFS=recentreFS(FS,kcentre);
%function newFS=recentreFS(FS,kcentre);
%
%generates new k-grid centred on kcentre and linearly interpolates cartE into
%this. output is same size and (cuboid) shape as input, just translated.
%assumes input cuboid is tileable

currkcentre=0.5*[FS.cuboidextents(1,1)+FS.cuboidextents(2,1) ...
    FS.cuboidextents(1,2)+FS.cuboidextents(2,2) ...
    FS.cuboidextents(1,3)+FS.cuboidextents(2,3)];

trans_vec=kcentre-currkcentre;

newFS=FS;
newFS.cuboidextents=FS.cuboidextents+[trans_vec; trans_vec];
[newFS.cartX newFS.cartY newFS.cartZ]=ndgrid(newFS.cuboidextents(1,1):FS.dx:newFS.cuboidextents(2,1),...
    newFS.cuboidextents(1,2):FS.dy:newFS.cuboidextents(2,2),...
    newFS.cuboidextents(1,3):FS.dz:newFS.cuboidextents(2,3));
wrappedcoords=wrapcoords(newFS.cartX,newFS.cartY,newFS.cartZ,FS);

for bandnum=1:length(FS.Wien2k_bandnums)
    newFS.cartE{bandnum}=interp3(permute(FS.cartX,[2 1 3]),permute(FS.cartY,[2 1 3]),permute(FS.cartZ,[2 1 3]),...
        permute(FS.cartE{bandnum},[2 1 3]),permute(wrappedcoords.X,[2 1 3]),permute(wrappedcoords.Y,[2 1 3]),permute(wrappedcoords.Z,[2 1 3]));
    newFS.cartE{bandnum}=permute(newFS.cartE{bandnum},[2 1 3]);
end