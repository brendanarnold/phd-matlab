function newcoords=wrapcoords(inputX,inputY,inputZ,FS);
%wraps input coordinates to volume for which E is given in FS
%assumes cuboid volume in FS is tileable

newcoords.X=ones(size(inputX)).*(FS.cuboidextents(1,1))+mod(inputX-FS.cuboidextents(1,1),ones(size(inputX)).*(FS.cuboidextents(2,1)-FS.cuboidextents(1,1)));
newcoords.Y=ones(size(inputY)).*(FS.cuboidextents(1,2))+mod(inputY-FS.cuboidextents(1,2),ones(size(inputY)).*(FS.cuboidextents(2,2)-FS.cuboidextents(1,2)));
newcoords.Z=ones(size(inputZ)).*(FS.cuboidextents(1,3))+mod(inputZ-FS.cuboidextents(1,3),ones(size(inputZ)).*(FS.cuboidextents(2,3)-FS.cuboidextents(1,3)));
