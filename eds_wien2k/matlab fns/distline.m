function dist=distline(r0,dirn,pts);
%dist2line(r0,dirn,pts);
%
%returns distance from points pts to line defined by r=r0+lamda*dirn

dirn=dirn/norm(dirn);
for ptnum=1:size(pts,1);
    pt=pts(ptnum,:);
    dist(ptnum,:)=norm(pt-(r0+dot(pt-r0,dirn)*dirn));
end
