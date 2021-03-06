function plist=boundingpoly(kpara,normal,bkparas,bnormals,incl_point);
%function plist=boundingpoly(kpara,normal,bkparas,bnormals,incl_point);
%
%returns a list of points in 3d space that are the vertices of intersection of a volume and a plane
%the volume is specified by planes defined by bkparas, bnormals, incl_point says what is interior of vol
%kpara, normal define the intersecting plane

%find coords of points in plane
%first calc r0's and d's of bounding lines i.e. intersection of bounding planes with intersection plane (see p55 lab book III)
bpt_cutoff=(5*max(abs(bkparas))); %don't consider intersections further away from incl_point than this
normal=normal/norm(normal); 
d=nan(numel(bkparas),3); r0=nan(numel(bkparas),3);
for planenum=1:numel(bkparas)
    bkpara=bkparas(planenum); bnormal=bnormals(planenum,:)/norm(bnormals(planenum,:));
    d(planenum,:)=cross(normal,bnormal);
    if norm(d(planenum,:))==0 | abs(kpara-bkpara)/norm(d(planenum,:))>bpt_cutoff
        %detects distant intersections of almost parallel planes
        continue;
    end
    r0(planenum,:)=((bkpara*dot(normal,bnormal)-kpara)*normal+(kpara*dot(bnormal,normal)-bkpara)*bnormal);
    denom=(dot(normal,bnormal)^2-1);
    if abs(denom)>=1e-30
        r0(planenum,:)=r0(planenum,:)/denom;
    else
        r0(planenum,:)=sign(r0(planenum,:))*sign(denom)*Inf;
    end
end
%now find intersection points of all pairs of lines
p=nan(numel(bkparas),numel(bkparas),3);
for planenum1=1:size(r0,1)
    for planenum2=planenum1+1:size(r0,1)
        lamda2=dot(r0(planenum1,:)-r0(planenum2,:),bnormals(planenum1,:));
        denom=dot(d(planenum2,:),bnormals(planenum1,:));
        if abs(denom)>=sqrt(realmin)
            lamda2=lamda2/denom;
        else
            lamda2=sign(lamda2)*sign(denom)*Inf;
        end
        p(planenum1,planenum2,:)=r0(planenum2,:)+lamda2*d(planenum2,:);
        if norm(permute(p(planenum1,planenum2,:),[1 3 2])-incl_point)>bpt_cutoff
            p(planenum1,planenum2,:)=nan;
        end
    end
end

%if true
%    tempplotlines(r0,d,p,kp);
%end
%eliminate intersections that are not vertices of bounding polygon
plist=reshape(p,numel(p)/3,3);
validpts=find(~any(isnan(plist) | isinf(plist),2));
outside=isoutside(bkparas,bnormals,incl_point,plist(validpts,:),1e-4);
plist(validpts(outside~=0),:)=nan;
plist=plist(~any(isnan(plist) | isinf(plist),2),:);

if size(plist,1)>0
%find angle within intersecting plane from incl_point to each p and sort points by this
dirn=d(find(~all(d==0,2) & ~any(isinf(d) | isnan(d),2),1),:);
dirn=dirn/norm(dirn);
for pnum=1:size(plist,1)
    cosangle(pnum)=dot(plist(pnum,:)-incl_point,dirn)/norm(plist(pnum,:)-incl_point);
    sinangle(pnum)=dot(plist(pnum,:)-incl_point,cross(dirn,normal))/norm(plist(pnum,:)-incl_point);
    angles(pnum)=atan2(sinangle(pnum),cosangle(pnum));
end
[sangles index]=sort(angles(:));
plist=plist(index,:);
end
