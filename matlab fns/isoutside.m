function outmatrix=isoutside(kparas,plane_normals,incl_point,inpoints,tol_r);
%function outmatrix=isoutside(kparas,plane_normals,incl_point,inpoints,tol_r);
%
%returns matrix same size as inpoints containing dist of each point from
%volume defined by plane_normals. incl_point is a point inside volume
%zero indicates point is inside
%tol_r is tolerance within which a point just outside is considered inside

npts=size(inpoints,1);
if npts==0
    return;
end

%this expands list of inpoints to include
if ~isempty(tol_r)
    offsets=[0 0 tol_r; 0 tol_r 0; tol_r 0 0; 0 0 -tol_r; 0 -tol_r 0; -tol_r 0 0; ...
        tol_r/1.73 tol_r/1.73 tol_r/1.73; -tol_r/1.73 tol_r/1.73 tol_r/1.73; tol_r/1.73 -tol_r/1.73 tol_r/1.73; tol_r/1.73 tol_r/1.73 -tol_r/1.73; ...
        tol_r/1.73 -tol_r/1.73 -tol_r/1.73; -tol_r/1.73 tol_r/1.73 -tol_r/1.73; -tol_r/1.73 -tol_r/1.73 tol_r/1.73; -tol_r/1.73 -tol_r/1.73 -tol_r/1.73;];
    orig_inpoints=inpoints;
    offset_inpoints(:,:)=repmat(inpoints,size(offsets,1),1)+kron(offsets,ones(npts,1));
    inpoints=cat(1,inpoints,offset_inpoints);
end
    
distmatrix=zeros(size(inpoints,1),size(plane_normals,1));
outmatrix=zeros(size(inpoints,1),1);
for plnum=1:size(plane_normals,1)
    %define matrix of intersections of lines from incl_point to vertices
    %with bounding planes      
    warning off MATLAB:divideByZero;
    IntMatrix(:,plnum)=((kparas(plnum)*plane_normals(plnum,1)-incl_point(1))*kparas(plnum)*plane_normals(plnum,1)+...
        (kparas(plnum)*plane_normals(plnum,2)-incl_point(2))*kparas(plnum)*plane_normals(plnum,2)+...
        (kparas(plnum)*plane_normals(plnum,3)-incl_point(3))*kparas(plnum)*plane_normals(plnum,3))./...
        ((inpoints(:,1)-incl_point(1))*kparas(plnum)*plane_normals(plnum,1)+...
        (inpoints(:,2)-incl_point(2))*kparas(plnum)*plane_normals(plnum,2)+...        
        (inpoints(:,3)-incl_point(3))*kparas(plnum)*plane_normals(plnum,3));
    warning on MATLAB:divideByZero;
    %value=0 indicates badly chosen incl_point, 0<value<1 indicates
    %point is on other side of plane, other values indicate point inside
    %plane    
end
for pnum=1:size(inpoints,1)
    absdists(pnum,1)=norm(inpoints(pnum,:)-incl_point);
end
outside=any(IntMatrix>0 & IntMatrix<1,2);

if ~isempty(tol_r)    
    indoffsets=(0:size(offsets,1))*npts;
    for n=1:npts
        outsidetol(n)=all(outside(n+indoffsets));
    end
    outside=outsidetol;
end

outmatrix=outside;
%the following distance calc is wrong
%if length(find(outside))>0
%    distmatrix(outside,:)=repmat(absdists(outside),1,size(plane_normals,1)).*(1-IntMatrix(outside,:));
%    distmatrix(distmatrix<=0)=nan;
%    outmatrix(outside)=min(distmatrix(outside,:),[],2);    
%end

    