function gen2dklist(bands,outfilename,plane,bkps,bnormals,dkx,dky,removepts,showslice);
%function gen2dklist(bands,outfilename,plane,bkps,bnormals,dkx,dky,removepts,showslice);
%
%generates klist file for W2k with kpts on a single plane
%plane is a struct with fields:
%       .xdir,.ydir specifying the two in plane dirns in cartesian coords
%       .kp, shortest distance of plane from origin (Angstrom^-1)
%       .inclpt is a point inside bounded area given by bkps, bnormals
%bkps, bnormals define the set of planes that bound the klist plane. cartesian coords
%dkx,dky is mesh spacing to aim for in Angst^-1

%should work for P B C F not H
%other possible ways of defining plane (not implemented):
%       either cell array 1st element 1x3 normal of plane in rlv coords, 2nd element dist of plane from orig
%       or 3x3 array of points on plane in rlv coords
%       or >3x4 array of points that bound region (error if not on one plane)

denom=1e4;
if isfield(plane,'xdir')
    plane.xdir=plane.xdir./norm(plane.xdir);
    plane.ydir=plane.ydir./norm(plane.ydir);        
    plane.normal=cross(plane.xdir,plane.ydir);
    plane.normal=plane.normal./norm(plane.normal);
end

%this for klist in the conventional (not primitive) form
%make sure origin and spanning vectors are rational fractions of a*_x,b*_y,c*_z
span.x=floor(denom*dkx*plane.xdir./abs(diag(bands.rec_latt_vecs)'))./(denom./abs(diag(bands.rec_latt_vecs)'));
span.y=floor(denom*dky*plane.ydir./abs(diag(bands.rec_latt_vecs)'))./(denom./abs(diag(bands.rec_latt_vecs)'));
orig=floor(denom*plane.kp*plane.normal./abs(diag(bands.rec_latt_vecs)'))./(denom./abs(diag(bands.rec_latt_vecs)'));
denom=denom/hcf(floor(abs([span.x span.y orig]*denom)));
%work back from rounded values to get new plane
plane.xdir=span.x/norm(span.x);
plane.ydir=span.y/norm(span.y);
plane.kp=norm(orig.*abs(diag(bands.rec_latt_vecs)'));
plane.normal=cross(plane.xdir,plane.ydir);

plist=boundingpoly(plane.kp,plane.normal,bkps,bnormals,plane.inclpt);
if isempty(plist)
    disp(['Slice plane does not intersect search volume']);
    hold on; plot3(kp(1),kp(2),kp(3),'^k');    
elseif showslice
    plotvol([],bkps,bnormals,[0 0 0]); hold on;
    xlabel('x'); ylabel('y'); zlabel('z');
    plot3(plist(:,1),plist(:,2),plist(:,3),'+');
    
    vertices=[]; faces=[];
    %make patch showing plane
    for ptnum=1:size(plist,1)
        vertices(end+1,:)=plist(ptnum,:);
        if ptnum>size(faces,2) & size(faces,1)>0
            %bug in matlab: should pad faces with <max vert with nans but
            %this causes extra edges to appear
            faces=[faces faces(:,size(faces,2))];
        end
        faces(1,ptnum)=size(vertices,1);
    end    
    hpatch=patch('Faces',faces,'Vertices',vertices);
    set(hpatch,'FaceColor',[1 0 0],'FaceAlpha',0.3,'EdgeColor',[0 0 0]);
end;

%express vertices of bounding polygon in in-plane coords rel. to point kp
for pnum=1:size(plist,1)
    planeplist(pnum,1)=dot(plist(pnum,:)-plane.kp*plane.normal,plane.xdir);
    planeplist(pnum,2)=dot(plist(pnum,:)-plane.kp*plane.normal,plane.ydir);
end

%define 2d array of plane kpts in units of span.x span.y
planeplist=planeplist./repmat([norm(span.x) norm(span.y)],[size(planeplist,1) 1]);
[mesh.x mesh.y]=ndgrid(floor(min(planeplist(:,1))):ceil(max(planeplist(:,1))),floor(min(planeplist(:,2))):ceil(max(planeplist(:,2))));

%coordinates of mesh points in 3D cartesian frame
planecoords.x=plane.kp*plane.normal(1)+mesh.x*span.x(1)+mesh.y*span.y(1);
planecoords.y=plane.kp*plane.normal(2)+mesh.x*span.x(2)+mesh.y*span.y(2);
planecoords.z=plane.kp*plane.normal(3)+mesh.x*span.x(3)+mesh.y*span.y(3);

nkpts=numel(planecoords.x);
disp(['n_kpts in rectangle=' num2str(nkpts)]);
if removepts
    outmatrix=isoutside(bkps,bnormals,plane.inclpt,[planecoords.x(1:end)' planecoords.y(1:end)' planecoords.z(1:end)'],1e-4);
    planecoords.x=planecoords.x(~outmatrix);
    planecoords.y=planecoords.y(~outmatrix);
    planecoords.z=planecoords.z(~outmatrix);
    disp(['after removing kpts outside bounds, n_kpts=' num2str(sum(~outmatrix))]);
    nkpts=numel(planecoords.x);
end
    
if showslice
    plot3(planecoords.x,planecoords.y,planecoords.z,'.b');
end

klistcoord=mod(round(denom*[planecoords.x(1:end)' planecoords.y(1:end)' planecoords.z(1:end)']./repmat(diag(bands.rec_latt_vecs)',[numel(planecoords.x) 1])),denom);

weight=1;
fid=fopen(outfilename,'w');
for n=1:nkpts
    fprintf(fid,'%10i%10i%10i%10i%10i%5.1f',n,klistcoord(n,1),klistcoord(n,2),klistcoord(n,3),denom,weight);
    if n==1
        fprintf(fid,'%5.1f%5.1f    %6i k, div: (%3i%3i%3i)\r\n',-7.0,1.5,nkpts,size(mesh.x,1),size(mesh.x,2),1);
    else
        fprintf(fid,'\r\n');
    end
end
fprintf(fid,'END\r\n');
fclose(fid);
    