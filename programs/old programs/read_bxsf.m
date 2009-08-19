function E=conv_bxsf(filename);
%conv bxsf file into .mat file with E data on cartesian grid
numheaderlines=13;

%read header info
infile=fopen(filename);
for n=1:numheaderlines
    headerline{n}=fgetl(infile);
end
fclose(infile);

%get useful data from header
FermiLevel=sscanf(headerline{2},'  Fermi Energy:     %f');
dims=sscanf(headerline{8},'%f %f %f');
spanning_vecs(1,:)=sscanf(headerline{10},'%f %f %f');
spanning_vecs(2,:)=sscanf(headerline{11},'%f %f %f');
spanning_vecs(3,:)=sscanf(headerline{12},'%f %f %f');
rec_latt_vecs(1,:)=spanning_vecs(1,:)*(dims(1)-1);
rec_latt_vecs(2,:)=spanning_vecs(2,:)*(dims(2)-1);
rec_latt_vecs(3,:)=spanning_vecs(3,:)*(dims(3)-1);

%read E(k) data
E=textread(filename,'%f ',dims(1)*dims(2)*dims(3),'headerlines',numheaderlines);
E=reshape(E,flipdim(dims,1)');
%re-order dims to have row<=>a-axis, column<=>b-axis, 'page'<=>c-axis
E=permute(E,[3 2 1]);

%calc Cartesian k for single rec latt unit cell
[a b c]=ndgrid(0:dims(1)-1,0:dims(2)-1,0:dims(3)-1);
kx=spanning_vecs(1,1).*a+spanning_vecs(2,1).*b+spanning_vecs(3,1).*c;
ky=spanning_vecs(1,2).*a+spanning_vecs(2,2).*b+spanning_vecs(3,2).*c;
kz=spanning_vecs(1,3).*a+spanning_vecs(2,3).*b+spanning_vecs(3,3).*c;

%make 2x2x2 tiling of rec latt unit cells
%first double-up cell || a
kx=[kx-rec_latt_vecs(1,1); kx(2:end,:,:)]; ky=[ky-rec_latt_vecs(1,2); ky(2:end,:,:)]; kz=[kz-rec_latt_vecs(1,3); kz(2:end,:,:)];
E=[E; E(2:end,:,:)];
%then along b
kx=[kx-rec_latt_vecs(2,1) kx(:,2:end,:)]; ky=[ky-rec_latt_vecs(2,2) ky(:,2:end,:)]; kz=[kz-rec_latt_vecs(2,3) kz(:,2:end,:)];
E=[E E(:,2:end,:)];
%finally along c
kx=cat(3,kx-rec_latt_vecs(3,1),kx(:,:,2:end)); ky=cat(3,ky-rec_latt_vecs(3,2),ky(:,:,2:end)); kz=cat(3,kz-rec_latt_vecs(3,3),kz(:,:,2:end));
E=cat(3,E,E(:,:,2:end));

%retain only every other element
%[kx,ky,kz,E]=reducevolume(kx,ky,kz,E,4);

%interpolate reciprocal lattice based array to cartesian array over bounding
%rectangle
maxX=max(reshape(kx,[],1)); maxY=max(reshape(ky,[],1)); maxZ=max(reshape(kz,[],1));
minX=min(reshape(kx,[],1)); minY=min(reshape(ky,[],1)); minZ=min(reshape(kz,[],1));
%find spacing dx,dy,dz of least dense Cartesian grid that includes all original data
%points
dxs=spanning_vecs(:,1); 
dx=min(abs(dxs(dxs~=0)));
dys=spanning_vecs(:,2); 
dy=min(abs(dys(dys~=0)));
dzs=spanning_vecs(:,3); 
dz=min(abs(dzs(dzs~=0)));
%note: only works if (ax,bx,cx)/min[(ax,bx,cx)~=0] are all integers

cartdims=round([(maxX-minX)/dx+1 (maxY-minY)/dy+1 (maxZ-minZ)/dz+1]);
[cartX cartY cartZ]=ndgrid(minX:dx:maxX,minY:dy:maxY,minZ:dz:maxZ);

starttime=clock;
cartE=interp3(permute(kx,[2 1 3]),permute(ky,[2 1 3]),permute(kz,[2 1 3]),permute(E,[2 1 3]),permute(cartX,[2 1 3]),...
    permute(cartY,[2 1 3]),permute(cartZ,[2 1 3]),'*linear');
cartE=permute(cartE,[2 1 3]);
%cartE=griddata3(reshape(kx,[],1),reshape(ky,[],1),reshape(kz,[],1),reshape(E,[],1),cartX,cartY,cartZ);
elapsedtime=etime(clock,starttime);
[filename ': interp from ' num2str(size(E)) ' to ' num2str(size(cartE)) ' took ' num2str(elapsedtime) 's.']

save([filename '.mat'],'spanning_vecs','rec_latt_vecs','dims','FermiLevel','cartX','cartY','cartZ',...
    'cartE','maxX','maxY','maxZ','cartdims');

