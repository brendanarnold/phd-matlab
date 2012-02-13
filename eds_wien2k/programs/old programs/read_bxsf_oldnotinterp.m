function E=read_bxsf(filename);
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
E=reshape(E,dims');

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

plothandle=findobj('Tag','FermiSurfacePlot');
if ishandle(plothandle)
    %do nothing
else
    plothandle=figure;    
end
plothandle;


isosurface(kx,ky,kz,E,FermiLevel);

%E=circshift(E,floor(dims/2));

