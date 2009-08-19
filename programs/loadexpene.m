function energy=loadexpene(fexpoutfilename);
%function energy=loadexpene(fexpoutfilename);
%loads Fourier interpolation coeffs from file produced by fexp fortran prog from wien2k output

fid=fopen(fexpoutfilename,'r');

hdr=fread(fid,1,'int32');
info=fread(fid,4,'int32');
hdr=fread(fid,1,'int32');
nx=info(1); ny=info(2); nz=info(3); nbands=info(4);

energy=zeros(nx,ny,nz,nbands);

disp(['Loading Fourier expanded energy bands ' num2str(nx) 'x' num2str(ny) 'x' num2str(nz)]);
for zn=1:nz
    for yn=1:ny
        for xn=1:nx
            hdr=fread(fid,1,'int32');
            energy(xn,yn,zn,:)=fread(fid,nbands,'float64');
            hdr=fread(fid,1,'int32');
        end
    end
end
disp(['Done']);

