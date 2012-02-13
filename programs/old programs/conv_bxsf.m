function E=conv_bxsf(filespec,outfilename,rec_latt_type);
%conv bxsf file into .mat file with E data on cartesian grid
%filespec either single filename (can incl wildcards) or cell array of
%filenames
%filename must contain 'spinup' or 'spindn' to indicate spin polarized data
%this version needs to be told rec latt type and does interpolation to s.c.
%Cartesian grid on basis of this
%structure of output:

%spanning_vecs(3,3) [Angstrom^-1]
%rec_latt_vecs(3,3) [Angstrom^-1]
%dims(3) [size of data matrix in bxsf file]
%band(bandnum)
%FermiLevel [Ryd? from bxsf]
%cartE{bandnum}(:,:,:) [meV] interpolated energy
%cartX,Y,Z(:,:,:) k-space coords of points at which cartE{bandnum}(:,:,:) is given
%[Angstrom^-1]
%cuboidextents(2,3) [Angstrom^-1] extents of interpolated data
%cartdims(3) [size of interpolated output data]
%dx,dy,dz [Angstrom^-1] spacing of interpolated data
%DOS_E{bandnum} array of energies that DOS and IDOS are specified at
%DOS{bandnum} in states per eV per unit cell per spin
%IDOS{bandnum} in states per unit cell per spin
%filename{bandnum} converted from

Bohr_radius=0.529177; %Angstrom
numheaderlines=13;

if iscell(filespec)
    filenames=filespec;
    numbands=size(filenames,2);
else
    matchinginputfiles=dir(filespec);
    numbands=size(matchinginputfiles,1);
    [inputpath,name,ext,versn]=fileparts(filespec);
    for filenum=1:numbands
        filenames{filenum}=[inputpath matchinginputfiles(filenum).name];
    end
end

%read header info from all band files first
allbandsgood=true;
for bandnum=1:numbands    
    %assign spin dirn on basis of filename
    spindirns(bandnum)=0;
    if ~isempty(findstr('spinup',filenames{bandnum}))
        spindirns(bandnum)=1;
    elseif ~isempty(findstr('spindn',filenames{bandnum}))
        spindirns(bandnum)=2;
    end
    
    infile=fopen(filenames{bandnum});
    if infile==-1
        disp([filenames{bandnum} ' not found.']);
        allbandsgood=false;
        break;
    end        
    for n=1:numheaderlines
        headerline{n}=fgetl(infile);
    end
    fclose(infile);

    %get useful data from header
    band_FermiLevel(bandnum)=sscanf(headerline{2},'  Fermi Energy:     %f');
    band_dims{bandnum}=sscanf(headerline{8},'%f %f %f');
    band_rec_latt_vecs{bandnum}(1,:)=2*pi/Bohr_radius*sscanf(headerline{10},'%f %f %f'); %k in Angstrom^-1
    band_rec_latt_vecs{bandnum}(2,:)=2*pi/Bohr_radius*sscanf(headerline{11},'%f %f %f');
    band_rec_latt_vecs{bandnum}(3,:)=2*pi/Bohr_radius*sscanf(headerline{12},'%f %f %f');
    Wien2k_bandnums(bandnum)=sscanf(headerline{13},' BAND: %f');
    %check header info for different bands is compatible
    if bandnum>1
        if band_FermiLevel(bandnum)~=band_FermiLevel(bandnum-1) | any(band_rec_latt_vecs{bandnum}~=band_rec_latt_vecs{bandnum}) ...
                | any(band_dims{bandnum}~=band_dims{bandnum-1})
            disp(['Not all band files compatible.']);
            disp(['filename,FermiLevel,(rec_latt_vecs),(dims)']);
            disp([filenames{bandnum} ',' num2str(band_FermiLevel(bandnum)) ',(' num2str(reshape(band_rec_latt_vecs{bandnum},1,9)) '),' ...
                num2str(reshape(band_dims{bandnum},1,3)) ')']);
            allbandsgood=false;
            break;
        end
    end
end

if allbandsgood & numbands>0
    %sort bands in order of Wien2k_bandnum then spindirns
    [dummy order1]=sort(spindirns);
    [dummy order2]=sort(Wien2k_bandnums(order1));
    order=order1(order2);
    Wien2k_bandnums=Wien2k_bandnums(order);
    spindirns=spindirns(order);
    filenames=filenames(order);
    
    %set non band-specific vars
    dims=band_dims{1};
    FermiLevel=band_FermiLevel(1);
    rec_latt_vecs=band_rec_latt_vecs{1};   
    spanning_vecs(1,:)=rec_latt_vecs(1,:)/(dims(1)-1);
    spanning_vecs(2,:)=rec_latt_vecs(2,:)/(dims(2)-1);
    spanning_vecs(3,:)=rec_latt_vecs(3,:)/(dims(3)-1);    

    for bandnum=1:numbands
        %read E(k) data
        E{bandnum}=textread(filenames{bandnum},'%f ',dims(1)*dims(2)*dims(3),'headerlines',numheaderlines);
        E{bandnum}=reshape(E{bandnum},flipdim(dims,1)');
        %re-order dims to have row<=>a-axis, column<=>b-axis, 'page'<=>c-axis
        E{bandnum}=permute(E{bandnum},[3 2 1]);

        [cartE{bandnum} dx dy dz cuboidextents]=interp2Cart(E{bandnum},rec_latt_vecs,rec_latt_type);

        minX=0; minY=0; minZ=0;
        cartdims=size(cartE);
        
        %calculate DOS and integrated DOS
        [DOS_E{bandnum} DOS{bandnum} IDOS{bandnum}]=calcDOS(0.001,cartE{bandnum},1);
        disp(['Done ' num2str(bandnum) '/' num2str(numbands) ': ' filenames{bandnum}]);
    end
    [cartX cartY cartZ]=ndgrid(cuboidextents(1,1):dx:cuboidextents(2,1),cuboidextents(1,2):dy:cuboidextents(2,2),...
    cuboidextents(1,3):dz:cuboidextents(2,3));

    save([outfilename '.mat'],'spanning_vecs','rec_latt_vecs','dims','FermiLevel','cartX','cartY','cartZ',...
        'Wien2k_bandnums','spindirns','cartE','cuboidextents','cartdims','dx','dy','dz','DOS_E','DOS','IDOS','filenames');
end