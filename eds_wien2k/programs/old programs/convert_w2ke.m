function convert_W2kE(pathcasename,bandnumlist,outfilename,rec_latt_type,isexpanded);
%function convert_W2kE(pathcasename,bandnumlist,outfilename,rec_latt_type,isexpanded);
%
%converts lapw1 output .energy file into cartesian array
%format of .energy:
%lines 1-4?
%columns of indexline: kx,ky,kz (fractions of rec latt vector), kpoint number, ?,
%number of bands, ?
%then list of band number, energy (Ryd)

%find if spin pol

Bohr_radius=0.529177; %Angstrom

if exist([pathcasename '.energyup'])==2 & exist([pathcasename '.energydn'])==2
    spinpol=true;
    dirnstr={'up' 'dn'};
else
    spinpol=false;
    dirnstr={''};
end


%get Fermi energy of final SCF cycle
%if false
for dirnnum=1:length(dirnstr)
    infile=fopen([pathcasename '.scf2' dirnstr{dirnnum}]);
    notatend=true; linestr=[];
    while notatend;
        prevlinestr=linestr;
        linestr=fgetl(infile);    
        if findstr(linestr,':FER')
            notatend=false;
            break;
        end
    end
    fclose(infile);
    FermiLevel(dirnnum)=str2num(linestr(findstr(linestr,'=')+1:end));

    disp('Reading IBZ energies for all bands...');
    dirnenergy=importenergy([pathcasename '.energy' dirnstr{dirnnum}]);
    if dirnnum==2 & (size(dirnenergy,1)~=size(energy,1))
        %minority spin might have fewer fully defined energy bands, if so discard corresponding
        %majority ones
        energy=energy(1:size(dirnenergy,1),:,:);
    end
    energy(:,dirnnum,:)=dirnenergy;
    numkpoints=size(energy,3);
end
energy(energy==0)=NaN;    
    
%find which bands cross Fermi Level
for dirnnum=1:length(dirnstr)
    FermiCrossing(:,dirnnum)=min(energy(:,dirnnum,:),[],3)<FermiLevel(dirnnum) & max(energy(:,dirnnum,:),[],3)>FermiLevel(dirnnum);
    FullyDefined(:,dirnnum)=sum(~isnan(energy(:,dirnnum,:)),3)==numkpoints;
end
disp(['... ' num2str(sum(FullyDefined)) ' fully defined bands read.']);
if spinpol
    disp(['Bands ' num2str(find(FermiCrossing(:,1)')) '(up), ' num2str(find(FermiCrossing(:,2)')) '(down) cross Fermi level.']);    
    numbands=sum(FullyDefined(:,1) & FullyDefined(:,2));
    if isempty(bandnumlist)
        %by default only keep bands that cross the Fermi surface
        bandnumlist=find(FermiCrossing(:,1) | FermiCrossing(:,2));
    end
else
    numbands=sum(FullyDefined(:,1));    
    disp(['Bands ' num2str(find(FermiCrossing(:,1)')) ' cross Fermi level.']);
    if isempty(bandnumlist)
        bandnumlist=find(FermiCrossing(:,1));
    end
end

disp(['Converting bands ' num2str(bandnumlist')]);

%end %skipped
%populate entire BZ from IBZ data

%read symm operations
disp('Reading symm operations, rec latt vecs from .outputkgen');
infile=fopen([pathcasename '.outputkgen'],'r');
numheaderlines=11;
for linenum=1:numheaderlines
    linestr=fgetl(infile);
end
if ~strcmp(linestr(5:26),'SYMMETRY MATRIX NR.  1')
    disp('.outputkgen not expected format');
end
notatend=true;
matgroupnum=1;
while notatend
    for rownum=1:3
        rowstr{rownum}=fgetl(infile);
        row(rownum,:)=sscanf(rowstr{rownum},'%f',12);
    end
    for colnum=1:4
        symmmat(:,:,(matgroupnum-1)*4+colnum)=[row(1,(colnum-1)*3+(1:3)); row(2,(colnum-1)*3+(1:3)); row(3,(colnum-1)*3+(1:3))];
    end
    linestr=fgetl(infile);
    if ~strcmp(linestr(5:19),'SYMMETRY MATRIX')
        notatend=false;
    else
        matgroupnum=matgroupnum+1;
    end
end    
numsymmmat=size(symmmat,3);
%rec latt vecs immediately follow
if ~strcmp(linestr(5:26),'G1        G2        G3')
    disp('Can''t find rec latt vecs');
else
    for vecnum=1:3
        linestr=fgetl(infile);
        %note vecs are in cols not rows
        rec_latt_vecs(:,vecnum)=2*pi/Bohr_radius*sscanf(linestr,'%f%f%f'); %result is in Angstrom^-1
    end
    %define matrix of rec latt vecs such that x-components are normalized
    %by a_x, y-components by b_y, z-components by c_z
    n_rec_latt_vecs=rec_latt_vecs./repmat(diag(rec_latt_vecs),1,3);
end

%.klist contains: dims; IBZ kpoint no., rec latt coord kx/(a_x/d),ky/(b_y/d),kz/(c_z/d), d=denominator,
%weight (i.e. no. of symmetry rel points incl IBZ one)
disp('Reading dims, Cartesian IBZ coords from .klist');
infile=fopen([pathcasename '.klist'],'r');
if infile==-1
    disp('File not found: case.klist');
    return;
end

%get dims from last columns of first line of .klist
firstlinestr=fgetl(infile);
startchar=findstr(firstlinestr,'(')+1;
endchar=findstr(firstlinestr,')')-1;
dims=sscanf(firstlinestr(startchar:endchar),'%f',3);
fclose(infile);
spanning_vecs=rec_latt_vecs./(repmat(dims,1,3));

%get IBZ kpoint coords
infile=fopen([pathcasename '.klist'],'r');
notatend=true; linenum=1;
while notatend
    linestr=fgetl(infile);
    if strcmp(linestr(1:3),'END')
        notatend=false;
        break;
    end
    IBZ_klist(linenum,:)=sscanf(linestr,'%f',6);
    linenum=linenum+1;
end
IBZ_numkpoints=linenum-1;
disp(['There are ' num2str(IBZ_numkpoints) ' point in .klist file']);

%get denominator of rec latt spanning vecs
%note this is only same as dims for cubic system
if length(unique(IBZ_klist(:,5)))~=1
    disp('Denominator changed?!');
else
    denominator=IBZ_klist(1,5);
end

%generate full unit. cell matrix
%convert inequiv kpoints to rec_latt_basis and then operate with symmetry
%matrices to give equiv points
disp('Generating k-mesh over whole unit cell');
unitcell_kpoints=zeros(0,4); BZ_kpointnum=1;
for IBZ_kpointnum=1:IBZ_numkpoints    
    for symmmatnum=1:numsymmmat     
        reclattbased_IBZpoint=(inv(n_rec_latt_vecs)*IBZ_klist(IBZ_kpointnum,2:4)');
        symmrel_kpoints(symmmatnum,1:4)=[IBZ_kpointnum (squeeze(symmmat(:,:,symmmatnum))*reclattbased_IBZpoint)'];        
    end
    %wrap back to first unit cell and discard duplicates
    symmrel_kpoints(:,2)=mod(symmrel_kpoints(:,2),dims(1));
    symmrel_kpoints(:,3)=mod(symmrel_kpoints(:,3),dims(2));
    symmrel_kpoints(:,4)=mod(symmrel_kpoints(:,4),dims(3));
    symmrel_kpoints=unique(symmrel_kpoints,'rows');

    unitcell_kpoints=[unitcell_kpoints; symmrel_kpoints];
%   following commented out as weights not given correctly in expanded lattice klist
%    if size(symmrel_kpoints,1)~=IBZ_klist(IBZ_kpointnum,6)
%        %no. of distinct equiv points is not as given in weight column in
%        %file
%        disp('No. of symm rel points not correct');
%    end
end
%generate list of bandnums
if spinpol
    Wien2k_bandnums=sort([bandnumlist; bandnumlist]);
    spindirns=repmat([1 2],1,length(bandnumlist));
else
    Wien2k_bandnums=bandnumlist;
    spindirns=ones(size(bandnumlist));
end

%unitcell_kpoints is 4 col list of all kpoints inside first unit cell
%cols: IBZ_kpoint_num (for energy lookup), n_a, n_b, n_c
%where n_a is multiple of bcc spanning vector i.e. n*a/dims(1)

if isexpanded
    disp('Populating k-mesh with E in selected bands');
    %expanded sc klist currently contains equivalent points so remove these
    [unitcellcoords index]=unique(unitcell_kpoints(:,2:4),'rows');
    unitcell_kpoints=[unitcell_kpoints(index,1) unitcellcoords];
    %generate Cartesian coords in bounding cuboid  
    switch rec_latt_type
        %exp_klist program inconsistent as dims in .klist file are new for 'hex' and old for 'bcc'
        case 'bcc'
            [cartX cartY cartZ]=ndgrid(0:dims(1)*2-1,0:dims(2)*2-1,0:dims(3)*2-1);
        case 'hex'
            [cartX cartY cartZ]=ndgrid(0:dims(1)+dims(2),0:dims(2),0:dims(3));
    end
    cuboidcoords=[reshape(cartX,[],1) reshape(cartY,[],1) reshape(cartZ,[],1)]*inv(n_rec_latt_vecs);
    cuboidcoords=[mod(cuboidcoords(:,1),dims(1)) mod(cuboidcoords(:,2),dims(2)) mod(cuboidcoords(:,3),dims(3))];
    cuboidunitcellserial=cuboidcoords(:,3)+dims(3)*cuboidcoords(:,2)+dims(3)*dims(2)*cuboidcoords(:,1);
    %define cuboidindex which maps to reshaped cartX,Y,Z and contains
    %indexes into unitcell_kpoints
    [dummy dummy cuboidindex]=unique(cuboidcoords,'rows');

    for relbandnum=1:length(Wien2k_bandnums)
        cartE{relbandnum}=energy(Wien2k_bandnums(relbandnum),spindirns(relbandnum),unitcell_kpoints(cuboidindex,1));
        cartE{relbandnum}=reshape(cartE{relbandnum},[2*dims(3) 2*dims(2) 2*dims(1)]);        
    end
    switch rec_latt_type
        case 'bcc'
            [cartX cartY cartZ]=ndgrid(0:dims(1)*2,0:dims(2)*2,0:dims(3)*2);
            cartX=cartX*abs(spanning_vecs(1,1));
            cartY=cartY*abs(spanning_vecs(1,1));
            cartZ=cartZ*abs(spanning_vecs(1,1));
        case 'hex'
            [cartX cartY cartZ]=ndgrid(0:dims(1)*2,0:dims(2),0:dims(3));
            cartX=cartX*spanning_vecs(1,1)/2;
            cartY=cartY*spanning_vecs(2,2);
            cartZ=cartZ*spanning_vecs(3,3);
    end
    cartE=addborders(cartE);
    cuboidextents=[min(reshape(cartX,[],1)) min(reshape(cartY,[],1)) min(reshape(cartZ,[],1));...
        max(reshape(cartX,[],1)) max(reshape(cartY,[],1)) max(reshape(cartZ,[],1))];
    dx=abs(cartX(2,1,1)-cartX(1,1,1));
    dy=abs(cartY(1,2,1)-cartY(1,1,1));
    dz=abs(cartZ(1,1,2)-cartZ(1,1,1));
end

%sort into raster order (i.e. bxsf order (a1,b1,c1), (a1,b1,c2) etc. in
%case of rec_latt_based array)
%maxs=max(unitcell_kpoints);
%sortterm=unitcell_kpoints(:,4)+unitcell_kpoints(:,3)*maxs(4)+unitcell_kpoints(:,2)*maxs(4)*maxs(3);
%[dummy order]=sort(sortterm);
%unitcell_kpoints(:,:)=unitcell_kpoints(order,:);

%disp('Populating k-mesh with E in selected bands');

%for relbandnum=1:length(Wien2k_bandnums)
%    energyBZ{relbandnum}=squeeze(energy(Wien2k_bandnums(relbandnum),spindirns(relbandnum),unitcell_kpoints(:,1)));
%    energyBZ{relbandnum}=reshape(energyBZ{relbandnum},flipdim(dims,1)');
%    energyBZ{relbandnum}=addborders(energyBZ{relbandnum});
%    if isexpanded
%        cartE{relbandnum}=energyBZ{relbandnum};
%        dx=rec_latt_vecs(1,1)/(2*dims(1));
%        dy=rec_latt_vecs(2,2)/(2*dims(2));
%        dz=rec_latt_vecs(3,3)/(2*dims(3));
%        cuboidextents=
%    else
%        [cartE(relbandnum) dx dy dz cuboidextents]=interp2Cart(energyBZ{relbandnum},rec_latt_vecs,rec_latt_type);
%    end
%    [DOS_E{relbandnum} DOS{relbandnum} IDOS{relbandnum}]=calcDOS(0.001,cartE{relbandnum},1);
%    disp(['Band ' num2str(relbandnum) '/' num2str(length(Wien2k_bandnums)) ' done.']);      
%end    

%[cartX cartY cartZ]=ndgrid(cuboidextents(1,1):dx:cuboidextents(2,1),cuboidextents(1,2):dy:cuboidextents(2,2),...
%    cuboidextents(1,3):dz:cuboidextents(2,3));

disp('Calculating DOS, IDOS in selected bands');
for relbandnum=1:length(Wien2k_bandnums)
    [DOS_E{relbandnum} DOS{relbandnum} IDOS{relbandnum}]=calcDOS(0.001,cartE{relbandnum},1);
    disp(['Band ' num2str(relbandnum) '/' num2str(length(Wien2k_bandnums)) ' done.']);      
end

cartdims=size(cartX);

BZfacenormals=getBZfacenormals(rec_latt_vecs,rec_latt_type);
searchvolnormals=[1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];

save([outfilename '.bands.mat'],'spanning_vecs','rec_latt_vecs','dims','FermiLevel','cartX','cartY','cartZ',...
    'Wien2k_bandnums','spindirns','cartE','cuboidextents','cartdims','dx','dy','dz','DOS_E','DOS','IDOS','BZfacenormals','searchvolnormals','pathcasename');
a=1;

function result=addborders(matrixorlist);
if ~iscell(matrixorlist)
    matrixlist{1}=matrixorlist;
else
    matrixlist=matrixorlist;
end
for mnum=1:length(matrixlist)
    result{mnum}=cat(1,matrixlist{mnum},matrixlist{mnum}(1,:,:));
    result{mnum}=cat(2,result{mnum},result{mnum}(:,1,:));
    result{mnum}=cat(3,result{mnum},result{mnum}(:,:,1));
end
if ~iscell(matrixorlist)
    nresult=result{1};
    result=nresult;
end
    
