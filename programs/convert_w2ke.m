function convert_W2kE(pathcasename,bandnumlist,ForceFermiLevel,interpmethod,outfilename,varargin);
%function convert_W2kE(pathcasename,bandnumlist,ForceFermiLevel,interpmethod,outfilename);
%
%reads W2k energy band files (.energy) output by lapw1 and converts to Matlab structure
%filename args without extensions
%bandnumlist contains bandnums to be converted. if empty only those crossing E_F are converted
%ForceFermiLevel: value to force EF (otherwise if empty, EF is read from Wien2k output file case.SCF2)
%interpmethod: 'linear' or 'nearest'
%requires files: case.klist, case.outputkgen, case.struct, case.SCF2(up/dn), case.energy(up/dn),
%case.output2(up/dn)

b.pathcasename=pathcasename;
readqtl=false;
plotbandranges=false;
load('constants');

pchgs=[]; pchglbls=[];

spinorbstr='';

if ~isempty(varargin)
    fexpfac=varargin{1};
end

if exist([b.pathcasename '.energyso'])==2
    spinorbinp=input('Spin-orbit energy file (case.energyso) is present. Use this file? (Y/N) ','s');
    if length(spinorbinp)>0 & strcmp(upper(spinorbinp(1)),'Y')
        spinorbstr='so';
    end
end
if exist([b.pathcasename '.energy' spinorbstr 'up'])==2 & exist([b.pathcasename '.energy' spinorbstr 'dn'])==2
    spinpol=true;
    dirnstr={'up' 'dn'};
    dirntextstr={'(spin up) ' '(spin dn) '};
else
    spinpol=false;
    dirnstr={''};
    dirntextstr={''};
end
if (spinpol==false & exist([b.pathcasename '.qtl'])==2) | (spinpol==true & exist([b.pathcasename '.qtlup'])==2 & exist([b.pathcasename '.qtlup'])==2)
     readqtl=true;
end
readqtl=false;%%%%%%%%%%%%%
if isempty(bandnumlist)
    userbandlist=false;
    if spinpol
        bandnumlist={[] []};
    else
        bandnumlist={[]};
    end
else
    userbandlist=true;
end

%get info from .struct file
infile=fopen([b.pathcasename '.struct']);
for lnum=1:4
    linestr{lnum}=fgetl(infile);
end
fclose(infile);
b.latt_type=linestr{2}(1);
b.latt_params=Bohr_radius*[str2num(linestr{4}(1:10)) str2num(linestr{4}(11:20)) str2num(linestr{4}(21:30))];
b.latt_angles=[str2num(linestr{4}(31:40)) str2num(linestr{4}(41:50)) str2num(linestr{4}(51:60))];

%read symm operations
[b.symmops b.rlvs]=getsymmops([b.pathcasename '.outputkgen']);
numsymmops=size(b.symmops,3);

%define matrix of rec latt vecs such that x-components are normalized appropriately
denom=max(abs(b.rlvs),[],1);
b.n_rlvs=b.rlvs./repmat(denom,3,1); %repmat(diag(b.rlvs)',3,1);

%.klist contains: dims; IBZ kpoint no., rec latt coord kx/(a_x/d),ky/(b_y/d),kz/(c_z/d), d=denominator,
%weight (i.e. no. of symmetry rel points incl IBZ one)
disp('Reading dims, Cartesian IBZ coords from .klist');
infile=fopen([b.pathcasename '.klist'],'r');
if infile==-1
    disp('File not found: case.klist');
    return;
end

%get dims from last columns of first line of .klist
firstlinestr=fgetl(infile);
startchar=findstr(firstlinestr,'(')+1;
endchar=findstr(firstlinestr,')')-1;
b.W2kdims=sscanf(firstlinestr(startchar:endchar),'%f',3);
fclose(infile);
b.spanning_vecs=b.rlvs./(repmat(b.W2kdims,1,3));

%get IBZ kpoint coords
infile=fopen([b.pathcasename '.klist'],'r');
notatend=true; linenum=1;
while notatend
    linestr=fgetl(infile);
    if strcmp(linestr(1:3),'END')
        notatend=false;
        break;
    end
    b.IBZ_klist(linenum,:)=sscanf(linestr,'%f',6);
    %use equivpts routine to check weight makes sense
    if b.latt_type=='H'
        [eqpts weight]=equivkpts(b.IBZ_klist(linenum,2:4),b.symmops,b.IBZ_klist(linenum,5));
    else
        [eqpts weight]=equivkpts(b.IBZ_klist(linenum,2:4)*inv(b.n_rlvs),b.symmops,b.IBZ_klist(linenum,5));
    end
    if weight~=b.IBZ_klist(linenum,6)
        disp(['IBZ pt no. ' num2str(linenum) ': weight in .klist does not agree with weight calc by equivpts']);
    end
    linenum=linenum+1;
end
IBZ_numkpoints=linenum-1;
b.LCMdivs=unique(b.IBZ_klist(:,5));
if numel(b.LCMdivs)~=1
    disp('Division unit changes within .klist?');
    return;
end
disp(['There are ' num2str(IBZ_numkpoints) ' points in .klist file']);

%get Fermi energy of final SCF cycle, bandranges from case.output2
for dirnnum=1:length(dirnstr)
    gotEFfromSCF2=false;
    if ~isempty(ForceFermiLevel)
        b.FermiLevel(dirnnum)=ForceFermiLevel(dirnnum);
        disp(['Forcing E_F(' dirnstr{dirnnum} ')=' num2str(b.FermiLevel(dirnnum))]);
    else
        infile=fopen([b.pathcasename '.scf2' dirnstr{dirnnum}]);
        notatend=true; linestr=[]; 
        while notatend;
            prevlinestr=linestr;
            linestr=fgetl(infile);    
            if findstr(linestr,':FER')
                gotEFfromSCF2=true;
                notatend=false;
                break;
            elseif linestr==-1
                notatend=false;
                b.FermiLevel=input('case.scf2 does not contain '':FER'' -> Enter E_F or [E_F_up E_F_dn] (Ryd): ')
            end
        end
        fclose(infile);
    end
    if gotEFfromSCF2
        b.FermiLevel(dirnnum)=str2num(linestr(findstr(linestr,'=')+1:end));
        disp(['From case.scf2: E_F(' dirnstr{dirnnum} ')=' num2str(b.FermiLevel(dirnnum))]);
    end
    bandranges{dirnnum}=getbandranges(b.pathcasename, dirnstr{dirnnum});
    if ~isempty(bandranges{dirnnum})
        nobandranges=false;
        numbands(dirnnum)=size(bandranges{dirnnum},1);
        if dirnnum==2 & numbands(2)>numbands(1)
            FermiCrossing=[FermiCrossing; zeros(numbands(2)-numbands(1),3)];
        elseif dirnnum==2 & numbands(2)<numbands(1)
            bandranges{2}=[bandranges{2}; zeros(numbands(1)-numbands(2),3)];
        end
        FermiCrossing(:,dirnnum)=bandranges{dirnnum}(:,1)<b.FermiLevel(dirnnum) & bandranges{dirnnum}(:,2)>b.FermiLevel(dirnnum);
        disp([dirntextstr{dirnnum} 'No. of bands: ' num2str(numbands(dirnnum)) '. Bands ' num2str(find(FermiCrossing(:,dirnnum)')) ' cross Fermi level.']);
        if isempty(bandnumlist{dirnnum})
           %by default only keep bands that cross the Fermi surface
           bandnumlist{dirnnum}=find(FermiCrossing(:,dirnnum));
           disp(['Band ranges read from case.output2. E_F crossing bands (' dirnstr{dirnnum} '): ' num2str(bandnumlist{dirnnum}')]);
        end
    else
        nobandranges=true;
    end
end
%read energies from case.energy files
for dirnnum=1:length(dirnstr)    
    if nobandranges
        dispstr=['Reading IBZ energies for ALL bands'];
    else
        dispstr=['Reading IBZ energies for bands (' dirnstr{dirnnum} '): ' num2str(bandnumlist{dirnnum}')];
    end
    b.energyfname{dirnnum}=[b.pathcasename '.energy' spinorbstr dirnstr{dirnnum}];
    disp([dispstr ' from ' b.energyfname{dirnnum}]);
    dirnenergy=importenergy(b.energyfname{dirnnum},bandnumlist{dirnnum},IBZ_numkpoints);
    numbands(dirnnum)=size(dirnenergy,1);   
    if nobandranges
        bandnumkey{dirnnum}=1:numbands(dirnnum);
    else
        bandnumkey{dirnnum}=bandnumlist{dirnnum};
    end
    if isempty(dirnenergy)
        return;
    else
        b.IBZenergy{dirnnum}=dirnenergy;
    end
end
for dirnnum=1:length(dirnstr)
    b.IBZenergy{dirnnum}(b.IBZenergy{dirnnum}==0)=nan;
end
clear dirnenergy;

%read partial charges from case.qtl
if readqtl
    for dirnnum=1:length(dirnstr)
        qtlfname=[b.pathcasename '.qtl' spinorbstr dirnstr{dirnnum}];
        [dpchgs b.pchglbls]=importQTL(qtlfname,bandnumlist{dirnnum},IBZ_numkpoints);
        b.pchgs{dirnnum}=dpchgs;
    end
end
    
    
if nobandranges & ~userbandlist %if case.output2 did not contain band ranges, we need to import all bands and find which cross E_F
    for dirnnum=1:length(dirnstr)
        bandranges{dirnnum}=[min(b.IBZenergy{dirnnum},[],2) max(b.IBZenergy{dirnnum}(:,:),[],2)];
        FermiCrossing(:,dirnnum)=bandranges{dirnnum}(:,1)<b.FermiLevel(dirnnum) & bandranges{dirnnum}(:,2)>b.FermiLevel(dirnnum);
        bandnumlist{dirnnum}=find(FermiCrossing(:,dirnnum));
        disp([dirntextstr{dirnnum} 'No. of bands: ' num2str(numbands(dirnnum)) '. Bands ' num2str(find(FermiCrossing(:,dirnnum)')) ' cross Fermi level.']);
    end
end

%plot band ranges
if plotbandranges
    symb={'+','x'};
    fig=findobj('Tag','bandranges');
    if ishandle(fig)
        figure(fig);
    else
        figure('Tag','bandranges');
    end
    hold on;
    for dirnnum=1:length(dirnstr)
        for bnum=1:size(bandranges{dirnnum},1)
            plot(bnum,bandranges{dirnnum}(bnum,1),symb{dirnnum},bnum,bandranges{dirnnum}(bnum,2),symb{dirnnum});
        end
    end
    line([0; size(bandranges{dirnnum},1)],[b.FermiLevel(dirnnum); b.FermiLevel(dirnnum)]);
end

%remove bands whose energies are not defined on all IBZ kpoints
for dirnnum=1:length(dirnstr)
    for bnum=numel(bandnumlist{dirnnum}):-1:1
        if numel(find(isnan(reshape(b.IBZenergy{dirnnum}(find(bandnumkey{dirnnum}==bandnumlist{dirnnum}(bnum)),:),[],1))))==0
            fulldef(bnum)=true;
        else
            fulldef(bnum)=false;
        end
    end
    bandnumlist{dirnnum}=bandnumlist{dirnnum}(fulldef);
    bandranges{dirnnum}=bandranges{dirnnum}(fulldef);
end

%populate entire BZ from IBZ data
%generate full unit. cell matrix
%convert inequiv kpoints to rec_latt_basis and then operate with symmetry
%matrices to give equiv points

%define matrix n_rlvs such that:
%(reclatt based coord)*matrix=(normalized cartesian coord)
b.rlv_mags=max(abs(b.rlvs),[],1);;
n_rlvs=b.rlvs./(repmat(b.rlv_mags,[3 1]));
[int_rlvs denom]=rat(n_rlvs,1e-5); %avoids rounding error
n_rlvs=int_rlvs./denom;

disp('Generating k-mesh over whole unit cell');
b.uc_klist=IBZ2unitcell(b.IBZ_klist,n_rlvs,b.rlv_mags,[],b.latt_type,b.LCMdivs,b.symmops,b.W2kdims);

if spinpol    
    b.Wien2k_bandnums=[bandnumlist{1}; bandnumlist{2}];
    b.spindirns=[ones(numel(bandnumlist{1}),1); 2*ones(numel(bandnumlist{2}),1)];
else
    b.Wien2k_bandnums=bandnumlist{1};
    b.spindirns=ones(numel(bandnumlist{1}),1);
end

for dirnnum=1:length(dirnstr)
    %generate list of bandnums
    if ~isempty(bandnumlist{dirnnum})
        disp(['Converting bands (' dirnstr{dirnnum} '): ' num2str(bandnumlist{dirnnum}')]);
    end
end
if (spinpol & isempty(bandnumlist{1}) & isempty(bandnumlist{2})) | (~spinpol & isempty(bandnumlist{1}))
    disp('No bands to convert!');
    return;
end    

disp('Calculating E on orthomesh over tileable cuboid.');
switch b.latt_type
    case {'F','B'}
        b.mcell=[2 2 2];
        % b.cuboidextents=[0 0 0; b.mcell.*abs(diag(b.rlvs)')];
        t=b.rlvs; t(t==0)=nan; b.cuboidextents=[0 0 0; b.mcell.*min(t,[],1)];
        %b.cuboidextents=[0 0 0; b.mcell.*min(abs(b.rlvs))];
    case {'C'}
        b.mcell=[2 2 1];
        b.cuboidextents=[0 0 0; b.mcell.*abs(diag(b.rlvs)')];
    case 'P'
        b.mcell=[1 1 1];
        b.cuboidextents=[0 0 0; b.mcell.*abs(diag(b.rlvs)')];
    case 'H'        
        b.mcell=[1 2 1];
        b.cuboidextents=[0 0 0; b.mcell.*abs(diag(b.rlvs)')];
end
%make energy matrix with bands in columns, all up then all dn
energy2=[];
for dirnnum=1:length(dirnstr)
    for bnum=1:numel(bandnumlist{dirnnum})
        energy2(end+1,:)=squeeze(b.IBZenergy{dirnnum}(find(bandnumkey{dirnnum}==bandnumlist{dirnnum}(bnum)),:));
    end
end
disp(['Using submesh shift [' num2str(b.IBZ_klist(1,2:4)./(b.LCMdivs./b.W2kdims')) '] (rec latt spanning vecs)']);
b.cart_submesh_shift=b.IBZ_klist(1,2:4)./(b.LCMdivs./b.W2kdims')*b.n_rlvs;

if strcmp(lower(interpmethod),'fourier')
    newdims=(b.mcell.*b.W2kdims'+1-1)*fexpfac+1;
    for dirnnum=1:numel(dirnstr)
        if dirnnum==1
            [b.cartE(b.spindirns==dirnnum) b.cartX b.cartY b.cartZ]=fexp(b,b.pathcasename,b.energyfname{dirnnum},newdims,b.mcell);
        else
            b.cartE(b.spindirns==dirnnum)=fexp(b,b.pathcasename,b.energyfname{dirnnum},newdims,b.mcell);
        end
    end
    b.cartdims=newdims;
    disp(['Fourier interp [' num2str(size(b.cartX)) '] to [' num2str(newdims) '].']);
else
    [b.cartX b.cartY b.cartZ b.cartE]=orthomeshe_griddata(energy2,b.uc_klist,b.W2kdims,b.LCMdivs,b.latt_type,b.n_rlvs,b.rlv_mags,b.cart_submesh_shift,b.mcell);
end

b.dx=abs(b.cartX(2,1,1)-b.cartX(1,1,1));
b.dy=abs(b.cartY(1,2,1)-b.cartY(1,1,1));
b.dz=abs(b.cartZ(1,1,2)-b.cartZ(1,1,1));
b.cuboidextents=[min([b.cartX(1:end)' b.cartY(1:end)' b.cartZ(1:end)']); ...
    max([b.cartX(1:end)' b.cartY(1:end)' b.cartZ(1:end)'])];
b.cartdims=size(b.cartX);

b.BZfacenormals=getBZfacenormals(b.latt_type);
b.searchvol=2*b.BZfacenormals*b.rlvs;
b.centresvol=1.1*b.BZfacenormals*b.rlvs;
b.ptinvol=[0 0 0];



save([outfilename '.bands.mat'], '-struct','b');
disp(['Done. Output file ' outfilename '.bands.mat']);
