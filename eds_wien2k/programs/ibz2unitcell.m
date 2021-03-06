function uc_klist=IBZ2unitcell(IBZklist,n_rec_latt_vecs,rl_vecmags,latt_params,latt_type,LCMdivs,symmops,dims);
%function uc_klist=IBZ2unitcell(IBZklist,n_rec_latt_vecs,rl_vecmags,latt_params,latt_type,LCMdivs,symmops);
%
%returns an (n x 4) matrix containing (IBZ serial no., na, nb, nc) of every point generated from IBZ
%klist that lies within the unit cell
%requires IBZklist with columns: serial, 'kx,ky,kz', LCMdivs, weight
%IBZklist must contain all points on an orthorhombic mesh
%rec_latt_vecs in Angst^-1
%LCMdivs is denominator of klist coords from .klist
%symmops contains symm operations required for IBZ->unit cell

%15/08/06 transpose needed in line 34 for H but not others???
%11/12/06 made some change and now needed for 'B'????
%05/03/06 changed line 21, num_extra_kpts=3 to =0
inv_n_rec_latt_vecs=inv(n_rec_latt_vecs);

switch latt_type
    case {'F','C'}     
        rl_IBZklist=IBZklist(:,2:4)*inv_n_rec_latt_vecs;
        num_extra_kpts=0;
    case {'B','P'}
        rl_IBZklist=IBZklist(:,2:4)*inv_n_rec_latt_vecs;
        num_extra_kpts=0;
    case 'H'
        %.klist coords are cartesian for most lattice types but not hex
        rl_IBZklist=IBZklist(:,2:4);
        num_extra_kpts=0;
end

uc_klist=[]; tetrah_pts=0;
for IBZkpn=1:size(rl_IBZklist,1)
    weight=IBZklist(IBZkpn,6);
    symmrel_kpts=[];
    if strcmp(latt_type,'H')
        for symmop=1:size(symmops,3)
            symmrel_kpts(symmop,:)=rl_IBZklist(IBZkpn,:)*symmops(:,:,symmop)'; %NB see note at top about transpose
        end
    else
        for symmop=1:size(symmops,3)
            symmrel_kpts(symmop,:)=rl_IBZklist(IBZkpn,:)*symmops(:,:,symmop)'; %NB see note at top about transpose
        end
    end
    symmrel_kpts=mod(symmrel_kpts,LCMdivs);       
    symmrel_kpts=unique(symmrel_kpts,'rows');
    extra_kpt_num=mod(IBZkpn-1,num_extra_kpts+1);
    if extra_kpt_num==0 & size(symmrel_kpts,1)~=weight
        disp(['IBZ k-point no. ' num2str(IBZkpn) ': Weight in .klist wrong?']);
    end
    if extra_kpt_num==0
        tetrah_pts=tetrah_pts+size(symmrel_kpts,1);
    end
    uc_klist=[uc_klist; [repmat(IBZkpn,[size(symmrel_kpts,1) 1]) symmrel_kpts]];    
end
%remove equivalent points
[uuc_klist uc_indices]=unique(uc_klist(:,2:4),'rows');
if size(uuc_klist,1)~=size(unique(uc_klist,'rows'),1)
    disp('>1 IBZ point for some point(s) in unit cell - either .klist contains some equivalent points or something is wrong.');
else
    disp([num2str(tetrah_pts) ' unit cell kpts generated from IBZ kpts. ' num2str(dims(1)) '*' num2str(dims(2)) '*' num2str(dims(3)) '=' num2str(prod(dims))]);
end
uc_klist=[uc_klist(uc_indices,1) uuc_klist];

