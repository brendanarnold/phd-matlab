function [cartX cartY cartZ indices uc_indices]=IBZ2mat2(IBZklist,rec_latt_vecs,latt_params,latt_type,LCMdivs,symmops);
%function [cartX cartY cartZ indices]=IBZ2mat(IBZklist,rec_latt_vecs);
%
%returns cartXYZ, coordinate matrices (Angst^-1) and indices, a matrix of IBZklist serial
%for each position in the new energy matrix
%requires IBZklist with columns: serial, 'kx,ky,kz', LCMdivs, weight
%IBZklist must contain all points on an orthorhombic mesh
%rec_latt_vecs in Angst^-1
%LCMdivs is denominator of klist coords from .klist
%symmops contains symm operations required for IBZ->unit cell

%define matrix n_rec_latt_vecs such that:
%(reclatt based coord)*matrix=(normalized cartesian coord)
n_rec_latt_vecs=rec_latt_vecs./(repmat(diag(rec_latt_vecs)',[3 1]));
[int_rec_latt_vecs denom]=rat(n_rec_latt_vecs,1e-5); %avoids rounding error
n_rec_latt_vecs=int_rec_latt_vecs./denom;
inv_n_rec_latt_vecs=inv(n_rec_latt_vecs);

switch latt_type
    case {'F','P','C','B'}
        LCMmult=[2 2 2];
        extra_kpts=[0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0];
        rl_IBZklist=IBZklist(:,2:4)*inv_n_rec_latt_vecs;
        %.klist coords are cartesian for these lattices but not for others
        ortho_mesh_cell=[1 1 1];
    case 'H'
        LCMmult=[1 2 1];
        extra_kpts=[];
        rl_IBZklist=IBZklist(:,2:4);
        %specify orthorh. mesh uc lengths in terms of tetrah. mesh lengths ax,by,cz
        ortho_mesh_cell=[1 2 1];               
end
ortho_meshcell_vol=prod(ortho_mesh_cell);
tetra_meshcell_vol=dot(n_rec_latt_vecs(1,:),cross(n_rec_latt_vecs(2,:),n_rec_latt_vecs(3,:)));
num_extra_kpts=size(extra_kpts,1);

%get size of tetrahedral mesh unit cell
tet_rl_IBZklist=rl_IBZklist(mod(0:size(rl_IBZklist)-1,num_extra_kpts+1)==0,:);
ua=unique(tet_rl_IBZklist(:,1)); ub=unique(tet_rl_IBZklist(:,2)); uc=unique(tet_rl_IBZklist(:,3));
if ~isequal(ua',linspace(ua(1),ua(end),numel(ua))) | ~isequal(ub',linspace(ub(1),ub(end),numel(ub))) | ...
        ~isequal(uc',linspace(uc(1),uc(end),numel(uc)))
    disp('Wien2k klist points are not on tetrahedral mesh.');
else
   meshunit=[ua(2)-ua(1) ub(2)-ub(1) uc(2)-uc(1)];   
end
ntetra_meshcells=prod(LCMdivs./meshunit);

uc_klist=[]; tetrah_pts=0;
for IBZkpn=1:size(rl_IBZklist,1)
    weight=IBZklist(IBZkpn,6);
    symmrel_kpts=[];
    for symmop=1:size(symmops,3)
        symmrel_kpts(symmop,:)=(symmops(:,:,symmop)*rl_IBZklist(IBZkpn,:)')';
    end
    symmrel_kpts=mod(symmrel_kpts,LCMdivs);       
    symmrel_kpts=unique(symmrel_kpts,'rows');
    extra_kpt_num=mod(IBZkpn-1,num_extra_kpts+1);
    %if extra_kpt_num~=0
        %this kpoint is one added to klist make an orthorhombic mesh
        %only keep those points that are in correct place in mesh 'unit cell'
    %    for pnum=size(symmrel_kpts,1):-1:1
    %        if ~any(all(repmat(mod(symmrel_kpts(pnum,:)./meshunit,1),size(extra_kpts,1),1)==extra_kpts,2));
                %symmrel_kpts(pnum,:)=[];
    %        end
    %    end
    %else
        if size(symmrel_kpts,1)~=weight
            disp('Weight in .klist wrong?');
        end
        tetrah_pts=tetrah_pts+size(symmrel_kpts,1);
    %end
    uc_klist=[uc_klist; [repmat(IBZkpn,[size(symmrel_kpts,1) 1]) symmrel_kpts]];    
end
%remove equivalent points
[uuc_klist uc_indices]=unique(uc_klist(:,2:4),'rows');
if size(uuc_klist,1)~=size(unique(uc_klist,'rows'),1)
    disp('Something wrong. >1 IBZ point for some point(s) in unit cell');
end
uc_klist=[uc_klist(uc_indices,1) uuc_klist];

%define uc_indices to have same size as unit cell tetra-mesh and contain index to IBZ point
uc_indices=reshape(uc_klist(:,1),[LCMdivs./meshunit([3 2 1])]);
uc_indices=permute(uc_indices,[3 2 1]);

%remove points that do not lie on orthorhombic mesh specified by ortho_mesh_cell
ortho_meshuc_pts=mod(uc_klist(:,2:4)*n_rec_latt_vecs,repmat(meshunit.*ortho_mesh_cell,size(uc_klist,1),1));
on_orthomesh=all(repmat([0 0 0],size(ortho_meshuc_pts,1),1)==ortho_meshuc_pts,2);
numon_orthomesh=numel(find(on_orthomesh));
if numon_orthomesh~=ntetra_meshcells*tetra_meshcell_vol/ortho_meshcell_vol
    disp(['No. of kpoints on ortho-mesh is wrong (' num2str(numon_orthomesh) ')']);
    return;
end
disp([num2str(numon_orthomesh) ' out of ' num2str(size(uc_klist,1)) ' kpoints in unit cell are on ortho-mesh - discarding remainder']);
ortho_uc_klist=uc_klist(on_orthomesh,:);

%transform the unit cell reclatt list to cart coords and check it is on orthorhombic mesh
cart_uc_klist=ortho_uc_klist(:,2:4)*n_rec_latt_vecs;
uc_xs=unique(cart_uc_klist(:,1)); uc_ys=unique(cart_uc_klist(:,2)); uc_zs=unique(cart_uc_klist(:,3));
if ~isequal(uc_xs,linspace(uc_xs(1),uc_xs(end),length(uc_xs))') | ...
        ~isequal(uc_ys,linspace(uc_ys(1),uc_ys(end),length(uc_ys))') | ~isequal(uc_zs,linspace(uc_zs(1),uc_zs(end),length(uc_zs))')
    disp('.klist points not on orthorhombic mesh');
end
dx=uc_xs(2)-uc_xs(1); dy=uc_ys(2)-uc_ys(1); dz=uc_zs(2)-uc_zs(1);

%generate cart grids in normed units
[cartX cartY cartZ]=ndgrid(0:dx:dx*(LCMdivs/meshunit(1)-1),...
    0:dy:dy*(LCMdivs/meshunit(2)-1),...
    0:dz:dz*(LCMdivs/meshunit(3)-1));

%make matrix of IBZ indices over cuboid
%convert ortho X,Y,Z grid in cuboid to rec latt coords
cub_rl=[cartX(1:end)' cartY(1:end)' cartZ(1:end)']*inv_n_rec_latt_vecs;

%cub_rl=cub_rl.*repmat([2 1 1],size(cub_rl,1),1);
cub_rl=mod(cub_rl,LCMdivs);

indA=1+cub_rl(:,1)/meshunit(1);
indB=1+cub_rl(:,2)/meshunit(2);
indC=1+cub_rl(:,3)/meshunit(3);

indices=uc_indices(sub2ind(size(uc_indices),indA(1:end)',indB(1:end)',indC(1:end)'));

%for n=1:size(cub_rl,1)
%    matchabc=mod(cub_rl(n,:),LCMdivs);
%    matchnum=find(all(repmat(matchabc,size(uc_klist,1),1)==uc_klist(:,2:4),2));
%    indices(n)=uc_klist(matchnum,1);
%end

abscartunits=(diag(rec_latt_vecs)')/LCMdivs;
cartX=cartX*abscartunits(1); cartY=cartY*abscartunits(2); cartZ=cartZ*abscartunits(3);
