function [cartX cartY cartZ cartE]=orthomeshe_griddata(IBZdata,uc_klist,dims,LCMdivs,latt_type,n_rec_latt_vecs,rl_vecmags,cart_submesh_shift,mcell);

%interpolate data given on IBZ (energy or partial charge) onto Cartesian mesh
%of size mcell.*[a*_x b*_y c*_z] which is smallest tileable cuboid

%IBZdata is numbands*num_IBZ_kpts array containing energies for bands in rows (all spin up then
%all spin dn)
%uc_klist is n*4 array containing all kpts in full zone with cols [IBZ_ptno. k_a* k_b* k_c*]
%dims is no. of kpt divisions along a*, b*, c* as given in case.klist
%LCMdivs is lowest common multiple of the intervals for a*, b*, c* in case.klist units
%latt_type - see below
%n_rec_latt_vecs is a 3*3 array containing normalized rec latt vecs
%rl_vecmags contains the magnitudes in Angst^-1 of the rlv's
%cart_submesh_shift is the submesh shift used by W2k

%defines orthorhombic mesh spacing in terms of spanning vectors a_x/dims(1),b_y/dims(2),c_z/dims(3)

switch latt_type
    case {'F','B'}
        mcell=[2 2 2];
        dx=1; dy=1; dz=1;
    case {'C'}
        mcell=[2 2 1];
        dx=1; dy=1; dz=1;
    case 'P'
        mcell=[1 1 1];
        dx=1; dy=1; dz=1;
    case 'H'
        mcell=[1 2 1];
        dx=0.5; dy=1; dz=1;
end

%first generate cartesian coords of orthorhombic kpoints over matrix cell,
%in units of spanning vector components
%this is 'destination' grid for interp
[cartX cartY cartZ]=ndgrid(cart_submesh_shift(1)+(0:dx:mcell(1)*dims(1)),cart_submesh_shift(2)+(0:dy:mcell(2)*dims(2)),cart_submesh_shift(3)+(0:dz:mcell(3)*dims(3)));

%now express the points in the destination grid in terms of W2k coords wrapped back to the unit cell
%use fact that every W2k mesh point coincides with an orthomesh point
%convert cartXYZ to wrapped uc_klist coords and make index column
mc_W2k_klist=[cartX(1:end)' cartY(1:end)' cartZ(1:end)']*inv(n_rec_latt_vecs);
mc_W2k_klist=mc_W2k_klist;
mc_W2k_klist=[mod(mc_W2k_klist(:,1),dims(1))*LCMdivs/dims(1) ...
    mod(mc_W2k_klist(:,2),dims(2))*LCMdivs/dims(2) ...
    mod(mc_W2k_klist(:,3),dims(3))*LCMdivs/dims(3)];
[ismem index]=ismember(mc_W2k_klist,uc_klist(:,2:4),'rows');
%ismem has same number elements as dest grid and is true at each point that is equiv to a W2k point
%index contains the index into uc_klist for each point that is equiv to a W2k point
disp(['Interping to orthomesh. Pts in matrix cell=' num2str(size(mc_W2k_klist,1)) ', ' num2str(numel(find(ismem))) ' coincide with W2k mesh pt']);
%relate back to Cart coords
cart_mc_W2k_klist=[uc_klist(index(ismem),1) cartX(ismem) cartY(ismem) cartZ(ismem)];
%this is source grid for interp

%interpolate Wien2k IBZ data to matrix cell
cartX=permute(cartX,[2 1 3]); cartY=permute(cartY,[2 1 3]); cartZ=permute(cartZ,[2 1 3]);
for bandnum=1:size(IBZdata,1)
    cartE{bandnum} = griddata3(cart_mc_W2k_klist(:,2),cart_mc_W2k_klist(:,3),cart_mc_W2k_klist(:,4),...
        IBZdata(bandnum,cart_mc_W2k_klist(:,1)),cartX,cartY,cartZ,'linear');
    cartE{bandnum}=reshape(cartE{bandnum},size(cartX));
    cartE{bandnum}=permute(cartE{bandnum},[2 1 3]); 
end
cartX=permute(cartX,[2 1 3]); cartY=permute(cartY,[2 1 3]); cartZ=permute(cartZ,[2 1 3]);

%scale cart coords back to Anstrom^-1
cartX=cartX/max(cartX(1:end))*rl_vecmags(1)*mcell(1);
cartY=cartY/max(cartY(1:end))*rl_vecmags(2)*mcell(2);
cartZ=cartZ/max(cartZ(1:end))*rl_vecmags(3)*mcell(3);
