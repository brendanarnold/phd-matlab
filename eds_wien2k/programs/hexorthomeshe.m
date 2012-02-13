function [cartX cartY cartZ cartE]=hexorthomeshe(energy,uc_klist,dims,LCMdivs,latt_type,n_rec_latt_vecs,rl_vecmags,interpmethod);

%for hex crystal, interpolate Wien2k energy data from hexagonal mesh within unit cell
%'uc_klist' contains k coords of hexagonal points over unit cell, and reference to unique number in
%IBZ, for energy lookup in 'energy'
%in 'H' xtal system, these coords are reciprocal lattice vector based coords

%matrix cell is cuboid of dims (a*_x,2b*_y,c*z), all symmetry related cells on border included
%first generate cartesian coords of orthorhombic kpoints over matrix cell,
%in units of spanning vector components

[cartX cartY cartZ]=ndgrid(0:0.5:(dims(1)*2-1+1)*0.5,0:dims(2)*2-1+1,0:dims(3)-1+1);

%generate vector coords and uc_klist index for hexagonal mesh points within matrix cell
%first double up uc_klist along b* dirn
mc_hexpts=[uc_klist; [uc_klist(:,1:2) uc_klist(:,3)+LCMdivs uc_klist(:,4)]];
%map coords to those used for cartX,Y,Z
mc_hexpts(:,2)=mc_hexpts(:,2)/(LCMdivs/dims(1));
mc_hexpts(:,3)=mc_hexpts(:,3)/(LCMdivs/dims(2));
mc_hexpts(:,4)=mc_hexpts(:,4)/(LCMdivs/dims(3));
cart_mc_hexpts=[mc_hexpts(:,1) mc_hexpts(:,2)+0.5*mc_hexpts(:,3) mc_hexpts(:,3) mc_hexpts(:,4)];

%wrap cart coords back into matrix cell
cart_mc_hexpts(:,2)=mod(cart_mc_hexpts(:,2),dims(1));

%add symm related points around borders of matrix cell (so interpolation works)
x_border=cart_mc_hexpts(:,2)==0;
cart_mc_hexpts=[cart_mc_hexpts; ...
    [cart_mc_hexpts(x_border,1) cart_mc_hexpts(x_border,2)+dims(1) cart_mc_hexpts(x_border,3:4)]];
y_border=cart_mc_hexpts(:,3)==0;
cart_mc_hexpts=[cart_mc_hexpts; ...
    [cart_mc_hexpts(y_border,1:2) cart_mc_hexpts(y_border,3)+2*dims(2) cart_mc_hexpts(y_border,4)]];
z_border=cart_mc_hexpts(:,4)==0;
cart_mc_hexpts=[cart_mc_hexpts; ...
    [cart_mc_hexpts(z_border,1:3) cart_mc_hexpts(z_border,4)+dims(3)]];

%interpolate Wien2k energies to matrix cell
cartX=permute(cartX,[2 1 3]); cartY=permute(cartY,[2 1 3]); cartZ=permute(cartZ,[2 1 3]);
for bandnum=1:size(energy,1)
    cartE{bandnum} = griddata3(cart_mc_hexpts(:,2),cart_mc_hexpts(:,3),cart_mc_hexpts(:,4),...
        energy(bandnum,cart_mc_hexpts(:,1)),cartX,cartY,cartZ,interpmethod);
    cartE{bandnum}=reshape(cartE{bandnum},size(cartX));
    cartE{bandnum}=permute(cartE{bandnum},[2 1 3]); 
end
cartX=permute(cartX,[2 1 3]); cartY=permute(cartY,[2 1 3]); cartZ=permute(cartZ,[2 1 3]);



%scale cart coords back to Anstrom^-1
cartX=cartX/max(cartX(1:end))*rl_vecmags(1);
cartY=cartY/max(cartY(1:end))*rl_vecmags(2)*2;
cartZ=cartZ/max(cartZ(1:end))*rl_vecmags(3);
