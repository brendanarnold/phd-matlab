function repackso(inputfilename,outputfilename);

fs=load(inputfilename);

fsp.pathcasename=fs.pathcasename;
fsp.latt_type=fs.latt_type;
fsp.latt_params=fs.latt_params;
fsp.latt_angles=fs.latt_angles;
fsp.symmops=fs.symmops;
fsp.rlvs=fs.rlvs;
fsp.rec_latt_vecs=fs.rlvs;
fsp.n_rlvs=fs.n_rlvs;
fsp.W2kdims=fs.W2kdims;
fsp.spanning_vecs=fs.spanning_vecs;
fsp.IBZ_klist=fs.IBZ_klist;
fsp.LCMdivs=fs.LCMdivs;
fsp.FermiLevel=fs.FermiLevel;
fsp.energyfname=fs.energyfname;
fsp.rlv_mags=fs.rlv_mags;
fsp.uc_klist=fs.uc_klist;
fsp.mcell=fs.mcell;
fsp.cuboidextents=fs.cuboidextents;
fsp.cart_submesh_shift=fs.cart_submesh_shift;
fsp.cartX=fs.cartX;
fsp.cartY=fs.cartY;
fsp.cartZ=fs.cartZ;
fsp.dx=fs.dx;
fsp.dy=fs.dy;
fsp.dz=fs.dz;
fsp.cartdims=fs.cartdims;
fsp.BZfacenormals=fs.BZfacenormals;
fsp.searchvol=fs.searchvol;
fsp.centresvol=fs.centresvol;
fsp.ptinvol=fs.ptinvol;

cols=size(fs.cartE,2)/2;
energies=[];
for i=1:cols
  fsp.cartE(i)=fs.cartE(2*i-1);  
  fsp.Wien2k_bandnums(i)=fs.Wien2k_bandnums(2*i-1);
  fsp.spindirns(i)=fs.spindirns(2*i-1);
  energies=[energies;fs.IBZenergy{1}(i,:)];
end;
  fsp.IBZenergy{1}=energies;

save(outputfilename,'-struct','fsp');