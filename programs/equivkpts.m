function [equivpts weights]=equivkpts(kpts,symmops,LCMdivs);
%function result=equivkpts(kpts,symmops,LCMdivs);
%
%kpts is an n x 3 array containing kpts in rl coords
%symmops is a 3 x 3 x m array of symmetry operations (from outputkgen)
%LCMdivs is denominator in .klist coords
%generates symm-equivalent kpoints and counts how many in first unit cell
for n=1:size(symmops,3)
    equivpts(:,:,n)=kpts*squeeze(symmops(:,:,n)');
end
equivpts=mod(equivpts,LCMdivs);
for r=1:size(equivpts,1)
    [upts uindices]=unique(squeeze(permute(equivpts(r,:,:),[3 2 1])),'rows');
    weights(r)=numel(uindices);
    equivpts(r,:,:)=cat(3,equivpts(r,:,uindices),repmat([nan nan nan],[1 1 size(equivpts,3)-numel(uindices)]));
end