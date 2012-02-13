function [X,Y,Z,E]=tetraE(bands,bandnums,kps,nrms,inclpt,clip,plotpts)
%function [X,Y,Z,E]=tetraE(bands,bandnums,kps,nrms,inclpt,clip,plotpts)
%
%makes a matrix of E on original W2k mesh points 
%kspace volume covered is specified by bounding planes defined by nrms and kps

if isempty(bandnums)
    bandnums=1:numel(bands.Wien2k_bandnums);
end
if isempty(kps)
    [kps nrms]=normdir(nrms);
end

dirnnum=1; tol_r=1e-4;
%find Cart coords of smallest cuboid that includes all vertices of plotting vol defined by
%bplane_normals
[bfaces bverts]=vol_vertices(kps,nrms,inclpt);
rl_bverts=bverts*inv(bands.spanning_vecs);
exts=[min(rl_bverts(:,1)) min(rl_bverts(:,2)) min(rl_bverts(:,3)); ...
    max(rl_bverts(:,1)) max(rl_bverts(:,2)) max(rl_bverts(:,3))];

[a b c]=ndgrid(floor(exts(1,1)):ceil(exts(2,1)),floor(exts(1,2)):ceil(exts(2,2)),floor(exts(1,3)):ceil(exts(2,3)));

uc_abc=[mod(a(1:end)',bands.W2kdims(1))*bands.LCMdivs/bands.W2kdims(1) ...
    mod(b(1:end)',bands.W2kdims(2))*bands.LCMdivs/bands.W2kdims(2) ...
    mod(c(1:end)',bands.W2kdims(3))*bands.LCMdivs/bands.W2kdims(1)];
[ismem index]=ismember(uc_abc,bands.uc_klist(:,2:4),'rows');
if sum(ismem)~=numel(a)
    disp('Some points not on tetra mesh - something wrong.');
end
XYZ=[a(1:end)' b(1:end)' c(1:end)']*bands.spanning_vecs;
out=isoutside(kps,nrms,inclpt,XYZ,[]);
disp(['Removing ' num2str(sum(out)) ' pts outside vol. ' num2str(sum(~out)) ' pts remaining.']);
for bandnum=1:numel(bandnums)
    E{bandnum}=bands.IBZenergy{dirnnum}(bandnums(bandnum),bands.uc_klist(index,1));
    if clip
        E{bandnum}(out)=nan;
    end
    E{bandnum}=reshape(E{bandnum},size(a));
end
X=reshape(XYZ(:,1),size(a)); Y=reshape(XYZ(:,2),size(a)); Z=reshape(XYZ(:,3),size(a));

if plotpts
    plot3(X(~isnan(E{1}(1:end))),Y(~isnan(E{1}(1:end))),Z(~isnan(E{1}(1:end))),'.')
    hold on;
    plotvol(gca,kps,nrms,[0 0 0]);
end

  
