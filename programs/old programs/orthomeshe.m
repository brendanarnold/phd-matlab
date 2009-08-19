function [cartX cartY cartZ cartE]=orthomeshe(energy,uc_klist,dims,LCMdivs,latt_type,n_rec_latt_vecs,rl_vecmags);
%expects energy matrix to have one band per row

switch latt_type
    case {'F','C','P','B'}
        orthomesh_uc=[2 2 2];
end

[cartX cartY cartZ]=ndgrid(rl_vecmags(1)*(0:orthomesh_uc(1)*dims(1))/dims(1),rl_vecmags(2)*(0:orthomesh_uc(2)*dims(2))/dims(2),...
    rl_vecmags(3)*(0:orthomesh_uc(3)*dims(3))/dims(3));

switch latt_type   
case {'F','C','P','B'}  
    %generate unit-spaced Cartesian grid to cover whole orthomesh unit cell
    [cartX cartY cartZ]=ndgrid(0:orthomesh_uc(1)*dims(1),0:orthomesh_uc(2)*dims(2),0:orthomesh_uc(3)*dims(3));
    %set rl_cub to contain coords of above Cart grid in (spanning) rec_latt_vec units
    rl_cub=[cartX(1:end)' cartY(1:end)' cartZ(1:end)']*inv(n_rec_latt_vecs);
    rl_cub=mod(rl_cub,repmat(dims',size(rl_cub,1),1));
    rl_cub=rl_cub.*repmat(LCMdivs./dims',size(rl_cub,1),1);
    cub_serial=rl_cub(:,1)+rl_cub(:,2)*1e3+rl_cub(:,3)*1e6;
    uc_serial=uc_klist(:,2)+uc_klist(:,3)*1e3+uc_klist(:,4)*1e6;
    %all cub points must be in uc but not necessarily vice versa0
    [foundok index]=ismember(cub_serial,uc_serial);
    if ~all(foundok)
        disp('Problem!')
        return
    end   
    for bandnum=1:size(energy,1)
        cartE{bandnum}=energy(bandnum,uc_klist(index,1));
        cartE{bandnum}=reshape(cartE{bandnum},size(cartX));
    end
end
%calculate cart mesh that includes all points in unit cell

cartX=cartX/LCMdivs*rl_vecmags(1); cartY=cartY/LCMdivs*rl_vecmags(2); cartZ=cartZ/LCMdivs*rl_vecmags(3);