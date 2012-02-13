function [cartX cartY cartZ cartE]=orthomeshe(energy,uc_klist,dims,LCMdivs,latt_type,n_rec_latt_vecs,rl_vecmags);

%expects energy matrix to have one band per row
switch latt_type
    case 'H'
        orthomesh_uc=[1 2 1];
    case {'F','C','P','B'}
        orthomesh_uc=[2 2 2];
end

switch latt_type   
case 'H'
    %for hexagonal, interpolate orthorhombic mesh coords into sheared tetrahedral coords
    [a b c]=ndgrid(0:dims(1),0:dims(2),0:dims(3));
    a=permute(a,[2 1 3]); b=permute(b,[2 1 3]); c=permute(c,[2 1 3]);    
    [x y z]=ndgrid(0:orthomesh_uc(1)*dims(1),0:orthomesh_uc(2)*dims(2),0:orthomesh_uc(3)*dims(3));
    x=permute(x,[2 1 3]); y=permute(y,[2 1 3]); z=permute(z,[2 1 3]);
    cart=[x(1:end)' y(1:end)' z(1:end)']*[1 0 0; -0.5 1 0; 0 0 1];
    cartx=mod(cart(:,1),dims(1)); carty=mod(cart(:,2),dims(2)); cartz=mod(cart(:,3),dims(3));
    shearmethod=false;
    for bandnum=1:size(energy,1)
        uc_energy=energy(bandnum,uc_klist(:,1));
        uc_energy=reshape(uc_energy,dims(3),dims(2),dims(1));
        uc_energy=permute(uc_energy,[3 2 1]);
        uc_energy=addborders(uc_energy,false);
        if shearmethod
            cartE{bandnum}=interp3(a,b,c,permute(uc_energy,[2 1 3]),...
                cartx,carty,cartz);        
            cartE{bandnum}=reshape(cartE{bandnum},orthomesh_uc(2)*dims(2)+1,orthomesh_uc(1)*dims(1)+1,orthomesh_uc(3)*dims(3)+1);
            cartE{bandnum}=permute(cartE{bandnum},[2 1 3]);
        else
            rlmeshx=a*n_rec_latt_vecs(1,1)+b*n_rec_latt_vecs(2,1)+c*n_rec_latt_vecs(3,1);
            rlmeshy=a*n_rec_latt_vecs(1,2)+b*n_rec_latt_vecs(2,2)+c*n_rec_latt_vecs(3,2);
            rlmeshz=a*n_rec_latt_vecs(1,3)+b*n_rec_latt_vecs(2,3)+c*n_rec_latt_vecs(3,3);
            cartE{bandnum}=griddata3(rlmeshx,rlmeshy,rlmeshz,cartx,carty,cartz,permute(uc_energy,[2 1 3]));            
            cartE{bandnum}=reshape(cartE{bandnum},orthomesh_uc(2)*dims(2)+1,orthomesh_uc(1)*dims(1)+1,orthomesh_uc(3)*dims(3)+1);
            cartE{bandnum}=permute(cartE{bandnum},[2 1 3]);
        end
    end
case {'F','C','P','B'}  
    [cartX cartY cartZ]=ndgrid(0:orthomesh_uc(1)*dims(1),0:orthomesh_uc(2)*dims(2),0:orthomesh_uc(3)*dims(3));
    rl_cub=[cartX(1:end)' cartY(1:end)' cartZ(1:end)']*inv(n_rec_latt_vecs);
    rl_cub=mod(rl_cub,repmat(dims',size(rl_cub,1),1));
    rl_cub=rl_cub.*repmat(LCMdivs./dims',size(rl_cub,1),1);
    cub_serial=rl_cub(:,1)+rl_cub(:,2)*1e3+rl_cub(:,3)*1e6;
    uc_serial=uc_klist(:,2)+uc_klist(:,3)*1e3+uc_klist(:,4)*1e6;
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

[cartX cartY cartZ]=ndgrid(rl_vecmags(1)*(0:orthomesh_uc(1)*dims(1))/dims(1),rl_vecmags(2)*(0:orthomesh_uc(2)*dims(2))/dims(2),...
    rl_vecmags(3)*(0:orthomesh_uc(3)*dims(3))/dims(3));
