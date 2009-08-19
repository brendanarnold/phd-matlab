function result=calcTD(bandsdata,bandnums,undersample_factor,spinpol);
%function result=calcTD(bandsdata,bandnums,undersample_factor,spinpol);
%
%calcs various band parameters incl: Fermi velocity, plasma frequency
planckc=6.626e-34; el_charge=1.602e-19; epsilon0=8.854e-12; Ryd=13.61;
if spinpol
    spinpolfac=2;
else
    spinpolfac=1;
end

cuboidvol=abs(bandsdata.cuboidextents(2,1)-bandsdata.cuboidextents(1,1))*...
    abs(bandsdata.cuboidextents(2,2)-bandsdata.cuboidextents(1,2))*...
    abs(bandsdata.cuboidextents(2,3)-bandsdata.cuboidextents(1,3));
BZvol=dot(bandsdata.rec_latt_vecs(1,:),cross(bandsdata.rec_latt_vecs(2,:),bandsdata.rec_latt_vecs(3,:)));
if isempty(bandnums)
    bandnums=1:length(bandsdata.Wien2k_bandnums)
end

if ~isempty(undersample_factor)
    [bandsdata.cartX,bandsdata.cartY,bandsdata.cartZ,bandsdata.cartE{bands(1)}] = reducevolume(bandsdata.cartX,bandsdata.cartY,bandsdata.cartZ,bandsdata.cartE{bands(1)},undersample_factor);
    for bandnum=2:length(bands)
        bandsdata.cartE{bands(bandnum)} = reducevolume(bandsdata.cartE{bands(bandnum)},undersample_factor);
    end
else
    undersample_factor=1;
end

for bandnum=1:length(bandnums)
    x=permute(bandsdata.cartX,[2 1 3]);
    y=permute(bandsdata.cartY,[2 1 3]);
    z=permute(bandsdata.cartZ,[2 1 3]);
    v=permute(bandsdata.cartE{bandnums(bandnum)},[2 1 3]);
    [faces,verts]=isosurface(x,y,z,v,bandsdata.FermiLevel(bandsdata.spindirns(bandnums(bandnum))),'verbose','noshare');
    %[cfaces,cverts,fcolors]=clippatch(faces,verts,bplane_normals,[0 0 0],clipintersectedfaces);
    p=patch('Faces',faces,'Vertices',verts,'Visible','off');
    %computes face normals from data and face normals of p equal to these 
    vertnormals=isonormals(x,y,z,v,verts);
    if any(isnan(vertnormals))
        numnans=numel(find(isnan(vertnormals)));
        disp(['Warning! Vertnormals contains NaNs (' num2str(numnans) ' out of ' num2str(numel(vertnormals)) ' cartesian comps)']);
    end
    vertnormals(isnan(vertnormals))=0;
    %define matrix of components of Fermi vel (m/s) for each vertex
    vertVfs=2*pi/planckc*Ryd*el_charge/1e10*vertnormals;
    %Fermi vel for face is mean of Fermi vel for individual vertices (?)
    faceVfs=(vertVfs(faces(:,1),:)+vertVfs(faces(:,2),:)+vertVfs(faces(:,3),:))/3;
    facemodVf=sqrt(faceVfs(:,1).^2+faceVfs(:,2).^2+faceVfs(:,3).^2);
    
    %facenorms contains normals of patch polygons, with magnitude equal to poly area
    facenorms=cross(verts(faces(:,2),:)-verts(faces(:,1),:),verts(faces(:,3),:)-verts(faces(:,1),:))/2;
    faceareas=sqrt(facenorms(:,1).^2+facenorms(:,2).^2+facenorms(:,3).^2);
    
    result.Vf(bandnum)=sum(faceareas)/sum(faceareas./facemodVf);
    result.Sf(bandnum)=(BZvol/cuboidvol)*sum(faceareas);
    result.DOS(bandnum)=(1/2/pi^2/planckc/spinpolfac)*el_charge*(BZvol/cuboidvol)*(8*pi^3/(BZvol*1e30))*1e20*sum(faceareas./facemodVf);
    for axis1=1:3        
        for axis2=axis1:3
            result.TD(bandnum,axis1,axis2)=sum(faceVfs(:,axis1).*faceVfs(:,axis2)./facemodVf.*faceareas);           
        end
        %wp is plasma frequency (rad/s)
        result.wp(bandnum,axis1)=sqrt((1/2/pi^2/planckc/spinpolfac*((el_charge)^2)/epsilon0*result.TD(bandnum,axis1,axis1)*1e20*(BZvol/cuboidvol)));
        %Ep is plasma energy (eV)
        result.Ep(bandnum,axis1)=1/el_charge*planckc/2/pi*result.wp(bandnum,axis1);
    end
end