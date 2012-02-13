function expand_klist(pathcasename,outfilename,latt_type);
%function expand_klist(pathfilename,outfilename,latt_type);
%
%takes klist file as generated by Wien2k kgen and adds points to make complete orthorhombic mesh
%note latt_type is Bravais latt type not rec. latt type
%needs .klist, .outputkgen

[symmops rec_latt_vecs n_rec_latt_vecs cuboid]=getsymmops([pathcasename '.outputkgen']);

infile=fopen([pathcasename '.klist']);
[inpathstr,infilename,ext,versn] = fileparts(pathcasename);
outfile=fopen([inpathstr outfilename],'w');

linestr=fgetl(infile);
data=sscanf(linestr(1:35),'%f',6);
dims_start=findstr(linestr,'(')+1;
dims_end=findstr(linestr,')')-1;
dims=sscanf(linestr(dims_start:dims_end),'%f',3);
LCMdivs=vlcm(dims);

[nx ny nz]=ndgrid(0:LCMdivs/dims(1):cuboid(1)*LCMdivs-LCMdivs/dims(1),0:LCMdivs/dims(2):cuboid(2)*LCMdivs-LCMdivs/dims(2),...
    0:LCMdivs/dims(3):cuboid(3)*LCMdivs-LCMdivs/dims(3));
kpts=[reshape(nx,[],1) reshape(ny,[],1) reshape(nz,[],1)];
rl_kpts=kpts*inv(n_rec_latt_vecs);
rl_kpts=mod(rl_kpts,LCMdivs);

[uc_rl_kpts weights]=equivkpts(rl_kpts,symmops,LCMdivs);
%weights counts all equiv points. also calc weights2 which includes only points on orthomesh
weights2=zeros(size(uc_rl_kpts,1),1);
for r=1:size(uc_rl_kpts,1)
    for epn=1:size(uc_rl_kpts,3)
        ep=uc_rl_kpts(r,:,epn)*n_rec_latt_vecs;
        if any(isnan(ep))
            break;
        end
        if isequal(ep,round(ep))
            weights2(r)=weights2(r)+1;
        end
    end
end
%only keep inequiv points
uc_rl_kpts(isnan(uc_rl_kpts))=-1; %needed for uniqueness test since nan~=nan
[u_uc_rl_kpts u_indices]=unique(reshape(uc_rl_kpts,[size(uc_rl_kpts,1) 3*size(uc_rl_kpts,3)]),'rows');
u_uc_rl_kpts(u_uc_rl_kpts==-1)=nan;
u_uc_c_kpts=rl_kpts(u_indices,1:3)*n_rec_latt_vecs;
for kptnum=1:numel(u_indices)
    if kptnum==1
        remstr=linestr(36:end);
    else
        remstr='';
    end
    writeline(outfile,[kptnum u_uc_c_kpts(kptnum,:) LCMdivs weights(u_indices(kptnum))],remstr);
end
fclose(infile);
fprintf(outfile,'%s\n\n','END');
fclose(outfile);
disp([num2str(numel(u_indices)) ' IBZ kpoints in list. Total orthomesh weight in unit cell: ' num2str(sum(weights2(u_indices)))]);

function writeline(outfile,params,remstr);

newstr=sprintf('    % 6.0f% 5.0f% 5.0f% 5.0f% 5.0f% 5.1f ',params);
fprintf(outfile,'%s%s\n',newstr,remstr);


    