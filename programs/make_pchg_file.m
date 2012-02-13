function make_pchg_file(bands,filename,atoms,bandnums);
%function make_pchg_file(bands,filename,atomnos,bandnums);
%
%writes a file in the format of a .energy file but with partial chg info rather than energies
%so that fexp.f can read it and perform Fourier interpolation
%filename is name of source energy file, dest pchg file is filename.pchg
%bandnums contains indexes into bands.Wien2kbandnums for the bands to write

%atoms=[3 5 8]; for Y124
%atoms=[3 4 7 11 12]; for Y124

v=squeeze(sum(bands.pchgs{1}(:,:,atoms,1),3));

notatend=true; inbands=false; indexline=false;
fid=fopen(filename,'r');
fid2=fopen([filename '.pchg'],'w');
while notatend
    linestr=fgetl(fid);
    if linestr==-1
        notatend=false;
        continue;
    end
    if ~inbands & ~isempty(strfind(lower(linestr),'e')) & numel(linestr)>40
        indexline=true;
    end
    if indexline        
        data=textscan([linestr(1:19) ' ' linestr(20:38) ' ' linestr(39:end)],'%19f%19f%19f%10s%6f%6d%6f');
        if numel(data{4})<1
            debug=true;
        end
        kpt=str2num(data{4}{1});
        kptnbands=data{6};
        fprintf(fid2,[linestr '\r\n']);
%        fprintf(fid2,'%19.12f%19.12f%19.12f%10s%6i%6i%f',data);
        indexline=false;
        inbands=true;
        b=1;
        continue;
    elseif inbands 
        data=textscan(linestr,'%12d%f\r\n');
        [ismem ind]=ismember(b,bands.Wien2k_bandnums(bandnums));
        if ismem
            newv=v(ind,kpt);
        else
            newv=0;
        end
        fprintf(fid2,'%12i %19.16f\r\n',b,newv);
        b=b+1;
        if b==1+kptnbands
            indexline=true;
            continue;
        end
    else
        fprintf(fid2,[linestr '\r\n']);
    end
end
fclose(fid);
fclose(fid2);