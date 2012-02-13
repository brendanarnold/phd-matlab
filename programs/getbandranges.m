function bandranges=getbandranges(casename,dirnstr);

infile=fopen([casename '.output2' dirnstr],'r');
if infile==-1
    disp('case.output2 file not found');
    bandranges=[];
    return;
end

notatend=true; inbandranges=false;
while notatend
    linestr=fgetl(infile);
    if ~inbandranges & linestr==-1
        bandranges=0;
        notatend=false;       
    elseif inbandranges & (length(linestr)<5 | isempty(strfind(linestr,'band')))
        notatend=false;
    elseif inbandranges & length(linestr)>5 & ~isempty(strfind(linestr,'band'))        
        data=textscan(linestr,'%s%f%f%f%f');
        bandnum=data{2}; emin=data{3}; emax=data{4}; occ=data{5};
        bandranges(bandnum,:)=[emin emax occ];
    elseif ~inbandranges & length(linestr)>11 & ~isempty(strfind(linestr,'Bandranges'))
        inbandranges=true;
    end
end
fclose(infile);
    


