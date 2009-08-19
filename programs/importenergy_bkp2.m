function energy=importenergy(infilename,bandnums,numkpts);
%function energy=importenergy(infilename);
%
%reads Wien2k case.energy(up/dn) file and returns result in a 2d array
%energy(bandnum,kpointnum)
infile=fopen(infilename,'r');
if infile==-1
    disp('case.energy file not found');
    energy=[];
    return;
end
disp(['Importing energies from ' infilename]);

energy=zeros(length(bandnums),numkpts);
indexfieldpos=[0 19 38 57 67 73 79 84];
numfields=length(indexfieldpos)-1;
indexlinelengths=[84 87];
datalinelength=39;

notatend=true; kptnum=0; gotkptinfo=false; linenum=1;
while notatend
    linestr=fgetl(infile);
    if linestr==-1
        disp([num2str(numkpts) ' kpoints expected. Energy file only contained ' num2str(kptnum)]);
        energy=[];
        return;
    end
    %test whether line is index line for first kpoint
    if ismember(length(linestr),indexlinelengths)
        kptnum=kptnum+1;
        for fnum=1:numfields
            value=str2num(linestr(indexfieldpos(fnum)+1:indexfieldpos(fnum+1)));
            if numel(value)==1
                fieldval(fnum)=value;                
            end
        end
        numbands=fieldval(6);
        gotkptinfo=true;
    end
    if gotkptinfo
        cband=1;
        for bandnum=1:length(bandnums)
            if bandnums(bandnum)<=numbands
                fseek(infile,datalinelength*(bandnums(bandnum)-cband),'cof');
                linestr=fgetl(infile);
                %following if added to cope with header lines in middle of file when parallel run files
                %are merged
                if ismember(length(linestr),indexlinelengths)
                    checkbandnum=str2num(linestr(1:13));
                    if checkbandnum~=bandnums(bandnum)
                        disp('Problem!');                
                        return;
                    end
                    energy(bandnum,kptnum)=str2num(linestr(14:end));
                    cband=bandnums(bandnum)+1;                    
                end
            end
        end
        if kptnum==numkpts
            notatend=false;
        else
            status=fseek(infile,datalinelength*(numbands-bandnums(end)),'cof');
            if status==-1
                disp('Problem!');
            end
        end
        gotkptinfo=false;
    end
end



