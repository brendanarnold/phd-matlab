function energy=importenergy(infilename,bandnums,numkpts);
%function energy=importenergy(infilename);
%
%reads Wien2k case.energy(up/dn) file and returns result in a 2d array
%energy(bandnum,kpointnum)
%if bandnums is empty, loads all bands
infile=fopen(infilename,'r');
if infile==-1
    disp('case.energy file not found');
    energy=[];
    return;
end
disp(['Importing energies from ' infilename]);

loadallbands=false;
if isempty(bandnums)
    loadallbands=true;
else
    loadallbands=false;
end

energy=zeros(length(bandnums),numkpts);
formatinfo=load('formatinfo.mat');
numfields=length(formatinfo.indexfieldpos)-1;

notatend=true; kptnum=0; gotkptinfo=false; linenum=1;
while notatend
    linestr=fgetl(infile);
    if linestr==-1
        disp([num2str(numkpts) ' kpoints expected. Energy file only contained ' num2str(kptnum)]);
        energy=[];
        return;
    end
    %test whether line is index line for first kpoint
    if ismember(length(linestr),formatinfo.indexlinelengths)
        kptnum=kptnum+1;
        for fnum=1:numfields
            str_buffer=linestr(formatinfo.indexfieldpos(fnum)+1+formatinfo.skipleadchars(fnum):formatinfo.indexfieldpos(fnum+1));
            value=str2num(str_buffer);
            if numel(value)==1
                fieldval(fnum)=value;                
            end
        end
        IBZkptnum=fieldval(4);
        numbands=fieldval(6);
        gotkptinfo=true;
    end
    if gotkptinfo
        cband=1;
        linestr=fgetl(infile);
        datalinelength=length(linestr)+1;        
        if loadallbands    
            bandnums=1:numbands;
        end
        for bandnum=1:numel(bandnums)
            for n=1:(bandnums(bandnum)-cband-1) %note this goes wrong when reading band 1 i.e. cannot read -ve no. of lines
                linestr=fgetl(infile);
            end
            linestr=fgetl(infile);
            checkbandnum=str2num(linestr(1:13));
            if checkbandnum~=bandnums(bandnum)
                disp('Problem!');                
                return;
            end
            cband=bandnums(bandnum)+1;                    
            energy(bandnum,IBZkptnum)=str2num(linestr(14:end));    
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



