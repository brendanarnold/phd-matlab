function [pcharges coltitles]=importQTL(infilename,bandnums,numkpoints);
%function energy=importenergy(infilename);
%
%reads Wien2k case.energy(up/dn) file and returns result in a 2d array
%energy(bandnum,kpointnum)
infile=fopen(infilename,'r');
if infile==-1
    disp('case.qtl file not found');
    pcharges=[];
    return;
end
disp(['Importing partial charges from ' infilename '. Bands ' num2str(bandnums') '. Expecting ' num2str(numkpoints) ' kpoints.']);

%skip header lines
stillinheader=true; linenum=0;
while stillinheader
        linestr=fgetl(infile);
        linenum=linenum+1;
        %test whether line is index line for first kpoint        
        if length(linestr)>=9 & strcmp(linestr(1:9),' JATOM  1');
            stillinheader=false;
        end        
end

%read partial chg titles for each atom
notatend=true; atomnum=1;
while notatend
    titles{atomnum}=delimstr2cell(linestr(32:end),',');
    linestr=fgetl(infile);
    linenum=linenum+1;
    if ~strcmp(linestr(1:6),' JATOM')
        notatend=false;
    else
        atomnum=atomnum+1;
        datastart=ftell(infile);
    end
end
coltitles=cat(2,titles{:});
numatoms=atomnum;
for atomnum=1:numatoms
    numcols(atomnum)=length(titles{atomnum});
end
numcols(atomnum+1)=1;

[FL,p]=grep2('-s','BAND',infilename);
byteoffs=p.byteoffs; clear p,FL;

numbandsinfile=numel(byteoffs);
disp(['QTL file contains ' num2str(numbandsinfile) ' bands.']);
if isempty(bandnums)
    bandnums=1:numbandsinfile;
end

fseek(infile,byteoffs(1)-1,'bof');
linestr=fgetl(infile);
if ~strcmp(linestr,'BAND:   1')
    disp('Unexpected format in QTL file');
    return;
end

disp(['Num atoms=' num2str(numatoms) '. Num pchgs=' num2str(numel(titles{1})) ': ' sprintf('%s ',titles{1}{:})]);
%get data from bands in bandnums
for bandnum=1:length(bandnums)
    fseek(infile,byteoffs(bandnums(bandnum))-2,'bof');
    linestr=fgetl(infile);
    if ~strcmp(linestr(1:5),' BAND')
        disp('Unexpected format in QTL file');
        return;
    end
    for knum=1:numkpoints
    for atomnum=1:numatoms+1
        linestr=fgetl(infile);
        temp=cell2mat(textscan(linestr,'%f'));
        pcharges(bandnum,knum,atomnum,1:numcols(atomnum))=temp(3:end);
    end    
    end   
    disp(['Read band ' num2str(bandnums(bandnum))]);
end
