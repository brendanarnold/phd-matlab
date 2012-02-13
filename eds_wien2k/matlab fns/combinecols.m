function newdata=combinecols(filename,numhdrlines,delimstr,ccols);
%function newdata=combinecols(filename,numhdrlines,delimstr,ccols);
%
%reads delim text file, compresses all of ccols into single column, making data longer and narrower


if isempty(filename)
    [filename pathname]=uigetfile('*.*');
    filename=[pathname filename];
end
infile=fopen(filename,'r');

for hlnum=1:numhdrlines
    hdrlinestr(hlnum)=fgetl(infile);
end

coltitleslinestr=fgetl(infile);
coltitles=delimstr2cell(coltitleslinestr,delimstr);

notatend=true; rownum=1;
while notatend
    linestr=fgetl(infile);
    if linestr==-1
        notatend=false;
    else
        rowcelldata=delimstr2cell(linestr,delimstr);
        rowcelldata(strcmp(rowcelldata,''))={'NaN'};
        rowdata=str2num(char(rowcelldata))';
        data(rownum,:)=rowdata;
        rownum=rownum+1;
    end
end
ncols=size(data,2);
nrows=size(data,1);
disp(['Read table:' num2str(nrows) 'x' num2str(ncols)]);

%make new table with ccols combined into single column
newdata=[];
isccol=ismember(1:ncols,ccols);
isfirstccol=(1:ncols)==ccols(1);
for colnum=1:length(ccols)
    nonemptyrows=~isnan(data(:,ccols(colnum)));
    newrows=data(nonemptyrows,~isccol | isfirstccol);
    newrows(:,isfirstccol)=data(nonemptyrows,ccols(colnum));
    newdata=[newdata;  newrows;];
end    
    
%write data to .cmb file
[pathstr,name,ext,versn]=fileparts(filename);
outfilename=[pathstr name '.cmb' ext];
outfile=fopen(outfilename,'w');
for hlnum=1:numhdrlines
    fprintf(outfile,'%s\r\n',hdrlinestr(r));
end
for cnum=1:ncols
    if ~ismember(cnum,ccols)
        fprintf(outfile,'%s',coltitles{cnum});
    elseif cnum==ccols(1)
        fprintf(outfile,'Combined');
    end
    if cnum<ncols
        fprintf(outfile,delimstr);
    else
        fprintf(outfile,'\r\n');
    end
end
for rnum=1:size(newdata,1)
for cnum=1:size(newdata,2);
    if isnan(newdata(rnum,cnum))
        fprintf(outfile,' ');
    else
        fprintf(outfile,'%f',newdata(rnum,cnum));
    end
    if cnum==size(newdata,2)
        fprintf(outfile,'\r\n');
    else
        fprintf(outfile,'%s',delimstr);
    end
end
end
fclose(outfile);
disp(['Written ' num2str(size(newdata,1)) 'x' num2str(size(newdata,2)) ' data to ''' outfilename '''']);
    