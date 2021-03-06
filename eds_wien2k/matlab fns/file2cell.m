function output=file2cell(pathfilename,headerlines,readtitles,delimstr,showprogress);
%function output=file2cell(pathfilename,headerlines,readtitles,delimstr,showprogress);
%
%uses delimstr2cell to convert delimited text file to cell array, line by line
%damn slow!

infile=fopen(pathfilename,'r');
if infile==-1
    disp('File not found.');
    output=[];
    return;
end

for r=1:headerlines
    linestr=fgetl(infile);
end
if readtitles
    titlelinestr=fgetl(infile);
    output.titles=delimstr2cell(titlelinestr,delimstr);
end

notatend=true; rownum=1;
while notatend
    linestr=fgetl(infile);
    if linestr==-1
        notatend=false;
        continue;
    end
    rowchar=char(delimstr2cell(linestr,delimstr)');
    for cnum=1:size(rowchar,1)
        fieldval=str2num(rowchar(cnum,:));
        if isempty(fieldval)
            rowdata(1,cnum)=nan;
        else
            rowdata(1,cnum)=fieldval;
        end
    end
    if rownum==1 | numel(rowdata)==size(output.data,2)
        output.data(rownum,:)=rowdata;
    else
        debug=true;
    end
    rownum=rownum+1;
    if showprogress & rownum==(1000*round(rownum/1000))
        disp(['Done ' num2str(rownum) ' rows.']);
    end
end