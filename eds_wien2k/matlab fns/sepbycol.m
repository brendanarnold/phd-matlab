function result=sepbycol(fileorvar,titlesorhdrlines,numcols,delim,col1,col2,labels);
%function result=sepbycol(fileorvar,titlesorhdrlines,numcols,delim,col1,col2,labels);
%
% reads data from pathfilename and writes a new .sep file in which values
% of col2 are separated into different cols according to value of col1

fileout=false;
if isstr(fileorvar)
    fileout=true;
    pathfilename=fileorvar;
    headerlines=titlesorhdrlines;
    infile=fopen(pathfilename);
    [pathstr,name,ext,versn] = fileparts(pathfilename);
    outfile=fopen([pathstr name '.sep' ext],'w');

    for r=1:headerlines
        linestr=fgetl(infile);
        fprintf(outfile,'%s\r\n',linestr);
    end

    %read col titles
    linestr=fgetl(infile);
    coltitles=delimstr2cell(linestr,delim);

    %read data
    notatend=true;
    linenum=1;
    while notatend
        linestr=fgetl(infile);
        if ischar(linestr) & length(linestr)>0
            data(linenum,:)=delimstr2cell(linestr,delim);
            linenum=linenum+1;
        else
            notatend=false;
        end
    end
else
    data=fileorvar;   
    coltitles=titlesorhdrlines;
    numcols=size(data,2);
end

%find unique values of col1
[colvalues ind1 ind2]=unique(data(:,col1));
numvalues=length(colvalues);
disp(['There are ' num2str(numvalues) ' distinct values of col ' num2str(col1)]);

%generate new col titles
if exist(coltitles) && ~isempty(coltitles)
    for counter=1:numvalues
        newlabels{counter}=[coltitles{col2} ' (' coltitles{col1} '=' num2str(colvalues{counter}) ')'];
    end
    if fileout
        %write col titles
        for colnum=1:numcols
            fprintf(outfile,'%s%s',coltitles{colnum},delim);
        end
        for colnum=1:numvalues
            fprintf(outfile,'%s',newlabels{colnum});
            if colnum==numvalues
                fprintf(outfile,'\r\n');
            else
                fprintf(outfile,'%s',delim);
            end
        end
    end
end

if fileout
    %write data
    numrows=size(data,1);
    for rownum=1:numrows    
        for colnum=1:numcols
            fprintf(outfile,'%s%s',data{rownum,colnum},delim);
        end
        for colnum=1:numvalues
            if colnum==ind2(rownum)
                splitdata(colnum)=data(rownum,col2);
            else
                splitdata(colnum)={' '};
            end
            fprintf(outfile,'%s',splitdata{colnum});
            if colnum==numvalues
                fprintf(outfile,'\r\n');
            else
                fprintf(outfile,'%s',delim);
            end
        end
    end
else
    if exist(coltitles) && ~isempty(coltitles)
        result.titles=[coltitles newlabels];
    end
    result.data=[data nan(size(data,1),numvalues)];
    for rownum=1:size(data,1)
        result.data(rownum,numcols+ind2(rownum))=data(rownum,col2);
    end
end