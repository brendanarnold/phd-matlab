function output=delimstr2cell(inputstr,delim);
%delimstr2cell(inputstr,delim);
%
%
if ~iscell(delim)
    delim=delim(1);
end


if strcmp(delim,' ')
    %treated differently as consecutive delims count as one, initial/final delims are ignored
    notatend=true;
    strcounter=1;
    remstr=inputstr;
    while notatend   
        [output{strcounter},remstr]=strtok(remstr,delim);
        strcounter=strcounter+1;
        if length(remstr)<1
            notatend=false;
        end
    end
else    
    delimpos=strfind(inputstr,delim);
    if isempty(delimpos)
        output{1}=inputstr;
    else
        fieldstarts=min([1 delimpos+1],length(inputstr));
        fieldends=max([delimpos-1 length(inputstr)],1);
        for n=1:length(delimpos)+1
            output{n}=inputstr(fieldstarts(n):fieldends(n));
        end
        if delimpos(1)==1
            output{1}='NaN';
        end
        if delimpos(end)==length(inputstr)
            output{length(delimpos)+1}='NaN';
        end
    end
end
        

