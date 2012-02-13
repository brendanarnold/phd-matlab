function [data tits]=clip2mat(hdrlines,cols);
%function [data tits]=clip2mat(hdrlines,cols);
%
%puts text contents of clipboard into matrix 'data'

str=clipboard('paste');
dfmtstr=[repmat(['%f\t'],1,cols-1) '%f\n'];
data=textscan(str,dfmtstr,1e7,'headerLines',hdrlines);
data=cell2mat(data);
if nargout>1
    tfmtstr=[repmat(['%s\t'],1,cols-1) '%s\n'];
    tits=textscan(str,tfmtstr,1);
end