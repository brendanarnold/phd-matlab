function mat2file(data,tits,fname);
%function mat2file(data,tits,fname);
%saves array 'data' to a tab delim text file, writes titles 'tits' if not empty

fid=fopen(fname,'w');
if ~isempty(tits)
for n=1:numel(tits)-1
    fprintf(fid,'%s\t',tits{n});
end
fprintf(fid,'%s\r\n',tits{end});
end

for n=1:size(data,1)
    for m=1:size(data,2)-1
        fprintf(fid,'%.8e\t',data(n,m));
    end
    fprintf(fid,'%.8e\r\n',data(n,end));
end
fclose(fid);
