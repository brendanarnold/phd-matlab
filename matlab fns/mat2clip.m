function mat2clip(m,titles);

ncols=0;
if ~iscell(m)
    m={m};
end

colind1=[]; colind2=[];
for n=1:numel(m)
    colind1(end+1:end+size(m{n},2))=n;
    colind2(end+1:end+size(m{n},2))=1:size(m{n},2);
    numrows(n)=size(m{n},1);
    ncols=ncols+size(m{n},2);
end
maxrows=max(numrows);

outstr='';
if nargin==2
    for a=1:numel(titles)
        outstr=[outstr titles{a}];
        if a<numel(titles)
            outstr=[outstr char(9)];    
        end
    end
    outstr=[outstr char(10)];
end

for l=1:maxrows
    for c=1:ncols
        if l<=size(m{colind1(c)},1) && ~isnan(m{colind1(c)}(l,colind2(c)))
            outstr=[outstr num2str(m{colind1(c)}(l,colind2(c)),'%1.6e')];
        else
            outstr=[outstr ''];
        end
        if c<ncols
            outstr=[outstr char(9)];    
        end
    end
    outstr=[outstr char(10)];
end

clipboard('copy',outstr);