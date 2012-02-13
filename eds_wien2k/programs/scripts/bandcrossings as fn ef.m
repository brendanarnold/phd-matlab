EFs=0.3:0.0001:0.5;

FC=zeros(numel(EFs),2);
for n=1:numel(EFs)
    FC(n,1)=numel(find(bandranges{1}(:,1)<EFs(n) & bandranges{1}(:,2)>EFs(n)));
    FC(n,2)=numel(find(bandranges{2}(:,1)<EFs(n) & bandranges{2}(:,2)>EFs(n)));
    fullb(n,1)=numel(find(bandranges{1}(:,2)<EFs(n)));
    fullb(n,2)=numel(find(bandranges{2}(:,2)<EFs(n)));
end

find(FC(:,1)==3 & FC(:,2)==4);