function result=addborders(matrixorlist,add);
if ~iscell(matrixorlist)
    matrixlist{1}=matrixorlist;
else
    matrixlist=matrixorlist;
end
for mnum=1:length(matrixlist)
    result{mnum}=cat(1,matrixlist{mnum},matrixlist{mnum}(1,:,:));
    result{mnum}=cat(2,result{mnum},result{mnum}(:,1,:));
    result{mnum}=cat(3,result{mnum},result{mnum}(:,:,1));
    if add
        result{mnum}(:,:,end)=result{mnum}(:,:,end)+(size(result{mnum},3)-1)*(result{mnum}(1,1,2)-result{mnum}(1,1,1));
        result{mnum}(:,end,:)=result{mnum}(:,end,:)+(size(result{mnum},2)-1)*(result{mnum}(1,2,1)-result{mnum}(1,1,1));
        result{mnum}(end,:,:)=result{mnum}(end,:,:)+(size(result{mnum},1)-1)*(result{mnum}(2,1,1)-result{mnum}(1,1,1));
    end
end
if ~iscell(matrixorlist)
    nresult=result{1};
    result=nresult;
end
    
