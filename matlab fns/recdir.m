function filelist=recdir(initpath,filename);
%function filelist=recdir(initpath,filename);
%
%same as dir but recurses subdirs

filelist=dir([initpath filename]);
tcelllist=struct2cell(filelist);
filelist=filelist(~strcmp(tcelllist(1,:),'.') & ~strcmp(tcelllist(1,:),'..'));
cfilelist=dir([initpath '*.*']);
tcelllist=struct2cell(cfilelist);
cfilelist=cfilelist(~strcmp(tcelllist(1,:),'.') & ~strcmp(tcelllist(1,:),'..'));
dirlist=cfilelist(find([cfilelist.isdir]));
for dirnum=1:length(dirlist)
    deeperfiles=recdir([initpath dirlist(dirnum).name '\'],filename);
    for fnum=1:length(deeperfiles)
        deeperfiles(fnum).name=[dirlist(dirnum).name '\' deeperfiles(fnum).name];
    end
    filelist=[filelist; deeperfiles];    
end
