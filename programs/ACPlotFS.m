function ACPlotFS(fs);

pvol=fs.BZfacenormals*fs.rlvs;
allcolors=[0 0 1;0 1 0;1 0 0;1 1 0;0 1 1;1 0 1];
allcolors=[allcolors;allcolors*0.5];
sheets=max(size(fs.Wien2k_bandnums));
colors=allcolors(1:sheets,:);

plot_FS(fs,[],fs.FermiLevel,colors,[],pvol,[],[],[],'linear',gca,true);
plotvol(gca,[],pvol,[0 0 0]);



