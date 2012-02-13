function plot_Wien2kspag(klistfilename,energyfilename,bandnums,spinpol,FermiLevels,Elimits,symm_labels,bandcharcol,exchangeplot,varargin);
%function plot_Wien2kspag(casename,bandnums,spinpol,FermiLevels,...
%    Elimits,symm_inds,symm_labels,bandcharcol,exchangeplot,varargin);
%
%plots Wien2k spaghetti calc: reads case(+.energy) file given by casename, plots bandnums,
%symm_inds give the indices of symm points in the kpoint array and symm_labels their names
%provide FermiLevels in Ryd, but plot is in eV

dirnstr={'up','dn'};
if ~spinpol
    dirnstr={''};
end

klistfile=fopen(klistfilename,'r');
if klistfile==-1
    disp(['File not found: ' klistfilename]);
    return;   
end
%get index numbers and labels of special kpoints from klist
symm_inds=[]; file_symm_labels=[];
notatend=true; linenum=1;
while notatend
    linestr=fgetl(klistfile);
    if strcmp(linestr(1:3),'END')
        notatend=false;
        num_kpts=linenum-1;
        continue;
    end
    label=sscanf(linestr(1:5),'%s');
    if length(label)>0
        symm_inds(end+1)=linenum;
        file_symm_labels{end+1}=label;
    end
    linenum=linenum+1;
end
disp([num2str(length(symm_inds)) ' special k-pts, ' num2str(num_kpts) ' total k-pts.']);
if isempty(symm_labels)
    symm_labels=file_symm_labels;
end


for dirnnum=1:length(dirnstr)
    E(:,dirnnum,:)=importenergy([energyfilename dirnstr{dirnnum}],bandnums,num_kpts);
    E(:,dirnnum,:)=13.6*(E(:,dirnnum,:)-FermiLevels(dirnnum));
    if dirnnum==1
        if isempty(bandnums)
            bandnums=1:size(E,1);
        end
        numbands=length(bandnums);
        numkpoints=size(E,3);
    end
    if ~isempty(bandcharcol)
        [QTL(:,dirnnum,:,:) QTLcoltitles]=importQTL([casename '.qtl' dirnstr{dirnnum}],bandnums,numkpoints);
    end
end

if ~isempty(bandcharcol)
    %prepare for band character plot
    circle=[cos(linspace(0,2*pi,15))' sin(linspace(0,2*pi,15))'];
    maxbandchar=max(reshape(QTL(:,:,:,bandcharcol),1,[]));
    minbandchar=min(reshape(QTL(:,:,:,bandcharcol),1,[]));
    aspectratio=(Elimits(2)-Elimits(1))/numkpoints;
    kstep=round(numkpoints/200);
    disp(['Min(max) band characters for column ' num2str(bandcharcol) ' (' QTLcoltitles{bandcharcol} '): ' num2str(minbandchar) '(' num2str(maxbandchar) ')']);
end

if ~isempty(bandcharcol)
    bandcharstr=['Band character %' QTLcoltitles{bandcharcol}];
else
    bandcharstr=[''];
end
titlestr=sprintf('%s%s',['Spaghetti plot of WIEN2k ' energyfilename],bandcharstr);
figure('Name',titlestr);
if exchangeplot 
    subplot(2,1,1);
end
hold on;

colours=propval(varargin,'colours');
styles=propval(varargin,'styles');
    
numdirnsforplot=1+(spinpol==1); curvenum=1;
for bandnum=1:numbands
for dirnnum=1:numdirnsforplot
    if ~isempty(colours)
        colchar=colours(curvenum);
    end
    if ~isempty(styles)
        stylechar=styles(curvenum);
    end
    style=[stylechar colchar];
    
    if ~isempty(bandcharcol)
        bandcharvals=(Elimits(2)-Elimits(1))/25*(reshape(QTL(bandnum,dirnnum,:,bandcharcol),1,[])-minbandchar)/(maxbandchar-minbandchar);
        for knum=1:kstep:numkpoints
            hplot(curvenum)=line([knum; knum],[E(bandnums(bandnum),dirnnum,knum)-bandcharvals(knum); ...
                E(bandnums(bandnum),dirnnum,knum)+bandcharvals(knum)],'Color',colchar);
        end
    else
        hplot(curvenum)=plot(1:numkpoints,reshape(E(bandnum,dirnnum,:),[],1),style);
    end
    legendstr{curvenum}=['Band ' num2str(bandnums(bandnum)) ', ' dirnstr{dirnnum}];   
    curvenum=curvenum+1;
end
end
legends=propval(varargin,'legends');
if ~isempty(legends)
    legend(hplot,legends,'Location','best');
end

ylabel('Energy (eV)');
ylim(Elimits);
xlim([1 numkpoints]);
set(gca,'Box','on');

currylim=ylim;
for labelnum=1:length(symm_inds)
    line([symm_inds(labelnum); symm_inds(labelnum)],[currylim(1); currylim(2)],'LineWidth',2,'Color','k');
    text(symm_inds(labelnum),currylim(1)-(currylim(2)-currylim(1))/15,symm_labels(labelnum),'Clipping','off');
end
line(xlim,[0 0],'LineWidth',1,'Color','k');

if exchangeplot
    subplot(2,1,2); hold on;
    xlim([1 numkpoints]);
    for bandnum=1:length(bandnums);
        colchar=colchars(mod(bandnums(bandnum),7)+1);
        hplotxc(bandnum)=plot(1:numkpoints,reshape(E(bandnums(bandnum),2,:)-E(bandnums(bandnum),1,:),[],1),['-' colchar]);
    end
    legend(hplotxc,legendstr,'Location','best');
    ylim('auto');
    currylim=ylim; ylim([0 currylim(2)]);
    currylim=ylim;
    for labelnum=1:length(symm_inds)
        line([symm_inds(labelnum); symm_inds(labelnum)],[currylim(1); currylim(2)],'LineWidth',2,'Color','k');
        text(symm_inds(labelnum),currylim(1)-(currylim(2)-currylim(1))/3,symm_labels(labelnum),'HorizontalAlignment','center');
    end
    line(xlim,[0 0],'LineWidth',1,'Color','k');
    xlim([1 numkpoints]);
    ylabel('Energy (eV)');
end
title(['SC without s-o on 50000kpt mesh, FSM state M=0.265\mu_B f.u.^{-1} (Zr3009a)']);
set(gca,'Units','centimeters');
cpos=get(gca,'Position');
set(gca,'Position',[cpos(1) cpos(2)+1 16 6]);
set(gca,'XTickLabel','');

