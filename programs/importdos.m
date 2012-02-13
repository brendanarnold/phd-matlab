function data=importDOS(pathfilename,isspinpol,makeplot);
%function data=importDOS(pathfilename,isspinpol,makeplot);
%
%takes Wien2k DOS output file and converts to matlab matrix (cell array with 2 elements if spinpol)

Palette=[1 0 0; 0 0 1; 0 1 0; 0 0.7 0.7; 0.7 0.7 0; 1 0 1; 0 0 0; 0 0.5 0];
if isspinpol
    dirnstr={'up'; 'dn'};
else
    dirnstr={''};
end

for spindirn=1:1+isspinpol
    infile=fopen([pathfilename dirnstr{spindirn}],'r');

    linestr1=fgetl(infile);
    linestr2=fgetl(infile);
    linestr3=fgetl(infile);
    titles{spindirn}=delimstr2cell(linestr3,' ');
    titles{spindirn}=titles{spindirn}(2:end-1);
    notatend=1; row=1;
    while notatend
        linestr=fgetl(infile);
        if linestr==-1
            notatend=false;
        else
            linevalues=delimstr2cell(linestr,' ');
            data{spindirn}(row,:)=str2num(char(linevalues))';            
            row=row+1;
        end
    end
    data{spindirn}(:,2:end)=data{spindirn}(:,2:end)/2; %unit cell^-1 -> f.u.^-1
    fclose(infile);
end

if any(cell2mat(makeplot))
    curvenum=1;
    for spindirn=1:1+isspinpol
        selcols=find(makeplot{spindirn});
        for n=1:length(selcols)            
            plot(data{spindirn}(:,1),data{spindirn}(:,1+selcols(n)),'Color',Palette(1+mod(curvenum-1,size(Palette,1)),:));
            hold on;
            legendstr{curvenum}=[dirnstr{spindirn} ': ' titles{spindirn}{1+selcols(n)}];
            curvenum=curvenum+1;
        end        
    end
    legend(legendstr,'Location','best');
    ylimits=ylim;
    line([0 0],[ylimits(1) ylimits(2)],'Color',[0 0 0]);
end