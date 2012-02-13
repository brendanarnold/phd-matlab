function DOStable=DOStable(bandsdata,Emin,Emax,deltaE,bandnums,pathfilename,makeplot);
%function DOStable=DOStable(bandsdata,Emin,Emax,deltaE,bandnums,pathfilename,makeplot);

dirnstr=['up'; 'dn'];
if isempty(bandnums)
    bandnums=1:length(bandsdata.Wien2k_bandnums);
end
spindirns=bandsdata.spindirns(bandnums);
FermiLevels=bandsdata.FermiLevel(spindirns);
if isempty(Emin)
    Emin=min(cell2mat(bandsdata.DOS_E));
end
if  isempty(Emax)    
    Emax=max(cell2mat(bandsdata.DOS_E));
end

Es=Emin:deltaE:Emax;
numbands=length(bandnums);
for bandnum=1:numbands
    DOStable(:,bandnum)=interp1(bandsdata.DOS_E{bandnums(bandnum)},bandsdata.DOS{bandnums(bandnum)},Es);
end
DOStable(:,numbands+1)=sum(DOStable(:,spindirns==1 | spindirns==0),2);
DOStable(:,numbands+2)=sum(DOStable(:,spindirns==2 | spindirns==0),2);
DOStable(:,numbands+3)=sum(DOStable(:,1:numbands),2);
DOStable=DOStable/13.6/2; %converts per Ryd to per eV and, for ZrZn2 per unit cell to per f.u.
Es=Es*13.6; FermiLevels=FermiLevels*13.6;
DOStable(isnan(DOStable))=0;
  
if ~isempty(pathfilename)
    outfile=fopen([pathfilename '.DOS.txt'],'w');
    fprintf(outfile,['DOS of ' bandsdata.pathcasename '. Units are states/eV/unit cell/2.\r\n']);
    fprintf(outfile,'Fermi level: %1.4f %1.4f\r\n',bandsdata.FermiLevel);
    coltitle{1}='E abs (eV)';
    coltitle{2}='E-E_F^up (eV)';
    coltitle{3}='E-E_F^dn (eV)';
    for c=1:numbands
        coltitle{c+3}=['band ' num2str(bandsdata.Wien2k_bandnums(c)) dirnstr(spindirns(c))];
    end
    coltitle{numbands+4}='total up';
    coltitle{numbands+5}='total dn';
    coltitle{numbands+6}='total';
    for c=1:numbands+6
        fprintf(outfile,'%s',coltitle{c});
        if c==(numbands+6)
            fprintf(outfile,'\r\n');
        else
            fprintf(outfile,',');
        end
    end
    for r=1:length(Es)
        linestr=sprintf('%1.3f,',DOStable(r,:));
        fprintf(outfile,'%1.4f,%1.4f,%1.4f,%s\r\n',Es(r),Es(r)-FermiLevels(1),Es(r)-FermiLevels(2),linestr(1:end-1));
    end
    fclose(outfile);
end
