function results=calcmoment(bandsdata,bandnums,FermiLevels)
%function results=calcmoment(bandsdata,bandnums,FermiLevels)
%
%calc the no. of carriers per unit cell of bands in bandsdata structure
%if bandsdata has spinpolarized bands (indicated by spindirns array), returns a
%moment; spindirns=1 +ve, 2 -ve
%uses IDOS(E) already in the structure
%if an array of bandnums is provided only entries in the array are used
%returns array of band indexes and band-moments with total moment in last
%column

if isempty(bandnums)
    bandnums=1:length(bandsdata.Wien2k_bandnums);
end
spindirns=bandsdata.spindirns(bandnums);
if isempty(FermiLevels)
    FermiLevels=bandsdata.FermiLevel(spindirns);
end

numbands=length(bandnums);
for bandnum=1:numbands
    if FermiLevels(bandnum)<min(bandsdata.DOS_E{bandnums(bandnum)})
        moment(bandnums(bandnum))=0;
        DOS(bandnums(bandnum))=0;
    else
        if FermiLevels(bandnum)>max(bandsdata.DOS_E{bandnums(bandnum)})
            moment(bandnums(bandnum))=1;
            DOS(bandnums(bandnum))=0;
        else
            moment(bandnums(bandnum))=interp1(bandsdata.DOS_E{bandnums(bandnum)},bandsdata.IDOS{bandnums(bandnum)},...
                FermiLevels(bandnum));
            DOS(bandnums(bandnum))=interp1(bandsdata.DOS_E{bandnums(bandnum)},bandsdata.DOS{bandnums(bandnum)},...
                FermiLevels(bandnum));
        end
    end
end
momentsign=(spindirns==0)+(spindirns==1)-(spindirns==2);
DOS=DOS/13.6/2; %converts per Ryd to per eV and, for ZrZn2 per unit cell to per f.u.

disp('band 1, band 2, ..., tot_up, tot_dn, total');
disp('moment, ...                        , moment');
disp('DOS(E_F), ...                      , DOS(E_F)');
results=[[bandnums -1 -1 -1]; ...
    [moment.*momentsign sum(moment(spindirns~=2).*momentsign(spindirns~=2)) sum(moment(spindirns==2).*momentsign(spindirns==2)) sum(moment.*momentsign)]; ...
    [DOS sum(DOS(spindirns~=2)) sum(DOS(spindirns==2)) sum(DOS) ]];
