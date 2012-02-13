function tDOS=totalDOS(bandsdata,bandnums,FermiLevels);


if isempty(bandnums)
    bandnums=1:length(bandsdata.Wien2k_bandnums)
end
if isempty(FermiLevels)
    FermiLevels=bandsdata.FermiLevel(bandsdata.spindirns(bandnums));
end

numbands=length(bandnums);
for bandnum=1:numbands
    bDOS(bandnum)=interp1(bandsdata.DOS_E{bandnums(bandnum)},bandsdata.DOS{bandnums(bandnum)},FermiLevels(bandsdata.spindirns(bandnums(bandnum))));
end

tDOS(1)=sum(bDOS(bandsdata.spindirns==1));
tDOS(2)=sum(bDOS(bandsdata.spindirns==2));