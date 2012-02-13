function [energies DOS IDOS vol]=DOS(Estepsize,cartE,interpiterations)
%function [energies DOS IDOS vol]=DOS(Estepsize,cartE,interpiterations)
%
%calculates the DOS and integrated DOS evaluated at Estepsize steps from
%Emin to Emax. vol is number of elements that cartE is divided into for calc

cartE=interp3(cartE,interpiterations);
cartE=cartE(1:end-1,1:end-1,1:end-1);
vol=numel(cartE);
Emax=max(reshape(cartE,[],1));
Emin=min(reshape(cartE,[],1));

energies=Emin:Estepsize:Emax;
numsteps=length(energies);

for stepnum=1:numsteps
    IDOS(stepnum)=length(find(cartE<energies(stepnum)));
end
IDOS=IDOS/vol;
DOS=[0 diff(IDOS)/Estepsize];