function [Ebins DOS IDOS]=DOS(Estepsize,cartE,convE,convcrit,interpmethod,makeplot,silent)
%function [Ebins DOS IDOS]=DOS(Estepsize,cartE,convE,convcrit,interpmethod,makeplot,silent)
%
%calculates the DOS and integrated DOS evaluated at intervals of Estepsize
%randomly samples kspace, and builds up distr fn vs energy
%cartE contains 3d energy matrix
%convE is energy at which convergence criterion is tested
%convcrit is either [frac std error on DOS(convE) for convergence] or [fracerror abserror]
%interpmethod is 'spline' or 'linear'
%makeplot is true or false
%silent=true means no console output

switch interpmethod
    case 'spline'
        nsam=1e4;
    case 'linear'
        nsam=1e6;
end
converged=false;
convfrac=convcrit(1);
cumnsam=0;
Ebins=Estepsize*floor(min(cartE(1:end))/Estepsize+0.5):Estepsize:max(cartE(1:end));
hcum=zeros(size(Ebins));
convbin=find(hist(convE,Ebins));

dims=size(cartE);
while ~converged
    X=rand(nsam,1)*(dims(1)-1)+1;
    Y=rand(nsam,1)*(dims(2)-1)+1;
    Z=rand(nsam,1)*(dims(3)-1)+1;
    E=interp3(cartE,X,Y,Z,interpmethod);    
    h=hist(E,Ebins);
    hcum=hcum+h;
    cumnsam=cumnsam+nsam;
    %convergence criterion based on statistics at convE e.g. E_F  
    pop=hcum(convbin);
    converged=pop>(1/convfrac^2);
    if numel(convcrit)==2
        converged=converged || pop<(Estepsize*cumnsam*convcrit(2))^2;
    end
    if ~silent
        disp(['nsam=' num2str(cumnsam) '. Pop convbin=' num2str(pop) '. Conv if pop>' num2str((1/convfrac^2)) ' or pop<' num2str((Estepsize*cumnsam*convcrit(2))^2)]);
    end
end

DOS=hcum/Estepsize/cumnsam;
IDOS=cumsum(DOS);
if makeplot
    hfig=findobj('Tag','DOSPlot');
    if isempty(hfig)
        figure('Tag','DOSPlot');
    else
        figure(hfig);
    end
    plot(Ebins,DOS,'b');
    hold on;
    line([convE convE],[0,max(DOS)],'Color',[0 0 0]);
    xlabel('E (Ryd)');
    ylabel('DOS (states/uc/spin/Ryd)');
end