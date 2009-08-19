function rawsquid(infile,scanpts,scansperpt);

fitrt='mf_flsqr';
fitfn='squid_Vofz';

info=file2cell([infile '.dc.raw'],30,true,',');

pos_col=find(strcmp(info.titles,'Position (cm)'));
sig_col=find(strcmp(info.titles,'Long Voltage'));


if mod(size(info.data,1),scanpts)~=0
    disp('Unexpected file length');
else 
    numscans=size(info.data,1)/scanpts;
end

scannum=1+floor((1:size(info.data,1))/scanpts);

fig=findobj('Tag','squidfig');
if ishandle(fig);
    figure(fig);
else
    figure('Tag','squidfig','Name','SQUID data fitter');
end
cla;

for n=1:numscans
    xs=info.data(scannum==n,pos_col);
    ys=info.data(scannum==n,sig_col);
    errs=ones(size(xs));
    fixed=[0 0 0 0];
    
    %intial guesses
    p(1)=1; p(2)=-4; p(3)=0; p(4)=1;
    [y, name, pnames, pin]=feval(fitfn,[],p,1);
    [p,dp]=feval(fitrt,xs,ys,errs,p,~fixed,fitfn);
    
    fitxs=linspace(min(xs)-(max(xs)-min(xs))/2,max(xs)+(max(xs)-min(xs))/2,100);
    fitvals=feval(fitfn,fitxs,p);
    plot(fitxs,fitvals);
    hold on;
    plot(xs,ys,'o');
end