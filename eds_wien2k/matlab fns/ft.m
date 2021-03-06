function [freq ftamp ftphase]=ft(x,y,padf,minx,maxx,polyord,maxf,doplot);
%function [freq ftamp ftphase]=ft(x,y,padf,minx,maxx,polyord,maxf);

if nargin==0
    disp(['Usage: freq ftamp ftphase]=ft(x,y,padf,minx,maxx,polyord,maxf)']);
    x=evalin('base','input(''X values (e.g. 1/field):'')');
    y=evalin('base','input(''Y values:'')');
    padf=evalin('base','input(''padding factor:'')');
end

if isempty(maxf)
    diffs=sort(abs(diff(x)));
    maxf=2*1/diffs(ceil(numel(diffs)*0.1));
end
delta_x=1/maxf/2;

[x order]=sort(x);
y=y(order);
inds=[find(x>=minx,1) find(x<=maxx,1,'last')];
inds(1)=max(1,inds(1)-1); inds(2)=min(inds(2)+1,numel(x));
x=x(inds(1):inds(2)); y=y(inds(1):inds(2));
%interp in x
int_x=(minx:delta_x:maxx)';
numpts=numel(int_x);
freq=2*linspace(0,maxf,numpts*(1+padf))';
window=sin((int_x-min(int_x))/(maxx-minx)*pi); %
%window=sin((int_x-min(int_x))/(maxx-minx)*pi).^2;
int_y=interp1(x,y,int_x);
if doplot
    subplot(2,2,1);
    plot(x,y,'-b',int_x,int_y,'-r');
    legend({'raw y vals','interp y vals'});
end

int_y=int_y.*window;

%subtract polynomial fit
if ~isempty(polyord)
    [pfit S mu]=polyfit(int_x,int_y,polyord);
    int_y=int_y-polyval(pfit,int_x,[],mu);
end

if doplot
    subplot(2,2,2);
    plot(int_x,int_y,'-b');
    ylabel('poly sub y');
end


pad_int_y=[int_y; repmat(zeros(size(int_y)),[padf 1])];
%calc FFT
ft=fft(pad_int_y);
ftamp=abs(ft);
ftphase=angle(ft);
half=ceil(numel(ft)/2);
ftamp=ftamp(1:half); ftphase=ftphase(1:half); freq=freq(1:half); ftcomplex=ft(1:half);

if doplot
    subplot(2,2,3);
    
    plot(freq,ftamp);
    xlim([0 3000]);
    ylabel('FT amp');
end

if doplot
    subplot(2,2,4); 
    phplot_Fs=linspace(0,maxf,numel(ftamp)*2);
    phplot_phs=linspace(-1,2,100);
    [phgrid_F phgrid_ph]=ndgrid(phplot_Fs,phplot_phs);
    %phase-shifted Fourier amp function (see PRL 93,166402)
    phgrid_K=real(interp1(freq,ftcomplex,phplot_Fs')*exp(i*pi*phplot_phs));
    phgrid_K(phgrid_K<0)=0;
    maxamp=max(phgrid_K(1:end));    
    levels=maxamp*[0 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.40 0.5 0.6 0.7 0.8 0.95 0.97 0.98 0.985 0.99 0.995 0.998 1];
    contourf(phgrid_F,phgrid_ph,phgrid_K,levels);
end   