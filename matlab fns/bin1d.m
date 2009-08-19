function [binx biny err]=bin1d(x,y,interval,varargin);
%[binx biny err]=function bin1d(x,y,interval,varargin);
%bins a 1d data set into bins that are interval wide

if numel(varargin)>=1
    binedge1=varargin{1}-interval/2;
else
    binedge1=min(x);
end

binedges=binedge1:interval:max(x);
binx=binedges(1:end-1)+interval/2;
for n=1:numel(binx)
    binys=y(x<binedges(n+1) & x>=binedges(n));
    if ~isempty(binys)
        err(n)=std(binys)/sqrt(numel(binys));
        biny(n)=mean(binys);
    else
        err(n)=nan;
        biny(n)=nan;
    end        
end

