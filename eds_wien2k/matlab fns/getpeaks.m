function [px py cf sgn]=getpeaks(x,y,k,varargin);
%function [px py cf sgn]=getpeaks(x,y,k,varargin);
%
%returns extrema of y(x) and curvature factors.
%criterion for maximum near y(n) is y(n-k)<=...<=y(n)>=...>=y(n+k)
%fits parabolas to points n-k ... n+k to get extrema locations and cf

if nargin>3
    mincf=varargin{1};
else
    mincf=[];
end

n=numel(y);
for j=[-k:k-1]
    sm(:,j+k+1)=diff(y(j+k+1:j+k+1+n-2*k));
end

pks=all([sm(:,1:k)>=0 sm(:,k+1:2*k)<=0],2);
trs=all([sm(:,1:k)<=0 sm(:,k+1:2*k)>=0],2);

extrema=find(pks|trs);
remove=zeros(size(extrema));
for extn=1:numel(extrema)
    sgn(extn)=1; 
    if trs(extrema(extn)) 
        sgn(extn)=-1;
    end
    extinds=extrema(extn)+(0:2*k);
    fit=polyfit(x(extinds),y(extinds),2);
    px(extn)=-fit(2)/2/fit(1);
    py(extn)=polyval(fit,px(extn));        
    cf(extn)=fit(2);
    if ~isempty(mincf) & cf(extn)<mincf
        remove(extn)=1;
    end
end    
px=px(~remove); py=py(~remove); cf=cf(~remove);
    



