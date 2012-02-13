function [px py cf sgn]=getpeaks2(x,y,Dx,dx,varargin);
%function [px py cf sgn]=getpeaks2(x,y,Dx,dx,varargin);
%
%returns extrema of y(x) and curvature factors.
%different method from getpeaks: tries running O(x^2) fits
%Dx is size of x interval to make fits over, dx is amount advanced between successive fits
if nargin>4
    mincf=varargin{1};
else
    mincf=[];
end

ranges=[min(x):dx:max(x)-Dx; min(x)+Dx:dx:max(x)];
nr=size(ranges,2);

for r=1:nr
    rinds=x>=ranges(1,r) & x<=ranges(2,r);
    fit{r}=polyfit(x(rinds),y(rinds),2);
    coeffs(r,:)=fit{r};
end

a=1;