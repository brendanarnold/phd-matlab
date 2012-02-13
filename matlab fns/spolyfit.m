function out=spolyfit(x,y,n,dx,dosc,ord,xmin,xmax)
%function spolyfit(x,y,n,ord)
%
%performs sliding polynomial fit over windows of 2n+1 points

doplot=true;
keep=x>=xmin & x<=xmax;
x=x(keep); y=y(keep);

if isempty(dx)
    n1=n; n2=n;
end

for m=1:numel(x)
    if ~isempty(dosc)
        dx=x(m)^2/660*dosc;
    end
    if ~isempty(dx)
        n1=find(abs(x(m)-x(m:-1:1))>=dx/2,1)-1;
        if isempty(n1)
            n1=m-1;
        end
        n2=find(abs(x(m)-x(m:end))>=dx/2,1)-1;
        if isempty(n2)
            n2=numel(x)-m;
        end
    end
    if m-n1<1 | m+n2<1
        debug=true;
    end
    
    fit=polyfit(x(m-n1:m+n2),y(m-n1:m+n2),ord);
    out(m,1)=mean(x(m-n1:m+n2));
    out(m,2:ord+2)=fit;
    out(m,ord+3)=mean(y(m-n1:m+n2));
    out(m,ord+4)=x(m-n1);
    out(m,ord+5)=x(m+n2);
    if doplot
        clf;
        plot(x(m-n1:m+n2),y(m-n1:m+n2),'+');
        fvalxs=linspace(x(m-n1),x(m+n2),100);
        fvals=polyval(fit,fvalxs);
        hold on;
        plot(fvalxs,fvals,'-r');
        a=1;
    end
end