function result=DebyeInt(xmax);

intfn=@(x) x.^4.*exp(x)./(exp(x)-1).^2;
result=quad(intfn,0,xmax);
