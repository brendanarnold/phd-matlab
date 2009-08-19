function output=derivbyfit(X,Y,xi,xf,n,filename);
%function derivbyfit(X,Y,xi,xf,n);
%
%fits polynomial to data in X,Y in ranges xi->xf to find nth deriv

X=X(1:min([numel(X) numel(Y)]));
Y=Y(1:min([numel(X) numel(Y)]));
for pt=1:numel(xi)
    rX=X(X>=xi(pt) & X<xf(pt));
    rY=Y(X>=xi(pt) & X<xf(pt));
    
    meanX(pt)=mean(rX);
    p=polyfit(rX,rY,n);
    fity(pt)=polyval(p,meanX(pt));
    nderiv(pt)=p(end-n);
end
output(:,1)=meanX;
output(:,2)=fity;
output(:,3)=nderiv';

if ~isempty(filename);
    outfile=fopen(filename,'w');
    fprintf(outfile,'X\tY\tnderiv\r\n');
    for row=1:size(output,1)
        fprintf(outfile,'%6e\t%6e\t%6e\r\n',output(row,:));
    end
end
fclose(outfile);