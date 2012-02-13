function o=int3nearest(X,Y,Z,v,xyz);
%function o=int3nearest(X,Y,Z,v,xyz);
%
%finds v at closest point in X,Y,Z to xi yi zi, X Y Z can be any mesh

xi=repmat(xyz(1),size(X)); yi=repmat(xyz(2),size(X)); zi=repmat(xyz(3),size(X));
dsqrd=(X-xi).^2+(Y-yi).^2+(Z-zi).^2;
ind=dsqrd==min(dsqrd(1:end));
o=v(ind);
