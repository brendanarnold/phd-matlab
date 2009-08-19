function result=insideBZ(k,n,rec_latt_vecs)
%returns true if point k is inside the 1st to nth BZ belonging to rec_latt_vecs
%only single kpoint at present

numpts=3; %matrix of distances will have (2*numpts+1)^3 elements

%make matrix of distances from k to rec_latt_points that are close to gamma
[a b c]=ndgrid(-numpts:numpts,-numpts:numpts,-numpts:numpts);
rec_latt_points.X=a*rec_latt_vecs(1,1)+b*rec_latt_vecs(2,1)+c*rec_latt_vecs(3,1);
rec_latt_points.Y=a*rec_latt_vecs(1,2)+b*rec_latt_vecs(2,2)+c*rec_latt_vecs(3,2);
rec_latt_points.Z=a*rec_latt_vecs(1,3)+b*rec_latt_vecs(2,3)+c*rec_latt_vecs(3,3);
distances=sqrt((rec_latt_points.X-k(1)).^2+(rec_latt_points.Y-k(2)).^2+(rec_latt_points.Z-k(3)).^2);

%gamma point is nth nearest neighbour
nthneighbour=length(find(distances<=distances(numpts+1,numpts+1,numpts+1)));
if nthneighbour<=n
    result=true;
else
    result=false;
end
    

