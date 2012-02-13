function testBZ;
%test BZ detection function insideBZ.m by plotting BZ surface

rec_latt_vecs=[-1 1 1; 1 -1 1; 1 1 -1];
cart2reclatt=inv(rec_latt_vecs);

[kgridX kgridY kgridZ]=ndgrid(-4:0.2:4,-4:0.2:4,-4:0.2:4);
outsideBZ=zeros(size(kgridX));

for z=1:size(kgridZ,3)
    for y=1:size(kgridY,2)
        for x=1:size(kgridZ,1)
            outsideBZ(x,y,z)=~insideBZ([kgridX(x,y,z) kgridY(x,y,z) kgridZ(x,y,z)],2,rec_latt_vecs);
        end
    end
end

