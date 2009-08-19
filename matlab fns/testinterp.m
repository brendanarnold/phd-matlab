function testinterp;

[x y z]=meshgrid([-1 1],[-1 1],[-1 1]);
E=ones(size(x));
E(:,:,1)=[1 2; 3 4];
E(:,:,2)=[5 6; 7 25];

[fx fy fz]=meshgrid(-1:0.1:1,-1:0.1:1,-1:0.1:1);
fE=interp3(x,y,z,E,fx,fy,fz);
[grdx grdy grdz]=gradient(fE);
[cx cy cz cfE]=reducevolume(fx,fy,fz,fE,[3 3 3]);

coneplot(fx,fy,fz,grdx,grdy,grdz,cx,cy,cz);
a=1;
