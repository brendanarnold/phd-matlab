function tempplotlines(r0,d,points,kp);

figure;
hold on;
lamda=-5:0.1:5;
for planenum=1:size(r0,1)
    xs=r0(planenum,1)+lamda*d(planenum,1);
    ys=r0(planenum,2)+lamda*d(planenum,2);
    zs=r0(planenum,3)+lamda*d(planenum,3);
    plot3(xs,ys,zs);
end
for planenum1=1:size(r0,1)
    for planenum2=planenum1+1:size(r0,1);
        plot3(points(planenum1,planenum2,1),points(planenum1,planenum2,2),points(planenum1,planenum2,3),'r+');
    end
end
plot3(kp(1),kp(2),kp(3),'g+');