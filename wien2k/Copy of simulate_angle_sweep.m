function [area]=simulate_angle_sweep(fs, angle, kz, FLshift, plotit, plot3d)

theta=angle*pi/180;
xf=2*pi/fs.latt_params(1);
yf=2*pi/fs.latt_params(2);
zf=2*pi/fs.latt_params(3);
%need to scale theta according to unit cell dimensions
%theta=atan(tan(theta)*fs.latt_params(1)/(fs.latt_params(3)));
limit=0.19*xf/cos(atan(tan(theta)*fs.latt_params(3)/(fs.latt_params(1))));
step=0.01;
[x y]=meshgrid(-limit:step:limit);
%fgrid=[x(1,:)' x(1,:)' zeros(size(x(1,:)',1),1)];
gsize=size(x,1);
fgrid=[reshape(x,gsize2,1) reshape(y,gsize2,1) zeros(gsize2,1)];

rotmat=[cos(theta) 0 sin(theta);0 1 0;-sin(theta) 0 cos(theta)];
rotgrid=(rotmat*fgrid')';
%made a 2D grid at angle theta now shift parallel to the desired level
rotgrid(:,3)=rotgrid(:,3)+zf*kz*ones(size(rotgrid,1),1);
%interpolate 3D data onto tilted 2D grid
egy=interp3(fs.cartX*xf,fs.cartY*yf,fs.cartZ*zf,fs.cartE,rotgrid(:,1),rotgrid(:,2),rotgrid(:,3));
%egy=interp3(fs.cartX,fs.cartY,fs.cartZ,fs.cartE,x,y,z,'spline');
%get energies back into required plaid format: using original x y grid
EI=griddata(fgrid(:,1),fgrid(:,2),egy,x,y);
%produce contour plot of this data
c=contour(x,y,EI,fs.FermiLevel+FLshift/13606.0);
close;
%hold
%c=contourc(egy,fs.FermiLevel);
%Assume just one contour, remove the index data at the start
cords=c(:,2:c(2,1));
% dump the cords for external plotting
%cordst=cords';
%save 'c:\transfer\contour.dat' '-ascii' 'cordst'
%plot the actual datapoints if desired
if plotit 
  plot(cords(1,:),cords(2,:),'.');
   set(gca,'DataAspectRatioMode','manual');
  set(gca,'DataAspectRatio',[1 1 1]);
  set(gcf,'PaperPositionMode','auto');
end;
%plot 3D
%calculate the area
format long;
area=10475.7*polyarea(cords(1,:),cords(2,:));
if plot3d
  plotcords=[cords(1,:);cords(2,:);zeros(1,size(cords,2))];
  plotcords=(rotmat*plotcords);
 
isosurface(fs.cartX*xf,fs.cartY*yf,fs.cartZ*zf,fs.cartE,fs.FermiLevel);
  hold;
  plot3(plotcords(1,:),plotcords(2,:),plotcords(3,:)+kz*zf,'-');
  set(gca,'DataAspectRatioMode','manual');
  set(gca,'DataAspectRatio',[1 1 1]);
  set(gcf,'PaperPositionMode','auto');
end;






