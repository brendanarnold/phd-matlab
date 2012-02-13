function [NX, NY, NZ]=rotate3axis(X,Y,Z,axis,angle);
%rotates in 3D about arbitrary axis
%note there may be significant cumulative error ~1% after rotating 100times
%without normalizing lines at end

%rotate so that axis vector is parallel to z
[NX NY NZ]=rotate3grid(X,Y,Z,axis,+1);

%rotate about axis (now z)
NNX=cos(angle)*NX-sin(angle)*NY;
NNY=sin(angle)*NX+cos(angle)*NY;

%reverse original rotation
[NX NY NZ]=rotate3grid(NNX,NNY,NZ,axis,-1);

%make sure length doesn't change
%for elnum=1:numel(NX)
%    scalefactor=sqrt((NX(elnum)^2+NY(elnum)^2+NZ(elnum)^2)/(X(elnum)^2+Y(elnum)^2+Z(elnum)^2));
%    NX(elnum)=NX(elnum)/scalefactor;
%    NY(elnum)=NY(elnum)/scalefactor;
%    NZ(elnum)=NZ(elnum)/scalefactor;
%end