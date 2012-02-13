function [NX,NY,NZ]=rotate3grid(X,Y,Z,vec,sense);
%returns coordinates NX,NY,NZ
%can view transformation either as rotating the input points in which case:
%sense=1 brings a point parallel to vec parallel to the z-axis
%sense=-1 brings a point parallel to the z-axis parallel to vec

vec=vec/sqrt(dot(vec,vec));
vec_cross_y=cross(vec,[0 1 0]);
if sqrt(dot(vec_cross_y,vec_cross_y))~=0
    vec_cross_y_hat=vec_cross_y/sqrt(dot(vec_cross_y,vec_cross_y));
end

%calculate rotation angle theta about Y axis to bring ZY coplanar with vec
%and phi, subsequent rotation angle about X' to bring Z' || vec

if (dot(vec_cross_y,vec_cross_y)~=0)
    costheta=dot(cross([0 1 0],vec_cross_y_hat),[0 0 1]);
    sintheta=dot(vec_cross_y_hat,[0 0 1]);
    cosphi=dot(cross([0 1 0],vec_cross_y_hat),vec);
    sinphi=-dot([0 1 0],vec);
else
    %no point rotating about Y if vec || or anti-|| to Y
    costheta=1;
    sintheta=0;
    
    %need to find whether vec is || or anti-|| to Y to get sign of phi
    cosphi=0;
    sinphi=-dot([0 1 0],vec);
end

RotMatrix=[costheta 0 -sintheta; 
    sinphi*sintheta cosphi sinphi*costheta;
    cosphi*sintheta -sinphi cosphi*costheta];

if (sense==-1)
    RotMatrix=inv(RotMatrix);
end

NX=RotMatrix(1,1)*X+RotMatrix(1,2)*Y+RotMatrix(1,3)*Z;
NY=RotMatrix(2,1)*X+RotMatrix(2,2)*Y+RotMatrix(2,3)*Z;
NZ=RotMatrix(3,1)*X+RotMatrix(3,2)*Y+RotMatrix(3,3)*Z;


