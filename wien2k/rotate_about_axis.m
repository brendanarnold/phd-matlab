function [xf, yf, zf] = rotate_about_axis(xi, yi, zi, axis, angle)
% Rotates cartesian co-ordiantes in <xi>, <yi> and <zi> by <angle> degrees about
% <axis>

% Make axis a unit vector
u = axis ./ norm(axis);
ux = u(1);
uy = u(2);
uz = u(3);

a = angle * pi / 180;

% From http://en.wikipedia.org/wiki/Rotation_matrix
R = [cos(a)+ux*ux*(1-cos(a)) ux*uy*(1-cos(a))-uz*sin(a) ux*uz*(1-cos(a))+uy*sin(a); ...
    uy*ux*(1-cos(a))+uz*sin(a) cos(a)+uy*uy*(1-cos(a)) uy*uz*(1-cos(a))-ux*sin(a); ...
    uz*ux*(1-cos(a))-uy*sin(a) uz*uy*(1-cos(a))+ux*sin(a) cos(a)+uz*uz*(1-cos(a))];

% Apply rotation
result = [xi yi zi] * R;

% Divvy the results back up again
xf = result(:,1);
yf = result(:,2);
zf = result(:,3);

end