function result=ucfacenormals(rl_vecs);
%returns the six faces of the conventional unit cell given 3x3 matrix rl_vecs:
%|a*_x a*_y a*_z|
%|b*_x b*_y b*_z|
%|c*_x c*_y c*_z|


a=rl_vecs(1,:); b=rl_vecs(2,:); c=rl_vecs(3,:);
acrossb=cross(rl_vecs(1,:),rl_vecs(2,:));
bcrossc=cross(rl_vecs(2,:),rl_vecs(3,:));
acrossc=cross(rl_vecs(1,:),rl_vecs(3,:));
result=[acrossb/norm(acrossb)*0.0001; acrossb/norm(acrossb)*norm(c)*dot(acrossb,c)/(norm(acrossb)*norm(c)); ...
    bcrossc/norm(bcrossc)*0.0001; bcrossc/norm(bcrossc)*norm(a)*dot(bcrossc,a)/(norm(bcrossc)*norm(a)); ...
    acrossc/norm(acrossc)*0.0001; acrossc/norm(acrossc)*norm(b)*dot(acrossc,b)/(norm(acrossc)*norm(b));];

    