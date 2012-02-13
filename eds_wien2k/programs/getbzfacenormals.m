function result=getBZfacenormals(latt_type);
%function result=getBZfacenormals(latt_type);
%
%note latt_type is type of real-space lattice
if strcmp(latt_type,'F')
    %define normals of BZ bounding planes in rec latt coords
    result=[1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1; 1 1 1; -1 -1 -1; ... %nn
        1 1 0; 1 0 1; 0 1 1; -1 -1 0; -1 0 -1; 0 -1 -1]/2; %nnn    
    %i.e. these are the bisected rec_latt_vecs
elseif strcmp(latt_type,'C')
    result=[1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1; -1 1 0; 1 -1 0]/2;
elseif strcmp(latt_type,'B')
    result=[0 0 1; 0 0 -1; 1 -1 0; -1 1 0; 1 0 0; -1 0 0; 0 1 0; 0 -1 0; 1 0 -1; -1 0 1; 0 1 -1; 0 -1 1; 1 1 -1; -1 -1 1]/2;
elseif strcmp(latt_type,'H')
    result=[1 0 0; 0 1 0; -1 1 0; -1 0 0; 0 -1 0; 1 -1 0; 0 0 1; 0 0 -1]/2;    
elseif strcmp(latt_type,'P')
    result=[1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1]/2;
end