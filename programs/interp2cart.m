function [cartE dx dy dz cuboidextents]=interp2Cart(E,rec_latt_vecs,rec_latt_type);
%
%
%interpolates energies on a rec latt based grid (e.g. from BXSF or Wien2k
%case.energy) onto a Cartesian grid

%decide how to interp to Cartesian grid according to rec_latt_type

if ~iscell(E)
    temp{1}=E;
    E=temp;
end
dims=size(E{1});
for bandnum=1:length(E)
    if (rec_latt_type=='bcc')
            num_bcc_pts=dims(1);
            num_sc_pts=num_bcc_pts*2-1;
            cartE{bandnum}=zeros(num_sc_pts,num_sc_pts,num_sc_pts); 
            for Z_ind=1:num_sc_pts
            for Y_ind=1:num_sc_pts
            for X_ind=1:num_sc_pts                   
                a_ind=(Y_ind+Z_ind)/2;
                b_ind=(X_ind+Z_ind)/2;
                c_ind=(Y_ind+X_ind)/2;
                a_ind=1+mod(a_ind-1,num_bcc_pts-1);
                %wrap point back to first rec latt unit cell if necessary
                %remember input bcc matrix includes all points shared between cells
                %hence num_bcc_pts-1 for correct wrapping
                b_ind=1+mod(b_ind-1,num_bcc_pts-1);
                c_ind=1+mod(c_ind-1,num_bcc_pts-1);
                if ([a_ind b_ind c_ind]==round([a_ind b_ind c_ind]))
                    %transfer points that lie on original bcc lattice directly
                    cartE{bandnum}(X_ind,Y_ind,Z_ind)=E{bandnum}(a_ind,b_ind,c_ind);
                else
                    %points in denser s.c. lattice that do not lie on original bcc
                    %lattice have non-integers values of 2 of (a_ind, b_ind, c_ind)
                    %The 2 nearest neighbours are at rounded values (up,down)
                    a1_ind=1+mod(floor(a_ind)-1,num_bcc_pts-1);
                    b1_ind=1+mod(floor(b_ind)-1,num_bcc_pts-1);
                    c1_ind=1+mod(floor(c_ind)-1,num_bcc_pts-1);
                    a2_ind=1+mod(ceil(a_ind)-1,num_bcc_pts-1);
                    b2_ind=1+mod(ceil(b_ind)-1,num_bcc_pts-1);
                    c2_ind=1+mod(ceil(c_ind)-1,num_bcc_pts-1);
                    nnE=0.5*(E{bandnum}(a1_ind,b1_ind,c1_ind)+E{bandnum}(a2_ind,b2_ind,c2_ind));
                    %4 nextnearest neighbours:
                    if a_ind==floor(a_ind)
                        nnninds=[a_ind+1 b_ind+0.5 c_ind+0.5; ...
                            a_ind-1 b_ind-0.5 c_ind-0.5; ...
                            a_ind b_ind+0.5 c_ind-0.5; ...
                            a_ind b_ind-0.5 c_ind+0.5];
                    elseif b_ind==floor(b_ind)
                        nnninds=[a_ind+0.5 b_ind+1 c_ind+0.5; ...
                            a_ind-0.5,b_ind-1,c_ind-0.5; ...
                            a_ind+0.5 b_ind c_ind-0.5; ...
                            a_ind-0.5 b_ind c_ind+0.5];
                    elseif c_ind==floor(c_ind)
                        nnninds=[a_ind+0.5,b_ind+0.5,c_ind+1; ...
                            a_ind-0.5,b_ind-0.5,c_ind-1; ...
                            a_ind+0.5,b_ind-0.5,c_ind; ...
                            a_ind-0.5,b_ind+0.5,c_ind];
                    end
                    nnninds=1+mod(nnninds-1,num_bcc_pts-1);
                    nnnE=0.25*(E{bandnum}(nnninds(1,1),nnninds(1,2),nnninds(1,3))+...
                        E{bandnum}(nnninds(2,1),nnninds(2,2),nnninds(2,3))+...
                        E{bandnum}(nnninds(3,1),nnninds(3,2),nnninds(3,3))+...
                        E{bandnum}(nnninds(4,1),nnninds(4,2),nnninds(4,3)));
                    cartE{bandnum}(X_ind,Y_ind,Z_ind)=(2*nnE+nnnE)/3;                    
                end 
            end
            end    
            end
            dx=abs(rec_latt_vecs(1,1)/(dims(1)-1));
            dy=dx; dz=dx;
            cuboidextents=[0 0 0; dx*(num_sc_pts-1) dy*(num_sc_pts-1) dz*(num_sc_pts-1)];
            cartdims=[num_sc_pts,num_sc_pts,num_sc_pts];
        elseif (rec_latt_type=='hex')
            orthodims=[dims(1)*2-1 2*dims(2)-1 dims(3)];    
            cartE{bandnum}=zeros(orthodims(1),orthodims(2),orthodims(3)); 
            for Z_ind=1:orthodims(3)
            for Y_ind=1:orthodims(2)
            for X_ind=1:orthodims(1)                   
                a_ind=1+(X_ind-Y_ind)/2;
                b_ind=Y_ind;
                c_ind=Z_ind;
                a_ind=1+mod(a_ind-1,dims(1)-1);
                %wrap point back to first rec latt unit cell if necessary
                %remember input matrix from bxsf includes all points shared between cells
                %hence dims-1 for correct wrapping
                b_ind=1+mod(b_ind-1,dims(2)-1);
                c_ind=1+mod(c_ind-1,dims(3)-1);
                if ([a_ind b_ind c_ind]==round([a_ind b_ind c_ind]))
                    %transfer points that lie on original hex lattice directly
                    cartE{bandnum}(X_ind,Y_ind,Z_ind)=E{bandnum}(a_ind,b_ind,c_ind);
                else
                    %points in denser orthorhombic lattice that do not lie on original
                    %hex lattice have non-integers values of a_ind
                    %The 2 nearest neighbours are at rounded values (up,down)
                    a1_ind=1+mod(floor(a_ind)-1,dims(1)-1);
                    a2_ind=1+mod(ceil(a_ind)-1,dims(1)-1);
        
                    cartE{bandnum}(X_ind,Y_ind,Z_ind)=0.5*(E{bandnum}(a1_ind,b_ind,c_ind)+E{bandnum}(a2_ind,b_ind,c_ind));
                end 
            end
            end    
            end
            dx=abs(rec_latt_vecs(1,1)/(dims(1)-1))/2;
            dy=abs(rec_latt_vecs(2,2)/(dims(2)-1));
            dz=abs(rec_latt_vecs(3,3)/(dims(3)-1));
            cuboidextents=[0 0 0; dx*(orthodims(1)-1) dy*(orthodims(2)-1) dz*(orthodims(3)-1)];
    end    
end
