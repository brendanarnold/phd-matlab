function [of,ov,ofcolor]=clippatch(f,v,bplane_normals,incl_point,clipintersectedfaces);

incl=ones(size(v,1),1);
for n=1:size(bplane_normals,1)
    ongammaside=((v(:,1)-bplane_normals(n,1))*bplane_normals(n,1)+...
        (v(:,2)-bplane_normals(n,2))*bplane_normals(n,2)+...
        (v(:,3)-bplane_normals(n,3))*bplane_normals(n,3))<=0;
    incl_gamma=(incl_point(1)-bplane_normals(n,1))*bplane_normals(n,1)+...
        (incl_point(2)-bplane_normals(n,2))*bplane_normals(n,2)+...
        (incl_point(3)-bplane_normals(n,3))*bplane_normals(n,3)<=0;
    if incl_gamma
        incl=incl & ongammaside;
    else
        incl=incl & (~ongammaside);
    end
end

disp([num2str(size(f,1)) ' faces, ' num2str(size(v,1)) ' vertices total.']);
%remove faces with all verts outside volume
outside=all(~incl(f),2);
inside=all(incl(f),2);
of=f(~outside,:);
ofcolor=repmat([1 0 0],size(of,1),1);
disp([num2str(size(of,1)) ' faces left after removing those wholly outside region']);

if clipintersectedfaces
inside=inside(~outside);

intersectedfaces=find(~inside);
disp([num2str(length(intersectedfaces)) ' faces are intersected by bounding planes']);

newvs=[];
for n=1:length(intersectedfaces);  
    vnew=[];
    for plnum=1:size(bplane_normals,1)        
        for vnum=0:2
            v1=v(of(intersectedfaces(n),1+mod(vnum,3)),:);
            v2=v(of(intersectedfaces(n),1+mod(vnum+1,3)),:);        
            denom=dot(v2-v1,bplane_normals(plnum,:));
            if denom~=0
                lamda=dot(bplane_normals(plnum,:)-v1,bplane_normals(plnum,:))/denom;
            else
                lamda=Inf;
            end
            if lamda>=0 & lamda<1
                %when some verts are on plane need to make sure don't
                %double count new vertices hence >=, <
                vnew(end+1,:)=v1+lamda*(v2-v1);                
                %disp(['Vertices ' num2str(1+mod(vnum,3)) ' and '
                %num2str(1+mod(vnum+1,3)) ' intersect with bplane ' num2str(plnum) ' at lamda=' num2str(lamda)]);                
            end
        end
        if size(vnew,1)==2
            %original face may be intersected by >1 plane but would need to iterate
            break;
        end
    end
    if size(vnew,1)~=2
        disp('Something wrong!');
    end
    newvs=[newvs; vnew];
end

orig_numverts=size(v,1);
v=[v; newvs];
for n=length(intersectedfaces):-1:1
    %which vertices are inside region?
    v_inside=find(incl(of(intersectedfaces(n),:)));    
        if length(v_inside)==1
        %if exactly one vertex is inside, just replace old face with a new
        %one containing both new vertices.
        of(intersectedfaces(n),:)=[of(intersectedfaces(n),v_inside) orig_numverts+(n-1)*2+2-(v_inside==1) orig_numverts+(n-1)*2+1+(v_inside==1)];
        ofcolor(intersectedfaces(n),:)=[0 1 1];
    elseif length(v_inside)==2        
        quadverts=[orig_numverts+(n-1)*2+1; orig_numverts+(n-1)*2+2; ...
            of(intersectedfaces(n),v_inside(1)); of(intersectedfaces(n),v_inside(2))];
        [uniq_verts ind1]=unique(v(quadverts,:),'rows');
        if length(uniq_verts)==3
            alteredface=[of(intersectedfaces(n),v_inside(1)) of(intersectedfaces(n),v_inside(2)) ...
                setdiff(quadverts(ind1),[of(intersectedfaces(n),v_inside(1)) of(intersectedfaces(n),v_inside(2))])];
            of(intersectedfaces(n),:)=alteredface;
            ofcolor(intersectedfaces(n),:)=[0 0 1];
        elseif length(uniq_verts)==4               
            p=setdiff([1 2 3],v_inside);
            alteredface=[of(intersectedfaces(n),1+mod(p,3)) of(intersectedfaces(n),1+mod(p+1,3)) orig_numverts+(n-1)*2+1+(p==1)];
            newface=[alteredface(3) orig_numverts+(n-1)*2+2-(p==1) alteredface(1)];
            of=[of(1:intersectedfaces(n)-1,:); alteredface; newface; ...
                of(intersectedfaces(n)+1:end,:)];            
            ofcolor=[ofcolor(1:intersectedfaces(n),:); [0 1 0]; ofcolor(intersectedfaces(n)+1:end,:)];
        end
    else
        disp('Something wrong!');
    end    
end

is_redundant=~ismember([1:length(v)]',of);
disp(['Removing ' num2str(length(find(is_redundant))) ' vertices']);
delta_index=cumsum(is_redundant);
of=of-delta_index(of);
v=v(~is_redundant,:);
disp('...done');
end
ov=v;

