function orbits=contareas(conts,lims,dims);
%function contareas=contareas(conts,lims,dims);
%
%takes output from contour fn and calcs areas etc. of all closed orbits

dkx=(lims(1,2)-lims(1,1))/dims(1);
dky=(lims(2,2)-lims(2,1))/dims(2);

currcol=1;
contnum=1;
notatend=size(conts,2)>0; orbits=zeros(0,5);
while notatend
    nverts=conts(2,currcol);
    verts.y=conts(1,currcol+1:currcol+nverts); %contourc interprets cols as x, rows as y
    verts.x=conts(2,currcol+1:currcol+nverts);
    if ~any(isnan(verts.x))      
        if any(abs(verts.y-lims(2,1))<dky/10 | abs(verts.y-lims(2,2))<dky/10 | abs(verts.x-lims(1,1))<dkx/10 | abs(verts.x-lims(2,1))<dkx/10)
           %find if contour is open in first BZ
            orbits(contnum,5)=1;        
        else
            orbits(contnum,5)=0;
            %contour is closed
        end
        %if contour open need to add a point at BZ corner to get area
        if orbits(contnum,5)==0
            orbits(contnum,1)=polyarea(verts.x,verts.y); %in Angstrom^-2
            orbits(contnum,2)=conts(1,currcol); %energy meV
            orbits(contnum,3)=mean(verts.x);
            orbits(contnum,4)=mean(verts.y);
            cverts{contnum}=verts;
        end
        contnum=contnum+1;
    end
    if currcol+nverts+1>size(conts,2)
        notatend=false;
    else
        currcol=currcol+nverts+1;
    end
end