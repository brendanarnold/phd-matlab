function [summary_results, vertices_results]=orbitareas(slice,bandsdata,E,plotslice,return_vertices);
%uses contourc to find orbits at energy E in slice and calculates areas
%result is a matrix with columns:
%1: area [Angstrom^-1]
%2: freq [T]
%3: centre kx [Angstrom^-1]
%4: centre ky [Angstrom^-1]
%5: centre kz [Angstrom^-1]
%note return of vertices not currently handled correctly for multiple bands

hbar=6.626e-34/(2*pi);
el_charge=1.6e-19;

if plotslice
    figure('Name',['Kpara=' num2str(slice.kpara)]);
    tE=slice.E{1};
    tE(slice.outsidebounds)=nan;
    contour(tE');
end
vertices_results=[];

orbits=cell(1,length(slice.E));
for bandnum=1:length(slice.E)
contours=contourc(slice.E{bandnum},[E(bandnum) E(bandnum)]);

currcol=1;
orbitnum=1;
notatend=size(contours,2)>0; orbits{1}=zeros(0,5);
while notatend
    num_vertices=contours(2,currcol);
    vertices.y=contours(1,currcol+1:currcol+num_vertices); %contourc interprets cols as x, rows as y
    vertices.x=contours(2,currcol+1:currcol+num_vertices);
    if any(interp2(permute(slice.outsidebounds,[2 1 3]),vertices.x,vertices.y)) ...
            | any(vertices.x==1 | vertices.x==size(slice.E{bandnum},1) | vertices.y==1 | vertices.y==size(slice.E{bandnum},2))
        %if any contour vertices are outsidebuonds OR if they touch matrix border, ignore this contour
    else
        %contour is complete
        orbits{bandnum}(orbitnum,1)=slice.inplanedx*slice.inplanedy*polyarea(vertices.x,vertices.y); %in Angstrom^-2
        orbits{bandnum}(orbitnum,2)=orbits{bandnum}(orbitnum,1)*(1e20)*(hbar/2/pi/el_charge);
        orbits{bandnum}(orbitnum,3)=interp2(permute(slice.X,[2 1 3]),mean(vertices.x),mean(vertices.y));
        orbits{bandnum}(orbitnum,4)=interp2(permute(slice.Y,[2 1 3]),mean(vertices.x),mean(vertices.y));
        orbits{bandnum}(orbitnum,5)=interp2(permute(slice.Z,[2 1 3]),mean(vertices.x),mean(vertices.y));    
        if plotslice
            hold on;
            plot(vertices.x,vertices.y,'-k','LineWidth',5);
            labelstr=sprintf('f=%2.2f, k_{centre}=%1.2f,%1.2f,%1.2f',orbits{bandnum}(orbitnum,2),...
                orbits{bandnum}(orbitnum,3),orbits{bandnum}(orbitnum,4),orbits{bandnum}(orbitnum,5));
            text(mean(vertices.x),mean(vertices.y),labelstr);
        end
        if return_vertices
            vertices_results{bandnum}{orbitnum}=vertices;
        end
        orbitnum=orbitnum+1;
    end
    if currcol+num_vertices+1>size(contours,2)
        notatend=false;
    else
        currcol=currcol+num_vertices+1;
    end
end
end

summary_results=orbits;


