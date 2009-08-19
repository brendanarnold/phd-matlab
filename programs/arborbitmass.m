function arborbitmass(bandsdata,deltaE,bandnums,Bdirn,kpara,inclpoint,params);
%function arborbitmass(bandsdata,deltaE,bandnums,Bdirn,kpara,inclpoint,params);
%
%calculates mass of orbits specified by band,Bdirn,kpara,intpoint(point that orbit surrounds)
%does not have to be extremal orbit
Ryd2Joules=13.6*1.6e-19;

if isempty(params)
    load('defsliceparams.mat');
    disp(['Using default parameters: largest slice=' num2str(params.num_pts) 'x' num2str(params.num_pts) ...
        ', no. of slices per B dirn=' num2str(params.num_slices) ', min_F=' num2str(params.minfreq)]);
end

[inplanedirns minwidth]=getinplanedirn(Bdirn,bandsdata);
dk=minwidth/(params.num_pts-1);

inplane_inclpt.x=dot(inplanedirns{1},(inclpoint-kpara*Bdirn));
inplane_inclpt.y=dot(inplanedirns{2},(inclpoint-kpara*Bdirn));
%get orbit areas
allorbits{1}=[];
for bandnum=1:length(bandnums)
    slice=sliceFS(Bdirn,inplanedirns,kpara,dk,bandsdata,bandnums(bandnum));
    [orbits vertices]=orbitareas(slice,bandsdata,bandsdata.FermiLevel(bandsdata.spindirns(bandnums(bandnum))),false,true);    
    [orbits2 vertices2]=orbitareas(slice,bandsdata,deltaE+bandsdata.FermiLevel(bandsdata.spindirns(bandnums(bandnum))),false,true);        

    f1found=false; f2found=false;
    for orb=1:size(orbits{1},1)
        %convert vertices from indices to in-plane coords
        abs_vertices{orb}=[interp2(permute(slice.X,[2 1]),vertices{1}{orb}.x',vertices{1}{orb}.y') ...
            interp2(permute(slice.Y,[2 1]),vertices{1}{orb}.x',vertices{1}{orb}.y') ...
            interp2(permute(slice.Z,[2 1]),vertices{1}{orb}.x',vertices{1}{orb}.y')];
        for rownum=1:numel(vertices{1}{orb}.x)
            ip_vertices{orb}.x(rownum)=dot(inplanedirns{1},(abs_vertices{orb}(rownum,:)-kpara*Bdirn));
            ip_vertices{orb}.y(rownum)=dot(inplanedirns{2},(abs_vertices{orb}(rownum,:)-kpara*Bdirn));
        end
           
        abs_vertices2{orb}=[interp2(permute(slice.X,[2 1]),vertices2{1}{orb}.x',vertices2{1}{orb}.y') ...
            interp2(permute(slice.Y,[2 1]),vertices2{1}{orb}.x',vertices2{1}{orb}.y') ...
            interp2(permute(slice.Z,[2 1]),vertices2{1}{orb}.x',vertices2{1}{orb}.y')];
        for rownum=1:numel(vertices2{1}{orb}.x)
            ip_vertices2{orb}.x(rownum)=dot(inplanedirns{1},(abs_vertices2{orb}(rownum,:)-kpara*Bdirn));
            ip_vertices2{orb}.y(rownum)=dot(inplanedirns{2},(abs_vertices2{orb}(rownum,:)-kpara*Bdirn));
        end

        if inpolygon(inplane_inclpt.x,inplane_inclpt.y,ip_vertices{orb}.x,ip_vertices{orb}.y)
            if f1found==true
                disp('>1 matching orbit for f1');
            end
            f1=orbits{1}(orb,2);
            f1found=true;
        end
        if inpolygon(inplane_inclpt.x,inplane_inclpt.y,ip_vertices2{orb}.x,ip_vertices2{orb}.y)
            if f2found==true
                disp('>1 matching orbit for f2');
            end
            f2=orbits2{1}(orb,2);
            f2found=true;
        end
    end
    mass(bandnum)=1.6e-19*((6.626e-34)/2/pi)/(9.1e-31)*(f2-f1)/(deltaE*Ryd2Joules);
    outstr=sprintf('F1=%4.1f T, F2=%4.1f T, m=%1.3fm_e',f1,f2,mass(bandnum));
    disp(outstr);
end