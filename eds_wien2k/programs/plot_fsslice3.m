function plot_FSslice3(haxes,normal,kpara,vol,bandsdata,FermiLevelShift,bandnums,bandcolors,justFermicontour,showlabels,plot3d,draworbits,params,varargin);
%function plot_FSslice3(haxes,normal,kpara,vol,bandsdata,FermiLevelShift,bandnums,bandcolors,justFermicontour,showlabels,plot3d,draworbits,params,varargin);

%
%makes contour plot of FS in plane perp. to normal at a point kpara along
%the direction of normal
normal=normal/norm(normal);

if isempty(params)
    load('defsliceparams.mat');
    disp(['Using default parameters: largest slice=' num2str(params.num_pts) 'x' num2str(params.num_pts) ...
        ', no. of slices per B dirn=' num2str(params.num_slices) ', min_F=' num2str(params.minfreq)]);
end

if isempty(bandnums)
    bandnums=1:length(bandsdata.Wien2k_bandnums);
end

for bandnum=1:numel(bandnums)
    FermiLevels(bandnum)=bandsdata.FermiLevel(bandsdata.spindirns(bandnums(bandnum)));
end
if ~isempty(FermiLevelShift)
    FermiLevels=FermiLevels+FermiLevelShift;
end

%getinplanedirn finds dirn in plane perp to normal and passing through
%gamma for which searchvol polygon is narrowest
[inplanedirns minwidth]=getinplanedirn(normal,bandsdata);
dk=minwidth/(params.num_pts-1);

%get orbit areas 
allorbits{1}=[];
slice=sliceFS(normal,inplanedirns,kpara,dk,bandsdata,bandnums);

[orbits verts]=orbitareas(slice,bandsdata,FermiLevels,false,true);    

%discard orbits with f outside limits
for bandnum=1:length(bandnums)
    if ~isempty(orbits{bandnum})
        keepmask=orbits{bandnum}(:,2)>=params.minfreq;
        orbits{bandnum}=orbits{bandnum}(keepmask,:);
        verts{bandnum}=verts{bandnum}(keepmask);
    end
end

%if axes handle isn't provided see if can find previously raised window
if isempty(haxes)
    plothandle=findobj('Tag','ContourPlot');
    if isempty(plothandle) | ~ishandle(plothandle)
        plothandle=figure('Tag','ContourPlot');
    end
    figure(plothandle);
    haxes=gca;
else
    axes(haxes);
end
cla;
hold on;

orb_h=nan(length(bandnums),1);
colours=propval(varargin,'colours');
styles=propval(varargin,'styles');

for bandnum=1:length(bandnums)
    if ~isempty(colours)
        colour=colours(1+mod(bandnum-1,length(colours)));
    else
        colour=bandcolors(bandnum,:);
    end
    if ~isempty(styles)
        style=styles(1+mod(bandnum-1,length(styles)));
    else
        stylestrs={'-' '-.'};
        style=stylestrs{bandsdata.spindirns(bandnums(bandnum))};
    end
    %draw orbits
    if draworbits
    for orbnum=1:size(orbits{bandnum},1)
        inplaneverts=[interp2(permute(slice.inplaneX,[2 1]),verts{bandnum}{orbnum}.x',verts{bandnum}{orbnum}.y') ...
            interp2(permute(slice.inplaneY,[2 1]),verts{bandnum}{orbnum}.x',verts{bandnum}{orbnum}.y')];           
        orb_h(bandnum)=plot(inplaneverts(:,1),inplaneverts(:,2),'Color',colour,'LineStyle',style,'LineWidth',1.5);
    end
    end

    if plot3d
        wrappedcoords=wrapcoords(slice.X,slice.Y,slice.Z,bandsdata);        
        contourslice(haxes,bandsdata.cartX,bandsdata.cartY,bandsdata.cartZ,bandsdata.cartE{bandnums(bandnum)},wrappedcoords.X,wrappedcoords.Y,wrappedcoords.Z,[FermiLevels FermiLevels]);
    else        
        [cdata hcont(bandnum)]=contour(haxes,slice.inplaneX,slice.inplaneY,slice.E{bandnum},[FermiLevels(bandnum) FermiLevels(bandnum)],'LineColor',colour);
        set(hcont(bandnum),'LineStyle',style);
    end
end

kp=normal*kpara;
%draw bndry of volume included in sliceFS i.e. bandsdata.searchvol
bplanes=bandsdata.searchvol;
[bkps bnrms]=normdir(bplanes);
plist=boundingpoly(kpara,normal,bkps,bnrms,[0 0 0]);
%express vertices of bounding polygon in in-plane coords rel. to point kp
for pnum=1:size(plist,1)
    planeplist(pnum,1)=dot(plist(pnum,:)-kp,inplanedirns{1});
    planeplist(pnum,2)=dot(plist(pnum,:)-kp,inplanedirns{2});
end
plot([planeplist(:,1); planeplist(1,1)],[planeplist(:,2); planeplist(1,2)],'-k');

xlims=xlim;

%draw BZ boundary in plane
bplanes=bandsdata.BZfacenormals(:,:)*bandsdata.rec_latt_vecs;
[bkps bnrms]=normdir(bplanes);
plist=boundingpoly(kpara,normal,bkps,bnrms,[0 0 0]);

%draw 3d BZ (need to calc in-plane coords of patch vertices)
bplane_normals=bandsdata.BZfacenormals*bandsdata.rec_latt_vecs;
for n=1:size(bplane_normals,1)
    bkparas(n)=norm(bplane_normals(n,:));
    bplane_normals(n,:)=bplane_normals(n,:)/norm(bplane_normals(n,:));
end

[BZfaces BZvertices]=vol_vertices(bkparas,bplane_normals,[0 0 0]);
for vnum=1:size(BZvertices,1)
    pBZvertices(vnum,1)=dot(BZvertices(vnum,:)-kp,inplanedirns{1});
    pBZvertices(vnum,2)=dot(BZvertices(vnum,:)-kp,inplanedirns{2});
    pBZvertices(vnum,3)=dot(BZvertices(vnum,:)-kp,cross(inplanedirns{1},inplanedirns{2}));
end
C=[0 1 0];
h.BZpatch=patch('Faces',BZfaces,'Vertices',pBZvertices,'FaceColor',C,'FaceAlpha',0.05,'EdgeColor',[0 0 0]);

%express vertices of bounding polygon in in-plane coords rel. to point kp
planeplist=[];
for pnum=1:size(plist,1)
    planeplist(pnum,1)=dot(plist(pnum,:)-kp,inplanedirns{1});
    planeplist(pnum,2)=dot(plist(pnum,:)-kp,inplanedirns{2});
    if showlabels
        labelstr=sprintf('[%1.2f,%1.2f,%1.2f]',plist(pnum,1:3));
%        text(planeplist(pnum,1)+(xlims(2)-xlims(1))/50,planeplist(pnum,2),labelstr);
    end
end
plot([planeplist(:,1); planeplist(1,1)],[planeplist(:,2); planeplist(1,2)],'-ob','MarkerFaceColor','b','MarkerSize',2);
daspect(gca,[1 1 1]);
legends=propval(varargin,'legends');
if ~isempty(legends)
    legend(hcont(ishandle(hcont)),legends(ishandle(hcont)));
end
normalstr=sprintf('%4.2f,%4.2f,%4.2f',normal);
kparastr=sprintf('%5.3f',kpara);
title([bandsdata.pathcasename '. Slice perp [' normalstr '], k_para=' kparastr],'Interpreter','none');

%-------------------------------------------------------------------------
function color=getbandcolor(bandnum,spindirn,palette);
color=palette(1+mod(spindirn+2*bandnum-2,8),:);

