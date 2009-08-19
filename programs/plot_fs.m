function [haxes,hFSpatches]=plot_FS(filenameORvar,bands,EF,bandcolours,bchar,vol,incl_point,centrek,undersample_factor,interpmethod,haxes,clipintersectedfaces);
%function [haxes,hFSpatches]=plot_FS(filenameORvar,bands,EF,bandcolours,bchar,vol,incl_point,centrek,undersample_factor,interpmethod,haxes,clipintersectedfaces);
%
%plots Fermi surface arising from 'bands' in the region 
%including incl_point that is bounded by bplane_normals (which are
%expressed in cartesian coords). by default, plots first BZ zone

%arguments
%   filenameorvar:          structure containing bandstructure info as output by convert_W2kE;
%   bands:                  numbers of bands to plot (empty means all)
%   EF:                     Fermi level (empty means use bands.FermiLevel)
%   bandcolours:            nx3 matrix of rgb colours
%   bchar:                  structure with fields x,y,z, and cell array of pchg matrices for bands
%                           x,y,z,pchg must cover the tileable cuboid
%   vol:                    volume of kspace in which to plot FS
%                           either 'uc' or 'bz' or an n x 3 matrix containing bplane_normals,
%                           normal vectors of planes bounding the volume to be plotted (Cart coords)
%   incl_point:             a point inside the volume to be plotted (Cart coords)
%   centrek:                k-point at centre of plot(?)
%   undersample_factor:     1 data point in undersample_factor points are missed along line in kspace
%                           (i.e. 1 means include all points)
%   haxes:                  handle of Matlab window to plot in
%   clipintersectedfaces:   if true, makes the edges neatly trimmed to plot region (can be slow)

%returns:
%   haxes:                  handle to Matlab window
%   hFSpatches:             Matlab patches of FS

if ~isstruct(filenameORvar)
    filename=filenameORvar;
    if ~exist([filename '.mat'])        
        filenotfound=true;
    end
    FS=load([filename '.mat']);
else
    FS=filenameORvar;
end

if isempty(bands)
    bands=1:numel(FS.Wien2k_bandnums);
end
if isempty(bandcolours)
    bandcolours=repmat([1 0 0],numel(bands),1);
end
if strcmp(lower(vol),'bz')
    bplane_normals=FS.BZfacenormals*FS.rlvs;;
elseif strcmp(lower(vol),'uc')
    bplane_normals=ucfacenormals(FS.rlvs);
    incl_point=[0.1 0.1 0.1]*FS.rlvs;
else
    bplane_normals=vol;
end
[bkparas bplane_normals]=normdir(bplane_normals);

if isempty(incl_point)
    incl_point=[0 0 0];
end

if ~isempty(undersample_factor) & undersample_factor>1
    [FS.cartX,FS.cartY,FS.cartZ,FS.cartE{bands(1)}] = reducevolume(FS.cartX,FS.cartY,FS.cartZ,FS.cartE{bands(1)},undersample_factor);
    if ~isempty(bchar)
        [bchar.x bchar.y bchar.z bchar.pchg{bands(1)}]=reducevolume(bchar.x,bchar.y,bchar.z,bchar.pchg{bands(1)},undersample_factor);
    end
    for bandnum=2:length(bands)
        FS.cartE{bands(bandnum)} = reducevolume(FS.cartE{bands(bandnum)},undersample_factor);
        if ~isempty(bchar)
            bchar.pchg{bands(bandnum)}=reducevolume(bchar.x,bchar.y,bchar.z,bchar.pchg{bands(bandnum)},undersample_factor);
        end
    end
elseif isempty(undersample_factor)
    undersample_factor=1;
end

%if axes handle isn't provided see if can find previously raised window
if isempty(haxes)
    plothandle=findobj('Tag','FermiSurfacePlot');
    if numel(plothandle)>1
        plothandle=plothandle(1);
    end
    if isempty(plothandle) | ~ishandle(plothandle)
        plothandle=figure('Tag','FermiSurfacePlot');
    end
    figure(plothandle);
    haxes=gca;
else
    axes(haxes);
end
titlestr=sprintf('Wien2k FS:%s, band %s (spin %s)',FS.pathcasename,num2str(FS.Wien2k_bandnums(bands)'),num2str(FS.spindirns(bands)'));
set(gcf,'Name',titlestr);

%find Cart coords of smallest cuboid that includes all vertices of plotting vol defined by bplane_normals
[bndryfaces bndryverts]=vol_vertices(bkparas,bplane_normals,incl_point);
exts=[min(bndryverts(:,1)) min(bndryverts(:,2)) min(bndryverts(:,3)); ...
    max(bndryverts(:,1)) max(bndryverts(:,2)) max(bndryverts(:,3))];
%round these exts out to the nearest x,y,z contained in FS.cartX,Y,Z
exts=exts./(repmat([FS.dx FS.dy FS.dz],2,1));
exts=[floor(exts(1,:)); ceil(exts(2,:))];
exts=exts.*(repmat([FS.dx FS.dy FS.dz],2,1));

if ~isempty(centrek)
    exts=exts+repmat(centrek-mean(exts,1),2,1);
end
    
[x y z]=ndgrid(exts(1,1):FS.dx*undersample_factor:exts(2,1),...
    exts(1,2):FS.dy*undersample_factor:exts(2,2),...
    exts(1,3):FS.dz*undersample_factor:exts(2,3));
wrappedcoords=wrapcoords(x,y,z,FS);    
wrappedcoords.X=permute(wrappedcoords.X,[2 1 3]); wrappedcoords.Y=permute(wrappedcoords.Y,[2 1 3]);
wrappedcoords.Z=permute(wrappedcoords.Z,[2 1 3]);

x=permute(x,[2 1 3]);
y=permute(y,[2 1 3]);
z=permute(z,[2 1 3]);

    
numbands=length(bands);
for bandnum=1:numbands            
    v=interp3(permute(FS.cartX,[2 1 3]),permute(FS.cartY,[2 1 3]),permute(FS.cartZ,[2 1 3]),...
        permute(FS.cartE{bands(bandnum)},[2 1 3]),wrappedcoords.X(1:end),wrappedcoords.Y(1:end),...
        wrappedcoords.Z(1:end),interpmethod);
    v=reshape(v,[size(wrappedcoords.X,1) size(wrappedcoords.X,2) size(wrappedcoords.X,3)]);
    if ~isempty(bchar)
        bchar.pchg{bands(bandnum)}=interp3(permute(bchar.x,[2 1 3]),permute(bchar.y,[2 1 3]),permute(bchar.z,[2 1 3]),...
            permute(bchar.pchg{bands(bandnum)},[2 1 3]),wrappedcoords.X(1:end),wrappedcoords.Y(1:end),...
            wrappedcoords.Z(1:end),interpmethod);
        bchar.pchg{bands(bandnum)}=reshape(bchar.pchg{bands(bandnum)},[size(wrappedcoords.X,1) size(wrappedcoords.X,2) size(wrappedcoords.X,3)]);
    end
    if isempty(EF)
        isoval=FS.FermiLevel(FS.spindirns(bands(bandnum)));
    else
        isoval=EF;
    end
    %compute vertices and faces of Fermi surface and plot the corresponding patch
    %graphics object 
    if ~isempty(bchar)
        pdata=isosurface(x,y,z,v,isoval,bchar.pchg{bands(bandnum)},'verbose','noshare');
    else
        pdata=isosurface(x,y,z,v,isoval,'verbose','noshare');
        pdata.facevertexcdata=repmat(bandcolours(bandnum,:),size(pdata.vertices,1),1);    
    end
    if ~isempty(bplane_normals)
        [pdata.faces,pdata.vertices,pdata.facevertexcdata]=clippatch(pdata.faces,pdata.vertices,pdata.facevertexcdata,repmat(bkparas',[1 3]).*bplane_normals,incl_point,clipintersectedfaces);
    end
    p=patch('Parent',haxes,'Faces',pdata.faces,'Vertices',pdata.vertices,'FaceVertexCData',pdata.facevertexcdata);
    %computes face normals from data and face normals of p equal to these 
    isonormals(x,y,z,v,p);
    set(p,'FaceColor','flat','edgecolor','none');
    %make scaling of plot isotropic
    if ~isempty(bchar)
        shading interp;
    end
    daspect(haxes,[1 1 1]);
    opengl software;
    %material shiny;
    %create light on right of initial view
    lights=findobj(gca,'Type','Light');
    delete(lights);
    light('Position',[1 0 0],'Style','infinite');
    light('Position',[-1 0 0],'Style','infinite');
    %set lighting algorithm
    lighting gouraud;    
    hFSpatches(bandnum)=p;
    %alpha(p,0.9);
end
xlabel('k_x'); ylabel('k_y'); zlabel('k_z');