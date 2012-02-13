function plot_spinsplitFS(filenameORvar,bands,bplane_normals,incl_point,undersample_factor,haxes,clipintersectedfaces);

if ~isstruct(filenameORvar)
    filename=filenameORvar;
    if ~exist([filename '.mat'])
        ['Need to make ' filename '.mat using conv_bxsf']
        filenotfound=true;
    end
    bandsdata=load([filename '.mat']);
else
    bandsdata=filenameORvar;
end

if isempty(bands)
    bands=1:8;
    numpanels=4;
else
    numpanels=1;
end

if ~isempty(undersample_factor)
    [bandsdata.cartX,bandsdata.cartY,bandsdata.cartZ,bandsdata.cartE{bands(1)}] = reducevolume(bandsdata.cartX,bandsdata.cartY,bandsdata.cartZ,bandsdata.cartE{bands(1)},undersample_factor);
    for bandnum=2:length(bands)
        bandsdata.cartE{bands(bandnum)} = reducevolume(bandsdata.cartE{bands(bandnum)},undersample_factor);
    end
else
    undersample_factor=1;
end

%if axes handle isn't provided see if can find previously raised window
if isempty(haxes)
    plothandle=findobj('Tag','FermiSurfacePlot');
    if isempty(plothandle) | ~ishandle(plothandle)
        plothandle=figure('Tag','FermiSurfacePlot');
    end
    figure(plothandle);
    haxes=gca;
else
    axes(haxes);
end
titlestr=sprintf('Wien2k FS:%s, band %2.0f (spin %1.0f)',bandsdata.pathcasename,bandsdata.Wien2k_bandnums(bands),bandsdata.spindirns(bands));
set(gcf,'Name',titlestr);
set(gca,'Position',[0 0 1 1]);

if ~isempty(bplane_normals)
    %make huge grid which we will chop bits off
    [x y z]=ndgrid(-max(bandsdata.cuboidextents(:,1))/2:bandsdata.dx*undersample_factor:max(bandsdata.cuboidextents(:,1))/2,...
        -max(bandsdata.cuboidextents(:,2))/2:bandsdata.dy*undersample_factor:max(bandsdata.cuboidextents(:,2))/2,...
        -max(bandsdata.cuboidextents(:,3))/2:bandsdata.dz*undersample_factor:max(bandsdata.cuboidextents(:,3))/2);
    wrappedcoords=wrapcoords(x,y,z,bandsdata);    
else
    x=bandsdata.cartX; y=bandsdata.cartY; z=bandsdata.cartZ;
end
x=permute(x,[2 1 3]);
y=permute(y,[2 1 3]);
z=permute(z,[2 1 3]);

bandcolor=[1 0.2 0.2; 0.2 0.2 1];
numbands=length(bands);

for panelnum=1:numpanels
%haxes=subplot(2,2,panelnum);

for bandnum=panelnum*2-1:panelnum*2
    if ~isempty(bplane_normals)      
        v=interp3(permute(bandsdata.cartX,[2 1 3]),permute(bandsdata.cartY,[2 1 3]),permute(bandsdata.cartZ,[2 1 3]),...
            permute(bandsdata.cartE{bands(bandnum)},[2 1 3]),reshape(wrappedcoords.X,[],1),reshape(wrappedcoords.Y,[],1),...
            reshape(wrappedcoords.Z,[],1));
        v=reshape(v,[size(wrappedcoords.X,2) size(wrappedcoords.X,1) size(wrappedcoords.X,3)]);
        %v(~incl)=10*bandsdata.FermiLevel(bandsdata.spindirns(bands(bandnum)));
    else
        v=permute(bandsdata.cartE{bands(bandnum)},[2 1 3]);
    end    
    %compute vertices and faces of Fermi surface and plot the corresponding patch
    %graphics object 
    [faces,verts]=isosurface(x,y,z,v,bandsdata.FermiLevel(bandsdata.spindirns(bands(bandnum))),'verbose','noshare');
    [cfaces,cverts,fcolors]=clippatch(faces,verts,[bplane_normals; [0.001-0.002*(mod(bandnum-1,2)==1) 0 0]],[0 0 0],clipintersectedfaces);
    fcolors=repmat(bandcolor(bandnum,:),size(cfaces,1),1);
    p=patch('Parent',haxes,'Faces',cfaces,'Vertices',cverts,'FaceVertexCData',fcolors);
    %computes face normals from data and face normals of p equal to these 
    isonormals(x,y,z,v,p);
    set(p,'FaceColor','flat','edgecolor','none');
    %make scaling of plot isotropic
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
[BZfaces BZvertices]=vol_vertices(bandsdata.BZfacenormals*bandsdata.rec_latt_vecs',[0 0 0]);      
C=[0 1 0];
h.BZpatch=patch('Faces',BZfaces,'Vertices',BZvertices,'FaceColor',C,'FaceAlpha',0.1,'EdgeColor',[0 0 0]);
xlabel('k_x'); ylabel('k_y'); zlabel('k_z');
symmpoints=[0 0 0; 1.01 0 0 ; 0.5 0.5 0.5; 1 0 -0.5 ; 1 0.25 -0.25]*abs(bandsdata.rec_latt_vecs(1,1));
symmlabels={'\Gamma','X','L','W','U'};
%lbloffs=[6 6; 0 10; 10 14; 6 -6; 6 -6];
lbloffs=[0.08 0 0; 0 0 0.12; 0.1 0 0.14; 0.08 0 0; 0.08 0 0];
hold on;
view(35,30);

for n=1:numel(symmlabels)
    plot3(symmpoints(n,1),symmpoints(n,2),symmpoints(n,3),'ok','MarkerFaceColor','black','MarkerSize',3,'Clipping','off');
    htxt=text(symmpoints(n,1)+lbloffs(n,1),symmpoints(n,2)+lbloffs(n,2),symmpoints(n,3)+lbloffs(n,3),symmlabels{n},'Clipping','off','FontWeight','bold','FontSize',12,'FontName','Times New Roman');    
    %set(htxt,'Units','pixels');
    %currpos=get(htxt,'Position');
    %set(htxt,'Position',[currpos(1)+lbloffs(n,1) currpos(2)+lbloffs(n,2) 0]);
    %set(htxt,'Units','data');
end
set(gca,'Visible','off');
xlim([-1 1.6]);
end

