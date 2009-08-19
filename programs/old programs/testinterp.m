function testinterp

a=ones(2,2,2);
a(:,:,2)=2;

%a=interp3(a,2);

plothandle=findobj('Tag','testPlot');
if ishandle(plothandle)
    %do nothing
else
    plothandle=figure('Tag','testPlot');    
end
plothandle;

[gridX,gridY,gridZ]=ndgrid(-0.5:0.5,-0.5:0.5,-0.5:0.5);

isosurface(permute(gridX,[2 1 3]),permute(gridY,[2 1 3]),permute(gridZ,[2 1 3]),a,1.5);

tgridX=(gridX/sqrt(2)+gridZ/sqrt(2));
tgridY=gridY;
tgridZ=(-gridX/sqrt(2)+gridZ/sqrt(2));
[tgridX tgridY tgridZ]=rotate3grid(gridX,gridY,gridZ,[0 1 1]);
%auntilt=interp3(gridX,gridY,gridZ,a,gridX,gridY,gridZ);
atilt=interp3(permute(gridX,[2 1 3]),permute(gridY,[2 1 3]),permute(gridZ,[2 1 3]),a,reshape(tgridX/2,[],1),reshape(tgridY/2,[],1),reshape(tgridZ/2,[],1));
atilt=reshape(atilt,2,2,2);

xlabel('X'); ylabel('Y'); zlabel('Z');
isosurface(permute(gridX,[2 1 3]),permute(gridY,[2 1 3]),permute(gridZ,[2 1 3]),atilt,1.5);

%[tgridX,tgridY,tgridZ]=rotate3grid(
