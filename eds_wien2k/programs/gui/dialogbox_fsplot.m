function dialogbox_FSplot(hcallingfig);

figure('Tag','DB_FSplot');
cpos=get(gcf,'Position');
set(gcf,'Units','pixels','Position',[cpos(1) cpos(2) 300 300],'Name','Fermi Surface Plot Settings');
uicontrol('Style','text','String','Bounding planes','Position);


%vars = evalin('base','who');