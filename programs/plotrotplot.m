function horbpoints=plotrotplot(rplotdata,selorbitnums,BackProjectRotplot,dMdB,B,params);
   
if BackProjectRotplot
    freq=[rplotdata.freq]-B*dMdB*[rplotdata.dFdM];
else
    freq=[rplotdata.freq];    
end
horbpoints=plot([rplotdata.angle]',freq/1e3',params.LineSpec,'MarkerSize',params.MarkerSize,'MarkerEdgeColor',params.MarkerEdgeColor);
set(horbpoints,'ButtonDownFcn','dHvA_GUI(''Axes_Rotplot_ButtonDown_Callback'',gcbo,[],guidata(gcbo))');
if any(selorbitnums)
    hselorbpoints=plot([rplotdata(selorbitnums).angle]',freq(selorbitnums)/1e3',params.LineSpec,'MarkerSize',1.5*params.MarkerSize,'MarkerEdgeColor',params.MarkerEdgeColor);
    set(hselorbpoints,'ButtonDownFcn','dHvA_GUI(''Axes_Rotplot_ButtonDown_Callback'',gcbo,[],guidata(gcbo))');
end
