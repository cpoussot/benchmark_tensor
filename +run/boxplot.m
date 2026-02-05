function h = boxplot(x, y, w, col)
y = sort(y,'ascend');
N = length(y);
[miny, maxy, meany, medy, iQ1, iQ3] = run.moustache_data(y);
%
if ~isnan(meany) 
    %
    hold on 
    patch('Vertices',[x-w/2, y(iQ1); x-w/2, y(iQ3); x+w/2 , y(iQ3);x+w/2, y(iQ1)],...
          'Faces',[1,2,3,4],'EdgeColor',col, 'FaceColor',col,'FaceAlpha',0.1, ...
          'LineWidth',2);
    %
    h=plot(x, meany, 'x', 'color',col);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h=plot([x-w/2, x+w/2], medy*[1,1],'color',col);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h=plot(x*[1,1],[miny, y(iQ1)],'color',col);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h=plot(x*[1,1], [y(iQ3), maxy],'color',col);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h=plot(x+[-w/4,w/4], miny*[1,1],'color',col);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h= plot(x+[-w/4,w/4], maxy*[1,1],'color',col);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %text(x+.25, maxy,sprintf('~~~%d',N), 'Rotation',0);
    %text(x, maxy,sprintf('~~%d',N), 'Rotation',90);
else
    hh = gca;
    %plot(x*[1,1],[min(hh.YLim(:)), max(hh.YLim(:))],'color','r');
    %plot(x+[-w/4,w/4],[min(hh.YLim(:)), max(hh.YLim(:))],'color','r');
    plot(x, max(hh.YLim(:)), '*', 'color',col);
    text(x, min(hh.YLim(:)*10),'NOT CONVERGED','Rotation',90,'Color',col,'FontWeight','bold');
    % %
    % YMAX = max(hh.YLim);
    % YMIN = min(hh.YLim);
    % patch('Vertices',[x-w/2, YMIN; x-w/2, YMAX; x+w/2 , YMAX;x+w/2, YMIN],...
    %     'Faces',[1,2,3,4],'EdgeColor',col, 'FaceColor',col,'FaceAlpha',0.1, ...
    %     'LineWidth',2);

end